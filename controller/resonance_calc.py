def grab_data(sample_formula, file_in):
    # Method for grabbing the isotopes information.
    import numpy as np

    f=open(file_in,'r')
    line1=f.readlines()
    whatsit=line1[0].split('<th>')
    data={}
    allkeys=[]

    for name in whatsit:
        name=name.strip()
        data[name] = []
        allkeys.append(name)

    for linex in line1[1:]:
        linexs=linex.split('<td>')
        for i,val in enumerate(allkeys):
            try:
                data[val].append(linexs[i].split()[0])
            except IndexError:
                print(linexs)

    f.close()

    Isotopes = np.array(data["Isotope"])

    (where_iso,conc_iso)=interpret_samplename(sample_formula,Isotopes)
    atconc=np.array(conc_iso)

    return atconc, data


def is_number(s):

    a=s.find('(')
    if a >0:
        s=s[0:a]

    s=s.replace("<i>i</i>","j")

    try:
        complex(s)
        return complex(s)
    except ValueError:
        return float('nan')


def interpret_samplename(sample,isotopelist):
    from numpy import array,arange,where

    samplesplit=sample.split()
    iso_in_sample=samplesplit[::2]
    conc_iso=samplesplit[1::2]
    conc_iso=array([is_number(val) for val in conc_iso])
    conc_iso=conc_iso.real

    where_iso=[]
    r=arange(len(isotopelist))
    for iso in iso_in_sample:
        if iso in isotopelist:
            where_iso.append(where(iso == isotopelist)[0][0])
        elif iso not in isotopelist:
            print(iso, "is a strange isotope. Please check")

    return where_iso,conc_iso


def rstd(name):
    """ python version of the IDL function rstd, opens a file,
    ignores all lines starting with # and reads as n column ascii
    """
    f=open(name,'r')
    lines=f.readlines()
    a=[]
    for line in lines:
        if line[0:1] != '#':
            zeil=line.split()
            a.append([float(element) for element in zeil])
    q=[b[0] for b in a]
    i=[b[1:] for b in a]
    f.close()
    return q,i


def resonancecalc(formula,atconc,data,rhos1,radius,sample_title):
    import pdb
    from numpy import array,arange,where, linspace,zeros,exp
    from numpy import interp
    import matplotlib.pyplot as plt


    isotopes=array(data["Isotope"])
    sig18=array(data["Abs_xs"])
    conc=array(data["conc"])
    lamda=linspace(0.1,6.1,601)
    xabstot=zeros(601)

    (where_iso,conc_iso)=interpret_samplename(formula,isotopes)
    si=zeros((len(where_iso),601))

    for counttype,val in enumerate(where_iso):

        iselement=conc[val] == '---'
        if iselement:
            natcomp=[]
            inisotopes=[]
            single_iso_element=(conc[val+1] == '---' or conc[val+1] == '100')
            if not single_iso_element:
                val1=val
                iniso=[]
                inisoconc=[]
                iniso_absxs=[]
                while conc[val1+1] != '---' and conc[val1+1] != '100':
                    iniso.append(isotopes[val1+1])
                    inisoconc.append(conc[val1+1])
                    iniso_absxs.append(sig18[val1+1].split('(')[0])
                    val1=val1+1
        if not iselement or single_iso_element:
            iniso=[isotopes[val]]
            inisoconc=['100']
            iniso_absxs=[sig18[val].split('(')[0]]
        for counter,thisiso in enumerate(iniso):
            nodata=False
            try:
                f1=open('./abslam/'+thisiso+'.dat','r')
            except IOError:
                nodata=True

            if not nodata:
                f1.close()
                (e,sigma)=rstd('./abslam/'+
                               thisiso+'.dat')
                l1 = array([(8.18145447496e-22/e1)**.5*1e10 for e1 in e])
                s1 = array([val[0] for val in sigma])
                f=interp(lamda,l1[::-1],s1[::-1])

                si[counttype,:]+=f*float(inisoconc[counter])/100.
                si[counttype,:]*=float(iniso_absxs[counter])/f[170]
                print(thisiso,f[170],iniso_absxs[counter],float(inisoconc[counter])/100.)
            elif True:
                si[counttype,:]+=float(iniso_absxs[counter])*lamda/1.8*float(inisoconc[counter])/100.

        xabstot+=si[counttype,:]*float(atconc[counttype])

    mur=xabstot*rhos1*radius
    ahkl=mur*0
    ahkl[(mur<10)]=exp(-1.7133*mur[(mur<10)]+.0927*mur[(mur<10)]**2)

    return lamda, 1-ahkl
