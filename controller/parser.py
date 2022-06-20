import json
import math

avogadro_const = 6.02E+23

f = open('elements.json', "r")
data = json.load(f)
f.close()


def validator(input_str):
    '''Validator for input chemical formula

    :param input_str: Input chemical formula string
    :return: Boolean for whether the input meets our requirement
    "rtype:" boolean
    '''

    capitals = sum(1 for c in input_str if c.isupper())

    count = 0

    for i in range(0, len(input_str)):
        if input_str[i] == " ":
            count += 1

    if capitals - count == 1:
        return True

    else:
        return False


def formula_parser(input_str: str):
    '''Parser for input chemical formula string

    :param input_str: Input chemical formula string
    :return: Dictionary for the processed formula
    :rtype: dict
    '''

    formula = input_str

    if validator(formula):
        parts = formula.split()

        y = 0
        comp_dict = {}

        while True:
            ele_name = ""
            i_tmp = 0
            for x in parts[y]:
                if x.isdigit() or x == ".":
                    break
                else:
                    ele_name += x
                    i_tmp += 1
            if i_tmp == len(parts[y]):
                z = 1.0
            else:
                try:
                    z = float(parts[y][i_tmp:])
                except ValueError:
                    return
            comp_dict[ele_name] = str(z)

            y += 1
            if y == len(parts):
                break

        return comp_dict
    else:
        return


def fz_ss_calculator(input_dict: dict):
    '''Calculator for Fzber-Ziman coefficient and self-scattering level

    :param input_dict: The processed chemical formula
    :return: Dictionary for the calculated factors
    :rtype: dict
    '''

    c = []
    d = []
    b = []
    b_c = []
    bcval = []
    bbar = []
    fz = []

    fzss_dict = {}

    # FZ coefficient and self-scattering level

    divisor = 0
    for count, key in enumerate(input_dict.keys()):
        d.append((data[key]['ss']))
        b.append(math.sqrt(d[count] / (4 * math.pi)))
        bcval.append(data[key]['bc'])
        c.append(input_dict[key])
        divisor += float(input_dict[key])

    bbar = [i / 10. for i in bcval]

    c = [float(item) / divisor for item in c]

    fz += [x * i for i, x in zip(bbar, c)]
    fz_val = (sum(fz)) ** 2

    ssl_val = 0.0
    for i in range(len(c)):
        b_c.append((c[i] * (b[i] ** 2)))
        ssl_val += b_c[i]

    fzss_dict['SelfScattering'] = ssl_val
    fzss_dict['FaberZiman'] = fz_val

    # Coefficients for partials

    fzss_dict['partials'] = {}

    all_ele = [key for key in input_dict.keys()]

    for count, key in enumerate(all_ele):
        for count1 in range(count, len(all_ele)):
            partial_name = key + "-" + all_ele[count1]
            if key != all_ele[count1]:
                partial_val = c[count] * c[count1]
                partial_val *= (bbar[count] * bbar[count1] * 2)
            else:
                partial_val = c[count] * c[count1] * bbar[count] * bbar[count1]
            fzss_dict['partials'][partial_name] = partial_val

    return fzss_dict
    
def absorption_calculator(input_dict: dict):

    sa_val = []
    ele_val = []
    abs_dict = {}
    abs_val = 0
    for i in input_dict:
        ele_val.append(input_dict[i])
        sa_val.append(data[i]['sa'])
    ele_int = [int(float(i)) for i in ele_val]
    for x, y in zip(ele_int, sa_val):
        abs_val += x * y
    abs_dict = abs_val

    return abs_dict

def pf_calculator(input_dict: dict):

    '''MassSampleCan = 40.762 g
       MassEmptyCan = 38.674 g
       CanInnerRadius = 0.3 cm
       SampleHeight = 4.6 cm
       FullNumDensity = 0.0905
       LevelAtHighQ = 0.452127
       SelfScattering = 0.370767'''
    pf_input = {}
    pf_dict = {}
    pf_input['Input'] = {}
    pf_dict['Output'] = {"NumDensity": "", "FullDensity": "",
                         "InitialPackingFrac": "", "SuggestedPackingFrac": ""}

    ele_val = []
    ele_wt = []

    key_list = ["MassSampleCan", "MassEmptyCan", "CanInnerRadius",
                "SampleHeight", "FullNumDensity", "LevelAtHighQ",
                "SelfScattering"]

    for i in key_list:
        pf_input['Input'][i] = input()
    mass = float(pf_input['Input']["MassSampleCan"])
    mass -= float(pf_input['Input']["MassEmptyCan"])
    volume = math.pi * (float((pf_input['Input']["CanInnerRadius"])) ** 2)
    volume *= float(pf_input['Input']["SampleHeight"])
    density = mass/volume

    for i in input_dict:
        ele_val.append(input_dict[i])
        ele_wt.append(data[i]['Atwt'])
    ele_int = [(float(i)) for i in ele_val]

    e11 = 0
    for x, y in zip(ele_int, ele_wt):
        e11 += x * y
    full_density = (float(pf_input['Input']["FullNumDensity"]) * e11)
    full_density /= (sum(ele_int) * avogadro_const * 1E-24)
    number_density = density / e11 * sum(ele_int) * avogadro_const * 1E-24
    init_pf = number_density / float(pf_input['Input']["FullNumDensity"])
    sugg_pf = init_pf * float(pf_input['Input']["LevelAtHighQ"])
    sugg_pf /= float(pf_input['Input']["SelfScattering"])
    dict_keys = [number_density, full_density, init_pf, sugg_pf]

    for i, j in zip(pf_dict['Output'], dict_keys):
        pf_dict['Output'][i] = j

    return pf_input, pf_dict
    
