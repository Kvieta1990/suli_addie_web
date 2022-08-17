from parser import *
from resonance_calc import *
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    user_input = "Sr1. Ti (16O)3."

    chem_dict = formula_parser(user_input)

    print("Parsed chemical formula -> ", chem_dict)
    print("FZ & Self-scattering level & Partials -> ", fz_ss_calculator(chem_dict))
    print("Absorption factor -> ", absorption_calculator(chem_dict))

    formula = ""
    for key, item in chem_dict.items():
        formula += f"{key} {item} "
    rhos1 = 1.0
    radius = 0.3
    sample_title = "test"
    atconc, data = grab_data(formula, "./nistall.dat")
    
    lamda, oneminusahkl = resonancecalc(formula,atconc,data,rhos1,radius,sample_title)
    
    plt.plot(lamda, oneminusahkl, '-',lamda,1. - np.zeros(601) - 1./np.exp(1),'--')
    plt.xscale('log')
    plt.show()
