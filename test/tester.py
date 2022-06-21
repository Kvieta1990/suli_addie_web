from parser import *

if __name__ == "__main__":

    user_input = "Sr1. Ti O3."

    chem_dict = formula_parser(user_input)

    print("Parsed chemical formula -> ", chem_dict)
    print("FZ & Self-scattering level & Partials -> ", fz_ss_calculator(chem_dict))
    print("Absorption factor -> ", absorption_calculator(chem_dict))
    print("Packing fraction values -> ", pf_calculator(chem_dict))
