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

    return fzss_dict
