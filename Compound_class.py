import re
import molecule
import numpy as np


def del_repetition(text: str):
    parsed = re.findall(r"([A-Z][a-z]*)(\d*)|(\()|(\))(\d*)", text)
    structure_dict = {}
    output = ''
    for element_details in parsed:
        element = element_details[0]
        if element not in structure_dict:
            structure_dict[element] = 0
            element_count = sum([int(x)
                                 for x in element_details[1:] if x != ""])
            if element_count > 0:
                structure_dict[element] += element_count
            else:
                structure_dict[element] = 1
        else:
            if element_details[1] != "":
                structure_dict[element] += int(element_details[1])
            else:
                structure_dict[element] += 1

    for element in structure_dict.items():
        if element[1] != 1:
            output += str(element[0]) + str(element[1])
        else:
            output += str(element[0])
    return(output)


# classe Compounds

class Compound:
    def __init__(self, Molecule, Adduct, Charge, Label):
        self.mol= Molecule
        self.add= Adduct
        self.charge= int(Charge)
        self.compound = (del_repetition(Molecule + Adduct))
        self.label = Label
    
    @property #Isotopic Mass Caculator
    def mass(self, test=True):
            molecola = molecule.Molecule(self.mol)
            mass = (molecola.isotipic_molecular_weight)
            return(mass)
    
    @property #Mass over Charge calculator
    def MoC(self):
        molecola = molecule.Molecule(self.compound)
        moz = round(((molecola.isotipic_molecular_weight +(self.charge * (-0.00054387)))) / abs(self.charge), 5)
        return(moz)