'''
Funzione per il calcolo del pattern isotopico inserendo la formula bruta della molecole.
Utilizzo la libreria pyisopatch per la lettura della molecola e per il calcolo delle masse.
'''
import itertools as it
from typing import Pattern
import numpy as np
import math

from numpy import NaN, testing
import molecule
import pandas as pd


# Funzione per raccogliere i dati (masse, ratios, num. atomi) relativi agli elementi
# presenti nella molecola iniziale.

def Generate_molecule(molecola: str):
    mol = molecule.Molecule(molecola)
    weights, ratios, atoms = [], [], []

    for elem in mol._as_elements:
        weights.append(elem.isotopic_weight)
        ratios.append(elem.isotopic_ratios)
        atoms.append(elem.count)

    w = np.array(weights, dtype=object)
    r = np.array(ratios, dtype=object)
    a = np.array(atoms)
    return(w, r, a)


# Funzione con la quale calcolo il coeff. multinomiale.

def Multinomial_coeff(combination):
    den = 1
    for i in combination:
        fact = math.factorial(i)
        den = den * fact
    return(math.factorial(sum(combination)) / den)



# Funzione per generare l'output: dopo aver generato l'output si fa un arrotondamento e
# poi si fitrano solo gli isotopi con frequenza maggiore dello 0.01%.


####################################################################################


def generate_output(raw_weights, raw_ratios, charge: int):
    raw_weights = np.reshape(raw_weights, (np.size(raw_weights), 1))
    raw_ratios = np.reshape(raw_ratios, (np.size(raw_ratios), 1))

    out = np.hstack((raw_weights, raw_ratios))
    if charge != 0:
        out[:, 0] = out[:, 0]-(charge*0.00054387)
        out[:, 0] /= abs(charge)
    
    # Isotopic Patter Normalization 
    out[:, 1] /= out[:, 1].max()
    # Sorting on m/z small -> large
    out = out[out[:, 1].argsort()]
    out = out[::-1]
    #elimino le intensita' minori di 0.00001 --> 10^-5 volte piu piccole rispetto al principale
    out = out[np.where(out[:,1] >= 0.00001)]
    return(out)
    

# Funzione per effettuare un merge dei picchi del pattern isotopico se la loro differenza
# in massa è inferiore al valore "soglia". Se la condizione è soddisfatta viene effettuata,
# una media pesata per il valore di m/z e vengono sommate le intensità.

# Merge del pattern isotopico, ma utilizzando la differenza in ppm

def PatternFilter(input_array, merge_threshold, ratio_threshold):
    output_array = input_array
    for i in range(len(input_array)-1):
        if abs(input_array[i, 0] - input_array[i+1, 0]) < merge_threshold:
            new_I = input_array[i, 1] + input_array[i+1, 1]
            new_mz = (input_array[i, 0]*input_array[i, 1] + input_array[i+1, 0]*input_array[i+1, 1])/ new_I

            output_array[i+1, 0], output_array[i+1, 1] = new_mz, new_I
            output_array[i, 0], output_array[i, 1] = NaN, NaN

    output_array = output_array[~np.isnan(output_array).all(1)]
    
    ### ELIMINO I PICCHI SOTTO SOGLIA DI INTENSITA
    output_array = output_array[np.where(output_array[:, 1] >= ratio_threshold)]
    ### RINORMALIZZO DOPO AVER UNITO I PICCHI
    #output_array[:, 1] /= output_array[:, 1].max()
    
    #output_array = output_array[output_array[:, 1].argsort()]
    #output_array = output_array[::-1]
    return output_array


# Funzione con la quale calcolo la probailità e la massa della molecola.

def calcolatore(Isotopic_ratio, Isotopic_weight, combination):
    probability = 1
    mass = 0
    for idx1, idx2 in zip(Isotopic_ratio, combination):
        probability = (idx1**idx2) * probability

    for idx1, idx2 in zip(Isotopic_weight, combination):
        mass = (idx1 * idx2) + mass

    return(mass, probability * Multinomial_coeff(combination))

#################################################
# FUNZIONE PER IL CALCOLO DEL PATTERN ISOTOPICO #
#################################################


def Patter_Calculator(molecola: str, charge: int, merge_threshold, ratio_threshold):
    w, r, a = Generate_molecule(molecola)
    ratio, weight = [], []

    for index, element in enumerate(a):
        M = len(r[index])
        atom = [i for i in range(element + 1)]
        comb = it.product(atom, repeat=M)

        Pattern = []
        i = 1
        for val in list(comb):
            if sum(val) == element:
                c = calcolatore(r[index], w[index], val)
                Pattern.append(c)
                i += 1

        ratio.append(np.array(Pattern)[:,1])
        weight.append(np.array(Pattern)[:, 0])

    raw_ratios = ratio[0]
    raw_weights = weight[0]

    for index2 in range(0, len(weight[:-1])):
        raw_weights = np.add.outer(raw_weights, weight[index2+1])
        raw_ratios = np.outer(raw_ratios, ratio[index2+1])

    New_Output = generate_output(raw_weights, raw_ratios, charge)
    Filtered = PatternFilter(New_Output, merge_threshold, ratio_threshold)  

    return Filtered
    
######################################
#### MAIN PROGRAM - Just for DEBUG ####


#Single molecule
#c =Patter_Calculator(molecola='FeC5N5H40Hp', charge=-1, merge_threshold=0.0005, ratio_threshold=0.001)


