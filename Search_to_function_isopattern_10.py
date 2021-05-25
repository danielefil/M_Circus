from mmap import ACCESS_COPY
from Compound_class import Compound
from pathlib import Path
from scipy import optimize
from collections import Counter
from IsoPattern3 import Patter_Calculator
import pandas as pd
import numpy as np
import itertools as it
import re
from progressbar import ProgressBar

## FILTRO ##

# GAUSSIAN FUNCTION (FOR NOISE ESTIMATION USING SCIPY.OPTIMIZE.CURVE_FIT)
def gaussian(x: np.array, mu: float, sd: float, a: float):
    return(a * np.exp(-0.5 * (((x - mu)**2) / (np.abs(sd)**2))))

# NOISE ESTIMATOR
def noiseEstimator(fit_parameters: list):
    noiseMean = fit_parameters[0]
    noiseSD = np.abs(fit_parameters[1])
    noiseEstimate = noiseMean + 3 * noiseSD
    return(noiseEstimate)

# GAUSSIAN FUNCTION (FOR NOISE ESTIMATION USING SCIPY.OPTIMIZE)
def gaussian_diff(params, x, y):
    mu=params[0]
    wid=params[1]
    amp=params[2]
    gauss_est = amp * np.exp(-0.5 * (((x - mu)**2) / (np.abs(wid)**2)))
    return(np.sum((y-gauss_est)**2))


def Apply_filterZ(spectra_df, parameters: list):
    _maxIntenseLimit = spectra_df.quantile(parameters[1])
    _filter = spectra_df['Intensity'] < _maxIntenseLimit[1]
    
    filtered = spectra_df[_filter].to_numpy() 

    counts, bins = np.histogram(filtered[:, 1], bins=parameters[2])
    c_bins = bins + ((bins[1]-bins[0])/2)

    # initize gaussian fit parameters
    p0 = [c_bins[np.argmax(counts)], 2000, np.max(counts)]

    try:
        fit_parameters, cov = optimize.curve_fit(gaussian, c_bins[:-1], counts, absolute_sigma=True, p0=p0, maxfev=100000)
        noise = noiseEstimator(fit_parameters)
    except RuntimeError:
        noise = 1

    _noisefilter = spectra_df['Intensity'] >= noise
    spectra_df_filtered = spectra_df[_noisefilter]

    return(spectra_df_filtered)

## ADDUCT GENERATOR

def adduct_generator(adducts: list, adduct_label: list, charge: int):
    _adductsList = []
    np.column_stack((adducts, adduct_label))
    for i in range(1, (abs(charge) + 1)):
        elem_comb = list(it.combinations_with_replacement(adducts, i))
        label_comb = list(it.combinations_with_replacement(adduct_label, i))
    # Subrutine per il calcolo del rapporto m/z per le combinazioni di addotti
        for i_Lc, i_Ec in zip(label_comb, elem_comb):
            if charge > 0:
                # [Formula, carica, label]
                _adductsList.append([i, ''.join(i_Ec), ''.join(i_Lc)])
            else:
                # [Formula, carica, label]
                _adductsList.append([-i, ''.join(i_Ec), ''.join(i_Lc)])
    # np array con gli addotti - [Charge, Elements, Labels]
    _adducts = np.array(_adductsList)
    return(_adducts)



#### RICERCA PICCHI ##########

def diff(spectra, Cmp_MoC, thold, mode):
    if mode == 'ppm':
        _diff = np.absolute((spectra - Cmp_MoC) * 1000000 / Cmp_MoC)
        min_diff = np.amin(_diff)
    elif mode == 'dalton':
        _diff = np.absolute(spectra - Cmp_MoC)
        min_diff = np.amin(_diff)

    if min_diff < thold:
        index_min = np.where(_diff == np.amin(_diff))
        return(index_min[0][0], min_diff)
    else:
        return(None)

#### RICERCA PATTERN ISOTOPICO ##########

def find_pattern(theo_pattern, spectra: np.array, _diff: int, search_mode: str):
    pattern_tbl = np.empty((0, 4), float)
    exp_p = 0
    for iso_molecule in theo_pattern:
        if search_mode == 'ppm':
            find_iso = diff(spectra[:, 0], iso_molecule[0], _diff, 'ppm')
        else:
            find_iso = diff(spectra[:, 0], iso_molecule[0], _diff, 'dalton')
        if find_iso != None:
            # lista output: m/z sperimentale, int. sperimentale, m/z teorico, ubt. teorico
            exp_p += 1
            ls = [spectra[find_iso[0], 0], spectra[find_iso[0], 1],
                  iso_molecule[0], iso_molecule[1]]
            pattern_tbl = np.vstack((pattern_tbl, np.array(ls)))
            ### Cancello dallo spettro il picco che ho trovato cosi non lo considero piu nel corso del calcolo del pattern isotopico
            spectra = np.delete(spectra, find_iso[0], axis=0)
        else:
            ls = [0, 0, iso_molecule[0], iso_molecule[1]]
            pattern_tbl = np.vstack((pattern_tbl, np.array(ls)))
    
    # Normalizzo il pattern isotopico solo se trovo dei picchi corrispondenti
    #if pattern_tbl[:, 1] != 0:
    #np.seterr(all='warn')
    pattern_tbl[:, 1] /= pattern_tbl[0, 1]
    #input('Press any key to continue...')


    findedpeaks = str(exp_p) + ' of ' + str(len(theo_pattern))
    accord = accordance(pattern_tbl)
    return(accord, pattern_tbl, findedpeaks)


# ACCORDANCE CALCULATOR
def accordance(pattern_tbl: np.array):
    if pattern_tbl.size == 0:
        _accordance = 0
    else:
        _accordance = sum(
            pattern_tbl[:, 3] * abs(pattern_tbl[:, 1] - pattern_tbl[:, 3]))
        _max = sum(pattern_tbl[:, 3]**2)
        _accordance = 100 - 100 * _accordance / (_max-1.0000001)
        if _accordance < 0:
            _accordance = 0

    return(round(_accordance, 1))


# ACCORDANCE CALCULATOR
def pattern_score(pattern_tbl: np.array):
    if pattern_tbl.size == 0:
        _accordance = 0
    else:
        _accordance = 1
        _pattern_tbl = np.copy(pattern_tbl[1:,:])
        _pattern_tbl[:, 1] = _pattern_tbl[:, 1]*(-1)
        mz_list = _pattern_tbl[:, 0].tolist() + _pattern_tbl[:, 2].tolist()
        i_list = _pattern_tbl[:, 1].tolist() + _pattern_tbl[:, 3].tolist()

        aray = np.column_stack((mz_list,i_list))
        aray = aray[~np.all(aray == 0, axis=1)]
        aray = aray[aray[:, 0].argsort()]

        _sum = aray[0,0]
        mz_delta = 0.001
        while _sum < aray[-1, 0]:
            _min = _sum
            _max = _min + mz_delta
            sub = aray[np.where((aray[:, 0] >= _min) & (aray[:, 0] <= _max))]
            _sum += mz_delta
            if sub.size != 0:
                _res = abs(sub[:,1].sum())
                _accordance *= 1-_res

    return(round(_accordance * 100, 0))

# SAVE PATTERN
def SaveOutput(output_list, input_filepath, Diff):
    colnames = ['PM', 'm/z(teorico)', 'm/z(trovato)', Diff, 'Intens. Ass.', 'Fromula Bruta', 'Addotto', 'Accordo', 'Score', 'Picchi']
    output_df = pd.DataFrame.from_records(output_list, columns=colnames)
    
    filename = Path(input_filepath).stem
    filepath = Path(input_filepath).parent
    outputstring = str(filepath)+'/RES/RESULTS_'+ filename +'.csv'

    filtro = output_df['Accordo'] > 0.1
    output_df[filtro].to_csv(outputstring, index=False)


############################################################################################

def search_peak(spectra_path: str, database_path: str, adduct_list: list, charge,  search_property: list, label_list: list, Filtering: list):
    # Leggo il file che contiene la lista di composti da cercare ## DATABASE ##
    df_db = pd.read_csv(database_path, sep='\t', dtype={'Formula': str, 'Mass': float})
    database = df_db.to_numpy()

    # Genero gli addotti
    adducts = adduct_generator(adduct_list, label_list, charge)
    CmpsToFind = []
    Iso_dict = {}


    # Genero i composti da cercare (Database + Adducts)
    for entries in database:
        for adduct in adducts:
            # adduct = [Charge, Elements, Labels]
            _Comp = Compound(entries[0], adduct[1], adduct[0], adduct[2])
            CmpsToFind.append(_Comp)
    


    # CmpsToFind contine tutte le info sui composti che devo cercare


    for spectrum_path in spectra_path:
        df_sp = pd.read_csv(spectrum_path, dtype={'Mass [m/z]': float, 'Intensity': float})
        
        
        # Filtraggio delle intensita basata sul Arthur's method
        if Filtering[0]:
            df_sp = Apply_filterZ(df_sp, Filtering)
        
        spectra = np.round(df_sp.to_numpy(), 5)
    
        # Genero liste per gli output
        output = []

        
        for Comp in CmpsToFind:
            if search_property[0] == 'ppm':
                # FIND 1 - Cerco i valori per i quali il ppm è minore di un valore
                find = diff(spectra[:, 0], Comp.MoC, search_property[1], search_property[0])
                if find != None:
                    if Comp.compound not in Iso_dict:
                        Iso_dict[Comp.compound] = Patter_Calculator(Comp.compound, Comp.charge, .0005, .001)
                    pattern_t = Iso_dict[Comp.compound]
                    # Find the isotopic pattern
                    accordance, pattern_tbl, findedrate = find_pattern(pattern_t, spectra, search_property[1], search_property[0])                    
                    score = pattern_score(pattern_tbl)
                    #                 PM     m/z(teorico)    m/z(trovato)            Diff         intens. Ass.        Formula Bruta  Addotto     Accordo   Score    Picchi
                    output.append([Comp.mass, Comp.MoC, float(spectra[find[0], 0]), find[1], float(spectra[find[0], 1]), Comp.mol, Comp.label, accordance, score, findedrate])
                else:
                    pass
                    
            elif search_property[0] == 'dalton':
                # FIND 2 - Cerco i valori per i quali la massa è compresa tra un intervallo
                find = diff(spectra[:, 0], Comp.MoC, search_property[1], search_property[0])
                if find != None:
                    if Comp.compound not in Iso_dict:
                        Iso_dict[Comp.compound] = Patter_Calculator(Comp.compound, Comp.charge, 'dalton', .001, .001)
                    pattern_t = Iso_dict[Comp.compound]
                    # Find the isotopic pattern
                    accordance, pattern_tbl, findedrate = find_pattern(pattern_t, spectra, search_property[1], search_property[0])
                    score = pattern_score(pattern_tbl)
                    #                 PM     m/z(teorico)    m/z(trovato)            Diff         intens. Ass.        Formula Bruta  Addotto     Accordo   Score    Picchi
                    output.append([Comp.mass, Comp.MoC, float(spectra[find[0], 0]), find[1], float(spectra[find[0], 1]), Comp.mol, Comp.label, accordance, score, findedrate])
                else:
                    pass                
            #else:
             #   print('Error !!!')
              #  break # da implementare errore?? Return('Error')
        
        SaveOutput(output, spectrum_path, Diff='Diff('+search_property[0]+')')
