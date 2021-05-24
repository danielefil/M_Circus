from pathlib import Path
import sys

## USER DEFINED LIBRARIES ##
# GUI layout
# Molecule class (for MM calculation) 
import molecule
import element
# Other functions
from Search_to_function_isopattern_10 import search_peak
import time


# Cartella contenente gli spettri
SP_folderPath = r'C:\Users\df426\Desktop\!!NEW\sp'
files = Path(SP_folderPath)

addotti = ['Hp', 'Cl'] #Lista degli addotti Hp --> per perdita H+, H
carica = -4 #Numero di cariche
# Percorso che contiene il database dei composti da cercare
DB_path = r'C:\Users\df426\Desktop\!!NEW\TEST.dat'
search_property = ['ppm', 40] #Modalita' di ricerca: ppm/dalton, tolleranza 
# Lista etichetta degli addotti: +H(+), -H(+)
addotti_label = ['-H(+)', '+Cl(-)']
filterValues = [False, 0, 0]


FilesList = []

for CurrentFile in files.glob("*.csv"):
    FilesList.append(CurrentFile)

start_time = time.time()
search_peak(FilesList, DB_path, addotti, carica,  search_property, addotti_label, Filtering=filterValues)
print(round(time.time() - start_time, 2), "seconds")

