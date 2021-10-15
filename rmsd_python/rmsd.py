from pymol import cmd
import pandas as pd
from matplotlib import cm, colors
import matplotlib.pyplot as plot
import os, glob
import numpy as np
import colorsys

# Esta funcion genera colores, util para graficas
def get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors


# Dejo la carpeta actual como directorio de trabajo, util para regresar aquí desde cualquier lugar
workingDir = os.getcwd()

# os.walk('.') almacena el arbol completo del directorio actual, para poder acceder a la lista de directorios y archivos
proteins = next(os.walk('.'))[1]

#print(proteins)
rmsdDict = {}
#Se recorre cada directorio de cada proteina
for i in proteins:
    #print(i)
    # Se entra al directorio de la proteina
    os.chdir(workingDir + "\\" + i)
    # se almacena el nombre de la proteina de referencia
    protRef = i + ".pdb"
    rmsdList = []
    # glob.glob permite buscar archivos con comodines y los almacena en una lista, en el for se recorre la lista y se calcula el RMSD
    # usando pymol
    for j in glob.glob("*_unrelaxed.pdb"):
        protModel = j[:-4]
        # print("Reference: {} and model: {}".format(protRef, protModel))

        # Funciones de pymol, para cargar las proteinas, limpiar y calcular el RMSD
        cmd.load(protRef)
        cmd.remove('solvent')
        cmd.remove('hetatm')
        cmd.remove('inorganic')
        cmd.remove('not polymer.protein')
        cmd.load(j)
        rmsd = cmd.cealign(protModel, i).get('RMSD')
        # print(rmsd)
        
        # Los rmsd se almacenan en una lista
        rmsdList.append(rmsd)
        cmd.quit
    
    # Los rmsd se almacenan en un diccionario
    rmsdDict[i] = rmsdList
    # print(rmsdDict)    
    os.chdir(workingDir)

# Se convierte el diccionario en un dataframe para usar pandas
df = pd.DataFrame(rmsdDict)
print(df)

# Como se van a graficar 12 proteinas se le pide a la funcion 12 colores diferentes
colors = get_colors(13)

# A continuación se grafican los datos, ojo que en la grafica se ignora una proteina:
# 5x5s: Es una predicción muy mala
# Los RMSD son muy buenos, a excepción de 5x5s
labels = list(df.loc[:, ~df.columns.isin(['5x5s'])].columns.values)
boxplot = plot.boxplot(df.loc[:, ~df.columns.isin(['5x5s'])], patch_artist=True, labels=labels, showfliers=False)
for patch, color in zip(boxplot['boxes'], colors):
    patch.set_facecolor(color)
# boxplot.plot()
plot.show()
