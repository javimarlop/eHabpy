# config.py
# Diccionario centralizado de configuración
# Las opciones REQUIRED las rellena 'python3 setup.py' a partir de ../inVars y pa/


GRASSDB = '~/grassdata' # REQUIRED: GRASS database folder (not the GRASS install folder)
GRASSLOC = 'ehab_guajares' # REQUIRED

ENV_VARS0 = ['precip,slope,ndwi,ndvimin,ndvimax,temp'] # REQUIRED
STUDY_AREA = 'perimetro_incendio' # REQUIRED
COL_ID = 'cat' # REQUIRED

RESOLUTION = 10 # REQUIRED
CLIP_TO_PA = False

ENV_VARS = { # REQUIRED
    'precip': {'file': 'precip.tif', 'nodata': -9999}
    ,'slope': {'file': 'slope.tif', 'nodata': -9999}
    ,'ndwi': {'file': 'ndwi.tif', 'nodata': -9999}
    ,'ndvimin': {'file': 'ndvimin.tif', 'nodata': -9999}
    ,'ndvimax': {'file': 'ndvimax.tif', 'nodata': -9999}
    ,'temp': {'file': 'temp.tif', 'nodata': -9999}
#    ,'ndvimax': {'file': 'ndvimax.tif', 'nodata': 65535.0}
#    ,'ndvimin': {'file': 'ndvimin.tif', 'nodata': 65535.0}
#    ,'herb': {'file': 'herb.tif', 'nodata': 255.0}
}

# REINICIAR DESDE CERO:
# True  = Borra los archivos de control, csv/shp/tiffs y los resultados anteriores de esta área
#         de estudio (results/<STUDY_AREA>) para empezar desde el principio.
# False = Continúa el procesamiento omitiendo las áreas ya completadas en csv/segm_done.csv
#         (un área interrumpida a medias se rehace desde cero, borrando antes sus restos).
#         Si csv/shp/tiffs son de otra área de estudio se empieza también desde cero.
# Los resultados de otras áreas de estudio (results/<otra área>) nunca se borran.
FORCE_RESTART = True

# DISTANCIA DE BUFFER OPCIONAL (en las unidades del mapa, ej. metros):
# 0 = Bounding box ajustado exactamente al área de estudio.
# > 0 = Extiende el bounding box de la región de trabajo en la distancia especificada.
PA_BUFFER = 500

# NÚMERO DE HFTs (corte del árbol de clasificación jerárquica en moranvar_plots.R):
# 1 = Corta en el percentil 25 de las alturas de fusión, con límite superior de HFTs.
# 2 = No agrupa segmentos: el número de HFTs lo fijan los límites superiores
#     (por superficie del área y nº de variables - 1).
# 3 = Corta el árbol al 25, 50 y 75 % de su altura máxima, sin límite superior; los
#     resultados de cada corte se guardan en results/<STUDY_AREA>/p25, p50 y p75.
SIMIL = 3

# ESTADÍSTICO DE AGREGACIÓN POR HFT:
# 'mean'   = Utiliza la media (numpy.mean)
# 'median' = Utiliza la mediana (numpy.median)
STAT_AGG = 'mean'
