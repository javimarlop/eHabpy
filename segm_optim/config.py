# config.py
# Diccionario centralizado de configuración


GRASSDB = '/Users/javier/grassdata' # REQUIRED
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
# True  = Borra los archivos de control y resultados anteriores para empezar desde el principio.
# False = Continúa el procesamiento omitiendo las áreas ya completadas en csv/segm_done.csv.
FORCE_RESTART = True

# DISTANCIA DE BUFFER OPCIONAL (en las unidades del mapa, ej. metros):
# 0 = Bounding box ajustado exactamente al área de estudio.
# > 0 = Extiende el bounding box de la región de trabajo en la distancia especificada.
PA_BUFFER = 500

# ESTADÍSTICO DE AGREGACIÓN POR HFT:
# 'mean'   = Utiliza la media (numpy.mean)
# 'median' = Utiliza la mediana (numpy.median)
STAT_AGG = 'mean'
