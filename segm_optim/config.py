# config.py
# Diccionario centralizado de variables ambientales
# Formato: 'nombre_variable': {'file': 'nombre_archivo.tif', 'nodata': valor_nulo}

ENV_VARS0 = ['precip,slope,ndwi,ndvimin,ndvimax,temp']

RESOLUTION = 10
CLIP_TO_PA = False
COL_ID = 'cat'

ENV_VARS = {
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
