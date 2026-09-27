# config.py
# Diccionario centralizado de variables ambientales
# Formato: 'nombre_variable': {'file': 'nombre_archivo.tif', 'nodata': valor_nulo}

ENV_VARS0 = ['precip,slope,ndwi,ndvimin,ndvimax,temp']


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
