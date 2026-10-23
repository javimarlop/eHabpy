# config.py
# Diccionario centralizado de configuración
# Las opciones REQUIRED las rellena 'python3 setup.py' a partir de ../inVars y pa/


GRASSDB = '~/grassdata' # REQUIRED: GRASS database folder (not the GRASS install folder)
GRASSLOC = 'Guajares_zona_quemada_dNBR2022_10m' # REQUIRED

ENV_VARS0 = ['Guajares_1_NDVImax_P95_S2_preincendio_2018_2021_sinAguaUrbanoCultivos_SIOSE_DERA_10m,Guajares_2_NDVImin_P05_S2_preincendio_2018_2021_sinAguaUrbanoCultivos_SIOSE_DERA_10m,Guajares_3_NDWI_Gao_mediana_S2_preincendio_2018_2021_sinAguaUrbanoCultivos_SIOSE_DERA_10m,Guajares_4_Pendiente_grados_CopernicusDEM_sinAguaUrbanoCultivos_SIOSE_DERA_10m,Guajares_5_TempMediaAnual_gradosC_CHELSA1981_2010_sinAguaUrbanoCultivos_SIOSE_DERA_10m,Guajares_6_PrecipAnual_mm_CHELSA1981_2010_sinAguaUrbanoCultivos_SIOSE_DERA_10m,Guajares_7_TWI_CopernicusDEM_sinAguaUrbanoCultivos_SIOSE_DERA_10m'] # REQUIRED
STUDY_AREA = 'Guajares_zona_quemada_dNBR2022_10m' # REQUIRED
COL_ID = 'cat' # REQUIRED

RESOLUTION = 10 # REQUIRED
CLIP_TO_PA = False

ENV_VARS = { # REQUIRED
    'Guajares_1_NDVImax_P95_S2_preincendio_2018_2021_sinAguaUrbanoCultivos_SIOSE_DERA_10m': {'file': 'Guajares_1_NDVImax_P95_S2_preincendio_2018-2021_sinAguaUrbanoCultivos-SIOSE-DERA_10m.tif', 'nodata': -9999}
    ,'Guajares_2_NDVImin_P05_S2_preincendio_2018_2021_sinAguaUrbanoCultivos_SIOSE_DERA_10m': {'file': 'Guajares_2_NDVImin_P05_S2_preincendio_2018-2021_sinAguaUrbanoCultivos-SIOSE-DERA_10m.tif', 'nodata': -9999}
    ,'Guajares_3_NDWI_Gao_mediana_S2_preincendio_2018_2021_sinAguaUrbanoCultivos_SIOSE_DERA_10m': {'file': 'Guajares_3_NDWI_Gao_mediana_S2_preincendio_2018-2021_sinAguaUrbanoCultivos-SIOSE-DERA_10m.tif', 'nodata': -9999}
    ,'Guajares_4_Pendiente_grados_CopernicusDEM_sinAguaUrbanoCultivos_SIOSE_DERA_10m': {'file': 'Guajares_4_Pendiente_grados_CopernicusDEM_sinAguaUrbanoCultivos-SIOSE-DERA_10m.tif', 'nodata': -9999}
    ,'Guajares_5_TempMediaAnual_gradosC_CHELSA1981_2010_sinAguaUrbanoCultivos_SIOSE_DERA_10m': {'file': 'Guajares_5_TempMediaAnual_gradosC_CHELSA1981-2010_sinAguaUrbanoCultivos-SIOSE-DERA_10m.tif', 'nodata': -9999}
    ,'Guajares_6_PrecipAnual_mm_CHELSA1981_2010_sinAguaUrbanoCultivos_SIOSE_DERA_10m': {'file': 'Guajares_6_PrecipAnual_mm_CHELSA1981-2010_sinAguaUrbanoCultivos-SIOSE-DERA_10m.tif', 'nodata': -9999}
    ,'Guajares_7_TWI_CopernicusDEM_sinAguaUrbanoCultivos_SIOSE_DERA_10m': {'file': 'Guajares_7_TWI_CopernicusDEM_sinAguaUrbanoCultivos-SIOSE-DERA_10m.tif', 'nodata': -9999}
}

# ESCENARIO: nombre de la subcarpeta de los resultados, results/<STUDY_AREA>/<SCENARIO>
# (letras, números, - y _), para guardar aparte los resultados de distintas opciones.
# run_hft.py lo muestra y permite cambiarlo antes de empezar (o --scenario NOMBRE), y avisa
# si ya hay resultados de ese escenario hechos con otras opciones. '' = 'base'.
SCENARIO = 'base'

# REINICIAR DESDE CERO (segmentation_pca_par.py):
# True  = Borra los archivos de control y csv/shp/tiffs y repite la segmentación desde el principio.
# False = Continúa una segmentación interrumpida omitiendo las áreas ya completadas en csv/segm_done.csv
#         (un área interrumpida a medias se rehace desde cero, borrando antes sus restos), y no
#         hace nada si ya terminó. Si desde entonces han cambiado las opciones de config.py,
#         palist.csv o los datos (run_state.json), o si csv/shp/tiffs son de otra área de
#         estudio, empieza desde cero.
# moranvar_plots.R mueve los resultados anteriores del escenario a
# results/<STUDY_AREA>/<SCENARIO>_previous; los de otros escenarios y otras áreas nunca se tocan.
FORCE_RESTART = False

# DISTANCIA DE BUFFER OPCIONAL (en las unidades del mapa, ej. metros; se redondea hacia arriba
# a un número entero de píxeles):
# 0 = Bounding box ajustado exactamente al área de estudio.
# > 0 = Extiende el bounding box de la región de trabajo en la distancia especificada.
PA_BUFFER = 500

# RADIO PARA UNIR SEGMENTOS PEQUEÑOS (en píxeles, así depende del tamaño de píxel):
# los segmentos con menos píxeles que minarea (raíz cuadrada del nº de píxeles con datos de
# la región) se unen al segmento más frecuente de los mayores que ellos a esta distancia o
# menos, en pasadas repetidas hasta que no cambia nada. Los que no tienen ningún segmento
# mayor tan cerca se unen después al segmento de minarea píxeles o más más cercano, esté a
# la distancia que esté. 5 = 5 x RESOLUTION (50 m con píxeles de 10 m).
SMALL_SEGM_RADIUS = 10

# NÚMERO DE HFTs (corte del árbol de clasificación jerárquica en moranvar_plots.R):
# 1 = Corta en el percentil 25 de las alturas de fusión, con límite superior de HFTs.
# 2 = No agrupa segmentos: el número de HFTs lo fijan los límites superiores
#     (por superficie del área y nº de variables - 1).
# 3 = Corta el árbol al 25, 50 y 75 % de su altura máxima, sin límite superior; los
#     resultados de cada corte se guardan en results/<STUDY_AREA>/<SCENARIO>/p25, p50 y p75.
SIMIL = 3

# ESTADÍSTICO DE AGREGACIÓN POR HFT:
# 'mean'   = Utiliza la media (numpy.mean)
# 'median' = Utiliza la mediana (numpy.median)
STAT_AGG = 'mean'
