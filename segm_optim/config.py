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

# ENTRADA DE LA SEGMENTACIÓN:
# True  = Segmenta los tres primeros ejes de una PCA estandarizada de las variables de
#         ENV_VARS0, calculada para cada área (necesita al menos 3 variables).
# False = Segmenta directamente las variables de ENV_VARS0, todas con el mismo peso, sin PCA
#         (las variables muy correlacionadas, p. ej. dos índices de vegetación, pesan más).
# En los dos casos i.segment escala cada capa a su rango en la región de trabajo del área,
# así que las unidades de las variables no influyen.
SEGM_PCA = True

# RADIO PARA UNIR SEGMENTOS PEQUEÑOS (en píxeles, así depende del tamaño de píxel):
# los segmentos con menos píxeles que minarea (raíz cuadrada del nº de píxeles con datos de
# la región) se unen al segmento más frecuente de los mayores que ellos a esta distancia o
# menos, en pasadas repetidas hasta que no cambia nada. Los que no tienen ningún segmento
# mayor tan cerca se unen después al segmento de minarea píxeles o más más cercano, esté a
# la distancia que esté. 5 = 5 x RESOLUTION (50 m con píxeles de 10 m).
SMALL_SEGM_RADIUS = 10

# ELECCIÓN DEL UMBRAL DE SEGMENTACIÓN (moranvar.py), entre los nueve de 0.1 a 0.9, con el I de
# Moran de las medias de los segmentos (MI, bajo = segmentos distintos de sus vecinos) y la
# varianza dentro de ellos (V, baja = segmentos homogéneos) de cada variable. Ver README, A5:
# 'f'              = Máximo de la función F de Johnson et al. (2015), como i.segment.uspo de
#                    GRASS, con MI y V reescalados a 0-1 entre los umbrales y THRESHOLD_ALPHA.
# 'gs'             = Mínimo de MI + V reescalados (Global Score de Espindola et al., 2006).
# 'gs_fixed'       = Mínimo de V / varianza total + |I de Moran|, sin reescalar entre los
#                    umbrales: no depende de qué umbrales se prueben.
# 'crossing'       = Primer cruce de las curvas de MI y V reescaladas menos un umbral (el método
#                    original de eHabitat+).
# 'crossing_first' = Primer cruce de las curvas, sin restar un umbral.
# Un número        = Umbral fijo (p. ej. 0.5).
# En la prueba con Guajares, 'f' con alfa 0.5 fue el más estable y eligió el umbral cuyos HFTs
# explican más varianza de los píxeles.
THRESHOLD_METHOD = 'f'

# PESO DE V FRENTE A MI EN 'f': 1 = el mismo peso; 0.5 = V pesa la mitad (umbrales más altos,
# menos segmentos); 2 = V pesa el doble (umbrales más bajos, más segmentos).
THRESHOLD_ALPHA = 0.5

# VECINOS PARA EL I DE MORAN:
# 'rook'  = Segmentos que comparten un borde.
# 'queen' = También los que solo se tocan en una esquina: con muchos píxeles sin datos entre los
#           segmentos deja menos segmentos sin vecinos.
MORAN_NEIGHBOURS = 'rook'

# CLASIFICACIÓN DE LOS SEGMENTOS DEL UMBRAL ELEGIDO EN HFTs (moranvar_plots.R), con el método de Ward:
# 'ward_area'        = Ward ponderado por el nº de píxeles de cada segmento, sobre su media (o mediana,
#                      STAT_AGG) de cada variable tipificada con la desviación típica de los píxeles: los
#                      HFTs reducen al mínimo la varianza de los píxeles dentro de ellos.
# 'ward_area_kmeans' = 'ward_area' y después k-medias ponderadas: cada segmento pasa al HFT de media más
#                      próxima. Mejora algo la separabilidad, pero los HFTs de los cortes de SIMIL = 3 ya no
#                      están exactamente anidados.
# 'ward'             = Ward sin ponderar sobre las medias y las varianzas de los segmentos reescaladas a
#                      0-1 (el método original de eHabitat+).
# En la prueba con Guajares 'ward_area' dio HFTs más separables que 'ward' con todos los umbrales. Su árbol
# tiene las uniones más altas al final, así que los cortes de SIMIL = 3 dan menos HFTs (README, A6).
CLASS_METHOD = 'ward_area'

# NÚMERO DE HFTs (corte del árbol de clasificación jerárquica en moranvar_plots.R):
# 1 = Corta en el percentil 25 de las alturas de fusión, con límite superior de HFTs.
# 2 = No agrupa segmentos: el número de HFTs lo fijan los límites superiores
#     (por superficie del área y nº de variables - 1).
# 3 = Corta el árbol al 10, 25, 50, 75 y 90 % de su altura máxima, sin límite superior,
#     de más HFTs (p10) a menos (p90); los resultados de cada corte se guardan en
#     results/<STUDY_AREA>/<SCENARIO>/p10, p25, p50, p75 y p90. Se sugiere el corte más
#     próximo al codo de la curva de R2 (r2_<id>_hfts.png), con _optim en sus archivos: más
#     HFTs añaden poca información.
SIMIL = 3

# ESTADÍSTICO DE AGREGACIÓN POR HFT:
# 'mean'   = Utiliza la media (numpy.mean)
# 'median' = Utiliza la mediana (numpy.median)
STAT_AGG = 'mean'
