# config.py
# Diccionario centralizado de configuración
# Las opciones REQUIRED las rellena 'python3 setup.py' a partir de ../inVars y pa/


GRASSDB = '~/grassdata' # REQUIRED: GRASS database folder (not the GRASS install folder)
GRASSLOC = 'Guajares_zona_quemada_dNBR2022_10m' # REQUIRED

ENV_VARS0 = ['gee_guajares_1_NDVI_S2_p95,gee_guajares_2_NDVI_S2_p5,gee_guajares_3_LSWI_S2_median,gee_guajares_4_slope,gee_guajares_5_temp,gee_guajares_6_precip'] # REQUIRED
STUDY_AREA = 'Guajares_zona_quemada_dNBR2022_10m' # REQUIRED
COL_ID = 'cat' # REQUIRED

RESOLUTION = 10 # REQUIRED
CLIP_TO_PA = False

ENV_VARS = { # REQUIRED
    'gee_guajares_1_NDVI_S2_p95': {'file': 'gee_guajares/1_NDVI_S2_p95.tif', 'nodata': -9999}
    ,'gee_guajares_2_NDVI_S2_p5': {'file': 'gee_guajares/2_NDVI_S2_p5.tif', 'nodata': -9999}
    ,'gee_guajares_3_LSWI_S2_median': {'file': 'gee_guajares/3_LSWI_S2_median.tif', 'nodata': -9999}
    ,'gee_guajares_4_slope': {'file': 'gee_guajares/4_slope.tif', 'nodata': -9999}
    ,'gee_guajares_5_temp': {'file': 'gee_guajares/5_temp.tif', 'nodata': -9999}
    ,'gee_guajares_6_precip': {'file': 'gee_guajares/6_precip.tif', 'nodata': -9999}
}

# CARPETA DE LAS VARIABLES DE ENTRADA: subcarpeta de ../inVars con los GeoTIFF de las variables,
# para probar otro conjunto de variables sin borrar ni mover las anteriores (p. ej. 'sin_twi' para
# ../inVars/sin_twi). '' = los GeoTIFF que están directamente en ../inVars.
# Al cambiarla, run_hft.py ejecuta antes setup.py, que importa en GRASS las variables de esa carpeta
# (con su nombre delante, p. ej. sin_twi_ndvi) y rellena ENV_VARS0, ENV_VARS y RESOLUTION, y después
# repite todos los pasos. Con SCENARIO = '' los resultados van al escenario con el nombre de la carpeta.
INVARS_FOLDER = 'gee_guajares'

# ESCENARIO: nombre de la subcarpeta de los resultados, results/<STUDY_AREA>/<SCENARIO>
# (letras, números, - y _), para guardar aparte los resultados de distintas opciones.
# run_hft.py lo muestra y permite cambiarlo antes de empezar (o --scenario NOMBRE), y avisa
# si ya hay resultados de ese escenario hechos con otras opciones.
# '' = el nombre de INVARS_FOLDER (otros caracteres pasan a _), o 'base' si está vacía.
SCENARIO = ''

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

# RECORTE DE LOS VALORES EXTREMOS (percentil, de 0 a menos de 50):
# 0   = Sin recorte.
# > 0 = Antes de la segmentación cada capa de i.segment (PC1-PC3 o las variables, según SEGM_PCA) y antes de
#       la clasificación cada variable se recortan a sus percentiles CLIP_PERCENTILE y 100 - CLIP_PERCENTILE
#       (p. ej. 0.5: los valores por debajo del percentil 0.5 pasan a valer el percentil 0.5, y los de por
#       encima del 99.5, el 99.5). i.segment escala cada capa a su rango y la clasificación tipifica cada
#       variable con su desviación típica, así que unos pocos píxeles extremos cambian cuánto pesa cada capa:
#       con el recorte, una variable con una cola larga (p. ej. TWI) no queda casi sin peso. Las gráficas y los
#       estadísticos de los HFTs usan siempre todos los valores originales.
CLIP_PERCENTILE = 0

# RADIO PARA UNIR SEGMENTOS PEQUEÑOS (en píxeles, así depende del tamaño de píxel):
# los segmentos con menos píxeles que minarea (raíz cuadrada del nº de píxeles con datos de
# la región) se unen al segmento más frecuente de los mayores que ellos a esta distancia o
# menos, en pasadas repetidas hasta que no cambia nada. Los que no tienen ningún segmento
# mayor tan cerca se unen después al segmento de minarea píxeles o más más cercano, esté a
# la distancia que esté. 5 = 5 x RESOLUTION (50 m con píxeles de 10 m).
SMALL_SEGM_RADIUS = 10

# ELECCIÓN DEL UMBRAL DE SEGMENTACIÓN (moranvar.py), entre los nueve de 0.1 a 0.9, con el I de
# Moran de las medias de los segmentos (MI, bajo = segmentos distintos de sus vecinos) y la
# varianza dentro de ellos (V, baja = segmentos homogéneos) de cada variable. Ver la wiki, How the HFT model works, A5:
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
# tiene las uniones más altas al final, así que los cortes de SIMIL = 3 dan menos HFTs (wiki, How the HFT model works, A6).
CLASS_METHOD = 'ward_area'

# NÚMERO DE HFTs (corte del árbol de clasificación jerárquica en moranvar_plots.R):
# 1 = Corta en el percentil 25 de las alturas de fusión, con límite superior de HFTs.
# 2 = No agrupa segmentos: el número de HFTs lo fijan los límites superiores
#     (por superficie del área y nº de variables - 1).
# 3 = Corta el árbol al 10, 25, 50, 75 y 90 % de su altura máxima, sin límite superior,
#     de más HFTs (p10) a menos (p90); los resultados de cada corte se guardan en
#     results/<STUDY_AREA>/<SCENARIO>/p10, p25, p50, p75 y p90. Se sugiere el corte más
#     próximo al codo de la curva de R2 (r2_<id>_hfts.png), con _optim en sus archivos: más
#     HFTs añaden poca información. hft_tree_<id>_<k>.png (y .csv) muestra en qué HFT de los cortes
#     más generales está cada HFT y en qué HFTs de los más detallados se divide.
SIMIL = 3

# ESTADÍSTICO DE AGREGACIÓN POR HFT:
# 'mean'   = Utiliza la media (numpy.mean)
# 'median' = Utiliza la mediana (numpy.median)
STAT_AGG = 'mean'

# ---------------------------------------------------------------------------------------------------------------------
# SIMILARIDAD DE LOS HFTs EN SU ECORREGIÓN Y MÉTRICAS DEL PAISAJE (similarity.py, después de run_hft.py).
# Usa las variables de ../inVars/<INVARS_FOLDER>/ecoregion (las mismas que las del área de estudio, con el mismo
# nombre) y los HFTs de un corte del escenario. Ver la wiki, Similarity and landscape metrics, y How the HFT model
# works, A12 y A13.

# ÁREA CON LA QUE SE COMPARAN LOS HFTs (su tamaño de píxel se detecta en sus GeoTIFF):
# ''    = Las ecorregiones del conjunto de variables de los HFTs (../inVars/<INVARS_FOLDER>/ecoregion).
# 'ruta' = Otra carpeta de ../inVars con las mismas variables (mismo nombre), p. ej. un área de comparación
#          descargada con get_vars_gee.py --area: 'gee_snevada/ecoregion'. Puede estar dentro o fuera del área de
#          estudio: si no la contiene, la media y la covarianza de los HFTs salen de su subcarpeta study_area (el área
#          de estudio con el píxel del área, que get_vars_gee.py descarga con ella).
ECOSIM_AREA = ''

# CORTE DEL ÁRBOL CUYOS HFTs SE COMPARAN CON LA ECORREGIÓN:
# ''    = El corte sugerido (_optim) de cada ID.
# 'p25' = Un corte concreto (p10, p25, p50, p75 o p90 con SIMIL = 3; con SIMIL 1 o 2 solo hay uno).
ECOSIM_CUT = ''

# TAMAÑO DE PÍXEL DE LA SIMILARIDAD (en las unidades del mapa, ej. metros):
# 0     = El píxel de las variables de la ecorregión (100 m con las de MonitorEO por defecto).
# > 0   = Un múltiplo entero de ese píxel (p. ej. 200, 500 o 1000): cada píxel grueso es la media de los píxeles
#         con datos de las variables de la ecorregión que contiene (necesita datos en al menos la mitad de su
#         superficie). Es más rápido y ocupa menos memoria, pero los HFTs pequeños se quedan con pocos píxeles
#         (con menos de max(5, nº de variables + 1) no se calculan). La media y la covarianza de cada HFT se
#         calculan siempre con los píxeles de este mismo tamaño, no con los de 10 m.
ECOSIM_RESOLUTION = 0

# UMBRAL DE LAS ÁREAS SIMILARES: estadístico de la similaridad de los píxeles del HFT (los que tienen
# similaridad mayor que 0); son similares los píxeles de la ecorregión con una similaridad igual o mayor.
# 'mean'   = La media (columna AveHFTSim, como ehab_optim.py).
# 'median' = La mediana (columna MedianHFTSim, como ehab_optim_median.py).
ECOSIM_STAT = 'mean'

# VECINOS DE LAS MANCHAS DE PÍXELES SIMILARES:
# 8 = También los que solo se tocan en una esquina (como eHabitat+): menos manchas, más grandes.
# 4 = Solo los que comparten un lado.
ECOSIM_NEIGHBOURS = 8

# FRACCIÓN MÍNIMA DE UN PÍXEL DE LA SIMILARIDAD CUBIERTA POR UN HFT para que el píxel sea de ese HFT (cada píxel es
# del HFT que más lo cubre, según los polígonos del corte). Los píxeles mixtos que no llegan no son de ningún HFT.
ECOSIM_MIN_COVER = 0.5

# POLÍGONOS DE LAS ECORREGIONES (para saber la ecorregión de cada HFT si hay varias, y para TerrEco del THDI):
# ''      = ../inVars/<INVARS_FOLDER>/ecoregion/ecoregions.gpkg (los descarga get_vars_gee.py).
# 'ruta'  = Un archivo propio (columna ECO_ID, o su primera columna de números enteros); relativa a la carpeta
#           del repositorio si no es absoluta.
ECOSIM_ECOREGIONS = ''

# POLÍGONOS DE LAS ECORREGIONES MARINAS (opcional, p. ej. MEOW), para MarEco del THDI (thdi.R):
# '' = Ninguno (MarEco = 0).
ECOSIM_MARINE = ''
