#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / PySAL 2.x (libpysal + esda) / GeoPandas (2024-2026)
#### Robust alignment for segm_id matching (fixes disappearing polygons)
#### License: CC BY-SA 3.0

import glob
import numpy as np
import os
import subprocess
import sys
import warnings
import libpysal
import esda
import geopandas as gpd

import hft_state

# Configuración dinámica desde config.py
try:
    from config import ENV_VARS, CLIP_TO_PA
    num_vars = len(ENV_VARS)
except ImportError:
    print("ERROR: No se encuentra config.py.")
    sys.exit(1)

# Elección del umbral y vecinos del I de Moran (ver THRESHOLD_METHOD en config.py)
THRESHOLD_METHOD = getattr(sys.modules['config'], 'THRESHOLD_METHOD', 'f')
THRESHOLD_ALPHA = float(getattr(sys.modules['config'], 'THRESHOLD_ALPHA', 0.5))
MORAN_NEIGHBOURS = getattr(sys.modules['config'], 'MORAN_NEIGHBOURS', 'rook')
METHODS = ('f', 'gs', 'gs_fixed', 'crossing', 'crossing_first')
FIXED = isinstance(THRESHOLD_METHOD, (int, float)) and not isinstance(THRESHOLD_METHOD, bool)
if not (THRESHOLD_METHOD in METHODS or (FIXED and 0.1 <= THRESHOLD_METHOD <= 0.9)):
    sys.exit(f"ERROR: THRESHOLD_METHOD in config.py must be {', '.join(map(repr, METHODS))} "
             f"or a threshold between 0.1 and 0.9, not {THRESHOLD_METHOD!r}")
if MORAN_NEIGHBOURS not in ('rook', 'queen'):
    sys.exit(f"ERROR: MORAN_NEIGHBOURS in config.py must be 'rook' or 'queen', not {MORAN_NEIGHBOURS!r}")

# Para si la segmentación o las estadísticas no están al día con config.py, palist.csv o los datos
hft_state.check_ready('moranvar')
hft_state.start('moranvar')

print(f"Modo de análisis espacial: {'Estricto a PA (CLIP_TO_PA=True)' if CLIP_TO_PA else 'Bounding Box completo (CLIP_TO_PA=False)'}")

idx_sumareas = 2 + (num_vars * 2)

ecox_list0 = np.genfromtxt('csv/segm_done.csv', dtype=int)
ecox_list = np.unique(np.atleast_1d(ecox_list0))
mx = len(ecox_list)

def clean_id_str(val):
    """Limpia cadenas e IDs numéricos quitando sufijos .0 sobrantes."""
    s = str(val).strip()
    if s.endswith('.0'):
        s = s[:-2]
    return s


def rescale(x, high_good=False):
    """Reescala cada columna (variable) a 0-1 entre los umbrales (filas); 0 si todos son iguales.
    high_good: 1 en el valor más bajo (para la varianza, 1 = segmentos más homogéneos)."""
    with warnings.catch_warnings(), np.errstate(invalid='ignore', divide='ignore'):
        warnings.simplefilter('ignore', RuntimeWarning)  # variable sin datos: nan
        rng = np.nanmax(x, axis=0) - np.nanmin(x, axis=0)
        out = (np.nanmax(x, axis=0) - x) / rng if high_good else (x - np.nanmin(x, axis=0)) / rng
    out[:, rng == 0] = 0
    return out


def nanmean_rows(x):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(x, axis=1)


def crossings(thr, x1, x2):
    """Cruces de dos curvas, interpolados linealmente sobre los umbrales."""
    above = (x1 > x2).astype(int)
    ip = np.where(np.diff(above) != 0)[0]
    s1, s2 = x1[ip + 1] - x1[ip], x2[ip + 1] - x2[ip]
    t = (x2[ip] - x1[ip]) / (s1 - s2)
    return thr[ip] + t * (thr[ip + 1] - thr[ip])


def choose_threshold(ks, mi, var, vartot):
    """Umbral elegido (índice k, 1-9) y puntuación de cada umbral (None con 'crossing' o un umbral fijo)
    según THRESHOLD_METHOD. mi, var y vartot: una fila por umbral de ks y una columna por variable."""
    m3, v3 = rescale(mi), rescale(var, high_good=True)
    mi_n, var_n = nanmean_rows(m3), nanmean_rows(v3)
    nearest = lambda k: int(ks[np.argmin(np.abs(ks - k))])  # umbral evaluado más próximo
    score = None
    if FIXED:
        return nearest(round(THRESHOLD_METHOD * 10)), score
    if THRESHOLD_METHOD in ('crossing', 'crossing_first'):
        xp = crossings(ks.astype(float), mi_n, var_n)
        if len(xp) == 0:
            print("WARNING: the Moran's I and variance curves do not cross: threshold 0.5 used")
            return nearest(5), score
        return nearest(np.round(xp[0]) - (1 if THRESHOLD_METHOD == 'crossing' else 0)), score
    if THRESHOLD_METHOD == 'gs':
        score = nanmean_rows(m3 + (1 - v3))  # menor = mejor
    elif THRESHOLD_METHOD == 'gs_fixed':
        score = nanmean_rows(var / vartot + mi)  # menor = mejor
    else:
        nv, nm = var_n, 1 - mi_n  # 1 = mejor en los dos
        a2 = THRESHOLD_ALPHA ** 2
        den = a2 * nm + nv
        with np.errstate(invalid='ignore', divide='ignore'):
            score = np.where(den > 0, (1 + a2) * nv * nm / den, 0)  # mayor = mejor
    if np.isnan(score).all():
        print("WARNING: no threshold could be scored: threshold 0.5 used")
        return nearest(5), score
    best = np.nanargmax(score) if THRESHOLD_METHOD == 'f' else np.nanargmin(score)
    return int(ks[best]), score


for pmx in range(0, mx):
    park = ecox_list[pmx]
    print('park id is:', park)

    csvname = 'csv/' + str(park) + '_movar_segmsum.csv'
    csvnamem = 'csv/' + str(park) + '_moran_segm.csv'
    csvnamev = 'csv/' + str(park) + '_var_segm.csv'
    csvnamem2 = 'csv/' + str(park) + '_moran_mean.csv'
    csvnamev2 = 'csv/' + str(park) + '_var_mean.csv'
    csvname2 = 'csv/' + str(park) + '_movar_thresholds.csv'
    csvname5 = 'csv/' + str(park) + '_movar_results.csv'
    csvthr = 'csv/' + str(park) + '_threshold.csv'

    # Los archivos de esta área se empiezan de cero: abajo se les añaden filas, y al ejecutarlo
    # otra vez se duplicarían
    for f in [csvname, csvnamem, csvnamev, csvnamem2, csvnamev2, csvname2, csvname5, csvthr] + \
            glob.glob('shp/park_segm_' + str(park) + '_[1-9]_diss.*'):
        if os.path.exists(f):
            os.remove(f)

    # Una fila por umbral evaluado y una columna por variable
    thr, nsegs, nisol, mors, varis, vartots = [], [], [], [], [], []

    for k in range(1, 10):
        shpname = 'shp/park_segm_' + str(park) + '_' + str(k) + '.shp'
        shpname2 = 'shp/park_segm_' + str(park) + '_' + str(k) + '_diss.shp'
        hriname = 'csv/park_' + str(park) + '_hri_results' + str(k) + '.csv'

        if not os.path.isfile(hriname) or not os.path.isfile(shpname):
            continue

        # 1. Cargar datos del CSV asegurando IDs formateados como string limpio
        raw_hri_ids = np.genfromtxt(hriname, skip_header=1, usecols=1, dtype=str)
        raw_hri_ids = np.atleast_1d(raw_hri_ids)
        hri_ids_all = [clean_id_str(x) for x in raw_hri_ids if str(x) != 'nan']
        # Todas las columnas como números ('None' -> nan)
        hri_all = np.atleast_2d(np.genfromtxt(hriname, skip_header=1))

        # 2. Cargar Shapefile y disolver geometrías por segm_id estandarizado
        gdf = gpd.read_file(shpname)
        gdf['segm_id'] = gdf['segm_id'].apply(clean_id_str)
        dissolved = gdf.dissolve(by='segm_id', as_index=False, aggfunc='first')

        # 3. Emparejamiento seguro: Filtrar solo las coincidencias reales en ambos datos
        shp_ids_set = set(dissolved['segm_id'])
        valid_idx = [idx for idx, hid in enumerate(hri_ids_all) if hid in shp_ids_set]

        if len(valid_idx) == 0:
            print(f"⚠️ Sin coincidencias entre {hriname} y {shpname}")
            continue

        hri_ids = [hri_ids_all[idx] for idx in valid_idx]
        hri = hri_all[valid_idx]
        sumareas = hri[:, idx_sumareas]

        # Reordenar GeoDataFrame exactamente igual que las filas filtradas del CSV
        dissolved = dissolved.set_index('segm_id').loc[hri_ids].reset_index()

        # Guardar el shapefile disuelto completo
        if not os.path.isfile(shpname2):
            dissolved.to_file(shpname2)

        # 4. Vecinos: segmentos que comparten un borde (rook) o también una esquina (queen)
        wclass = libpysal.weights.Queen if MORAN_NEIGHBOURS == 'queen' else libpysal.weights.Rook
        w = wclass.from_dataframe(dissolved, use_index=False, silence_warnings=True)
        print(f'segmentation threshold {k}: {w.n} segments, {len(w.islands)} without {MORAN_NEIGHBOURS} neighbours')
        if w.n <= 1:
            continue

        sareas = np.sum(sumareas)
        mi_k, var_k, vt_k = [], [], []
        for i in range(2, 2 + num_vars):
            y30 = hri[:, i]                      # media (o mediana) de la variable en cada segmento
            y330 = hri[:, i + (num_vars * 2) + 1]  # varianza x nº de píxeles de cada segmento
            if np.isfinite(y30).all():
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)  # variable constante: I = nan
                    mi_k.append(abs(esda.Moran(y30, w, permutations=0).I))
            else:
                mi_k.append(np.nan)
            wv = np.sum(y330) / sareas
            var_k.append(wv)
            # Varianza total de los píxeles: dentro de los segmentos + entre sus medias
            grand = np.sum(sumareas * y30) / sareas
            vt_k.append(wv + np.sum(sumareas * (y30 - grand) ** 2) / sareas)

        thr.append(k)
        nsegs.append(w.n)
        nisol.append(len(w.islands))
        mors.append(mi_k)
        varis.append(var_k)
        vartots.append(vt_k)
        with open(csvname2, 'a') as wb:
            wb.write(f'{k} {w.n} {sareas}\n')

    if len(thr) == 0:
        continue

    ks = np.asarray(thr)
    mors, varis, vartots = (np.asarray(x, dtype=float) for x in (mors, varis, vartots))
    print('list of MIs (one column per variable):\n', mors)
    print('list of variances:\n', varis)

    # I de Moran y varianza reescalados entre los umbrales (una fila por variable), y su media
    m3, v3 = rescale(mors), rescale(varis, high_good=True)
    tot = (m3 + v3) / 2
    for f, x in [(csvnamem, m3), (csvnamev, v3), (csvname, tot)]:
        with open(f, 'a') as wb:
            for j in range(num_vars):
                wb.write(','.join(map(str, x[:, j])) + ',\n')
    for f, x in [(csvname5, nanmean_rows(tot)), (csvnamem2, nanmean_rows(m3)), (csvnamev2, nanmean_rows(v3))]:
        with open(f, 'a') as wb:
            for k, val in zip(ks, x):
                wb.write(f"{k} {val}\n")

    # Elección del umbral (moranvar_plots.R lee la columna selected)
    best, score = choose_threshold(ks, mors, varis, vartots)
    method = f"{THRESHOLD_METHOD} (alpha = {THRESHOLD_ALPHA})" if THRESHOLD_METHOD == 'f' else THRESHOLD_METHOD
    print(f'Selected threshold for park {park}: 0.{best} (THRESHOLD_METHOD = {method}, {MORAN_NEIGHBOURS} neighbours)')
    fmt = lambda x: repr(float(x)) if np.isfinite(x) else 'NA'
    cols = [nanmean_rows(mors), nanmean_rows(varis / vartots), nanmean_rows(m3), nanmean_rows(v3),
            score if score is not None else np.full(len(ks), np.nan)]
    with open(csvthr, 'w') as wb:
        wb.write('threshold segments isolated moran var_rel moran_norm var_norm score selected\n')
        for r, k in enumerate(ks):
            wb.write(f"{k} {nsegs[r]} {nisol[r]} {' '.join(fmt(c[r]) for c in cols)} {int(k == best)}\n")

hft_state.done('moranvar')

# moranvar_plots.R calls hft_state.py with the same Python (EHAB_PYTHON)
r = subprocess.run(['Rscript', 'moranvar_plots.R'], env=dict(os.environ, EHAB_PYTHON=sys.executable))
if r.returncode != 0:
    sys.exit(f"ERROR: Rscript moranvar_plots.R failed (exit code {r.returncode}): see the messages above")
print("BATCH END")
