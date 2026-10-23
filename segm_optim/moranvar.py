#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / PySAL 2.x (libpysal + esda) / GeoPandas (2024-2026)
#### Robust alignment for segm_id matching (fixes disappearing polygons)
#### License: CC BY-SA 3.0

import glob
import numpy as np
import os
import subprocess
import sys
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

    # Los archivos de esta área se empiezan de cero: abajo se les añaden filas, y al ejecutarlo
    # otra vez se duplicarían
    for f in [csvname, csvnamem, csvnamev, csvnamem2, csvnamev2, csvname2, csvname5] + \
            glob.glob('shp/park_segm_' + str(park) + '_[1-9]_diss.*'):
        if os.path.exists(f):
            os.remove(f)

    for i in range(2, 2 + num_vars):
        print('variable is:', i)
        mors = []
        varis = []
        thr = []
        
        for k in range(1, 10):
            wpn = 0
            shpname = 'shp/park_segm_' + str(park) + '_' + str(k) + '.shp'
            shpname2 = 'shp/park_segm_' + str(park) + '_' + str(k) + '_diss.shp'
            hriname = 'csv/park_' + str(park) + '_hri_results' + str(k) + '.csv'
            
            if not os.path.isfile(hriname) or not os.path.isfile(shpname):
                continue

            # 1. Cargar datos del CSV asegurando IDs formateados como string limpio
            raw_hri_ids = np.genfromtxt(hriname, skip_header=1, usecols=1, dtype=str)
            raw_hri_ids = np.atleast_1d(raw_hri_ids)
            hri_ids_all = [clean_id_str(x) for x in raw_hri_ids if str(x) != 'nan']
            
            sumareas_all = np.genfromtxt(hriname, skip_header=1, usecols=(idx_sumareas))
            sumareas_all = np.atleast_1d(sumareas_all)

            y30_all = np.genfromtxt(hriname, skip_header=1, usecols=(i))
            y30_all = np.atleast_1d(y30_all)

            i2 = i + (num_vars * 2) + 1
            y330_all = np.genfromtxt(hriname, skip_header=1, usecols=(i2))
            y330_all = np.atleast_1d(y330_all)

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
            y30 = y30_all[valid_idx]
            sumareas = sumareas_all[valid_idx]
            y330 = y330_all[valid_idx]

            # Reordenar GeoDataFrame exactamente igual que las filas filtradas del CSV
            dissolved = dissolved.set_index('segm_id').loc[hri_ids].reset_index()

            # Guardar el shapefile disuelto completo
            if not os.path.isfile(shpname2):
                dissolved.to_file(shpname2)

            # 4. Crear matriz espacial directamente desde las geometrías emparejadas
            w = libpysal.weights.Rook.from_dataframe(dissolved, use_index=False)
            print('w.n es:', w.n)
            wpn = w.n
            
            if wpn > 1:
                print('segmentation threshold is:', k)
                sareas = np.sum(sumareas)
                thr.append(k)
                
                mi = esda.Moran(y30, w)
                mm = abs(mi.I)
                print('M.I. is:', mm)
                
                wv = np.sum(y330) / sareas
                print('Sum of the variance is:', wv)
                
                mors.append(mm)
                varis.append(wv)
                
                with open(csvname2, 'a') as wb:
                    outxt = str(k) + ' ' + str(w.n) + ' ' + str(sareas) + '\n'
                    wb.write(outxt)
                    
        print('list of MIs:', mors)
        print('list of variances:', varis)
        
        if len(mors) > 0:
            mors = np.asarray(mors, dtype=float)
            varis = np.asarray(varis, dtype=float)
            
            mors_range = max(mors) - min(mors)
            varis_range = max(varis) - min(varis)
            
            m3 = (mors - min(mors)) / mors_range if mors_range > 0 else np.zeros_like(mors)
            v3 = (max(varis) - varis) / varis_range if varis_range > 0 else np.zeros_like(varis)
            tot = (m3 + v3) / 2

            with open(csvnamem, 'a') as wb:
                wb.write(','.join(map(str, m3)) + ',\n')

            with open(csvnamev, 'a') as wb:
                wb.write(','.join(map(str, v3)) + ',\n')

            with open(csvname, 'a') as wb:
                wb.write(','.join(map(str, tot)) + ',\n')

    # Consolidación final de resultados
    if os.path.isfile(csvname2):
        thrs = np.genfromtxt(csvname2, skip_header=0, usecols=(0))
        thrs = np.atleast_1d(thrs)
        segs = np.unique(thrs)
        
        for d in range(len(segs)):
            try:
                thr_vals = np.genfromtxt(csvname, delimiter=',', skip_header=0, usecols=(d))
                tm = np.nanmean(thr_vals)
                with open(csvname5, 'a') as wb:
                    wb.write(f"{segs[d]} {tm}\n")
                    
                m3d = np.genfromtxt(csvnamem, delimiter=',', skip_header=0, usecols=(d))
                m3m = np.nanmean(m3d)
                with open(csvnamem2, 'a') as wb:
                    wb.write(f"{segs[d]} {m3m}\n")
                    
                v3d = np.genfromtxt(csvnamev, delimiter=',', skip_header=0, usecols=(d))
                v3m = np.nanmean(v3d)
                with open(csvnamev2, 'a') as wb:
                    wb.write(f"{segs[d]} {v3m}\n")
            except Exception as e:
                print(f"Error procesando columna {d} en resultados finales: {e}")

hft_state.done('moranvar')

# moranvar_plots.R calls hft_state.py with the same Python (EHAB_PYTHON)
r = subprocess.run(['Rscript', 'moranvar_plots.R'], env=dict(os.environ, EHAB_PYTHON=sys.executable))
if r.returncode != 0:
    sys.exit(f"ERROR: Rscript moranvar_plots.R failed (exit code {r.returncode}): see the messages above")
print("BATCH END")
