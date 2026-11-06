#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / PySAL 2.x (libpysal + esda) / GeoPandas (2024-2026)
#### Refactored for dynamic variables, PA support and Shapefile-CSV alignment
#### License: CC BY-SA 3.0

import numpy as np
import os
import sys
import libpysal
import esda
import geopandas as gpd

# Configuración dinámica desde config.py
try:
    from config import ENV_VARS, CLIP_TO_PA
    num_vars = len(ENV_VARS)
except ImportError:
    print("ERROR: No se encuentra config.py.")
    sys.exit(1)

print(f"Modo de análisis espacial: {'Estricto a PA (CLIP_TO_PA=True)' if CLIP_TO_PA else 'Bounding Box completo (CLIP_TO_PA=False)'}")

idx_sumareas = 2 + (num_vars * 2)

ecox_list0 = np.genfromtxt('csv/segm_done.csv', dtype=int)
ecox_list = np.unique(np.atleast_1d(ecox_list0))
mx = len(ecox_list)

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

            # 1. Cargar IDs de segmentos presentes en el CSV
            raw_hri_ids = np.genfromtxt(hriname, skip_header=1, usecols=1)
            raw_hri_ids = np.atleast_1d(raw_hri_ids)
            hri_ids = [str(int(x)) for x in raw_hri_ids if not np.isnan(x)]
            
            sumareas = np.genfromtxt(hriname, skip_header=1, usecols=(idx_sumareas))
            sumareas = np.atleast_1d(sumareas)

            # 2. Cargar Shapefile, estandarizar segm_id y disolver
            gdf = gpd.read_file(shpname)
            gdf['segm_id'] = gdf['segm_id'].apply(lambda x: str(int(float(x))) if str(x).replace('.','',1).isdigit() else str(x))
            dissolved = gdf.dissolve(by='segm_id', as_index=False, aggfunc='first')

            # 3. Alineación exacta: Filtrar y reordenar el GeoDataFrame según los IDs del CSV
            dissolved = dissolved.set_index('segm_id').reindex(hri_ids).reset_index()
            dissolved = dissolved.dropna(subset=['geometry'])

            if not os.path.isfile(shpname2):
                dissolved.to_file(shpname2)

            # 4. Crear matriz espacial directamente desde el GeoDataFrame alineado
            w = libpysal.weights.Rook.from_dataframe(dissolved, use_index=False)
            print('w.n es:', w.n)
            wpn = w.n
            
            if wpn > 1:
                print('segmentation threshold is:', k)
                sareas = np.sum(sumareas)
                thr.append(k)
                
                y30 = np.genfromtxt(hriname, skip_header=1, usecols=(i))
                y30 = np.atleast_1d(y30)
                
                # El cálculo de Moran ahora no falla por desajuste numérico
                mi = esda.Moran(y30, w)
                mm = abs(mi.I)
                print('M.I. is:', mm)
                
                i2 = i + (num_vars * 2) + 1
                y330 = np.genfromtxt(hriname, skip_header=1, usecols=(i2))
                y330 = np.atleast_1d(y330)
                
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

    # Consolidación final
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

os.system('Rscript moranvar_plots.R')
print("BATCH END")
