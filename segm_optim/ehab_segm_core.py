#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / GDAL 3+ (2024) - REFACTORED FOR DYNAMIC VARIABLES
#### License: CC BY-SA 3.0

from datetime import datetime
import numpy as np
import scipy.ndimage as nd
import os.path
import scipy
from scipy.linalg import cholesky, solve_triangular
from scipy.spatial import distance
from multiprocessing import cpu_count
import csv
import os
import sys
from osgeo import ogr, gdal

# Importar configuración global de variables
try:
    from config import ENV_VARS
except ImportError:
    print("ERROR: No se encuentra config.py. Asegúrate de crearlo en el mismo directorio.")
    sys.exit(1)

gdal.UseExceptions()
ogr.UseExceptions()

res = 10 
gmaps = 0
nwpath = ''
park = ''
global_maps = {}

def safelyMakeDir(d):
    try:
        os.makedirs(d)
        return True
    except OSError:
        if os.path.isdir(d):
            print(f"Directory exists: {d}. It will be used for the results.")
            return True
        else:
            print(f"Can't create directory: {d}.")
            return False

def _schedule(n, nproc):
    start = 0
    size = (n - start) // nproc
    while size > 100:
        yield slice(start, start + size)
        start += size
        size = (n - start) // nproc
    yield slice(start, n + 1)
    return

def _mahalanobis_distances_scipy(m, SI, X):
    n = X.shape[0]
    mahal = np.zeros(n)
    for i in range(X.shape[0]):
        x = X[i, :]
        mahal[i] = distance.mahalanobis(x, m, SI)
    return mahal


def initglobalmaps():
    global gmaps, global_maps
    indir = os.path.join(nwpath, '../inVars')
    print(f"Cargando mapas desde: {indir}")
    
    for var, props in ENV_VARS.items():
        filepath = os.path.join(indir, props['file'])
        try:
            ds = gdal.Open(filepath)
            global_maps[var] = {
                'ds': ds,
                'band': ds.GetRasterBand(1),
                'gt': ds.GetGeoTransform(),
                'nodata': props['nodata']
            }
            print(f"Cargada variable: {var} ({props['file']})")
        except Exception as e:
            print(f"Error al cargar {filepath}: {e}")
            sys.exit(1)
            
    print("Todas las variables globales han sido importadas dinámicamente")
    gmaps = 1


def ehabitat(ecor, nw, nwpathout, k, aggfun=np.mean):
    global nwpath, park
    if nw == '':
        nwpath = os.getcwd()
    else:
        nwpath = nw

    if gmaps == 0:
        initglobalmaps()
    outdir = nwpath if nwpathout == '' else nwpathout

    csvname1 = os.path.join(outdir, f'csv/ecoregs_done{k}.csv')
    if not os.path.isfile(csvname1):
        with open(csvname1, 'a') as wb:
            wb.write('None\n')

    csvname = os.path.join(outdir, f'csv/park_{park}_hri_results{k}.csv')
    if not os.path.isfile(csvname):
        # Generar cabeceras dinámicamente basadas en ENV_VARS
        vars_list = list(ENV_VARS.keys())
        mean_cols = " ".join([f"{v}pamean" for v in vars_list])
        var_cols = " ".join([f"{v}pavar" for v in vars_list])
        var2_cols = " ".join([f"{v}pavar2" for v in vars_list])
        header = f"ecoregion segm_id {mean_cols} {var_cols} sumpamask {var2_cols}\n"
        with open(csvname, 'a') as wb:
            wb.write(header)

    eco_csv = f'csv/park_{park}_{ecor}.csv'
    ecoparksf = os.path.join(nwpath, eco_csv)
    
    if not os.path.isfile(ecoparksf):
        return

    pa_list0 = np.genfromtxt(ecoparksf, dtype=int)
    pa_list = np.unique(pa_list0)
    
    for px in range(len(pa_list)):
        pa = pa_list[px]
        print(f"Procesando PA: {pa}")

        outfile = os.path.join(outdir, f'csv/park_{pa}_hri_results{ecor}.csv')
        pa_infile = f'tiffs/pa_{pa}.tif'
        pa4 = os.path.join(nwpath, pa_infile)

        dropcols = np.zeros(len(ENV_VARS), dtype=int)
        
        if not os.path.isfile(outfile) and os.path.isfile(pa4):
            src_ds_pa = gdal.Open(pa4)
            par = src_ds_pa.GetRasterBand(1)
            pa_mask0 = par.ReadAsArray(0, 0, par.XSize, par.YSize).astype(np.int32)
            pa_mask = pa_mask0.flatten()
            ind = pa_mask > 0
            
            sum_pa_mask = sum(pa_mask[ind])
            print(f"sum_pa_mask: {sum_pa_mask}")
            
            if sum_pa_mask <= 0:
                continue
                
            gt_pa = src_ds_pa.GetGeoTransform()
            
            # Usar la primera variable como referencia para el chequeo de límites original
            ref_var = list(ENV_VARS.keys())[0]
            gt_ref = global_maps[ref_var]['gt']
            xoff_test = int((gt_pa[0] - gt_ref[0]) / res)
            yoff_test = int((gt_ref[3] - gt_pa[3]) / res)
            
            if xoff_test > 0 and yoff_test > 0:
                results = {}
                
                # Bucle dinámico para procesar todas las capas ambientales
                for idx, var in enumerate(ENV_VARS.keys()):
                    band_info = global_maps[var]
                    gt_var = band_info['gt']
                    nodata = band_info['nodata']
                    
                    xoff = int((gt_pa[0] - gt_var[0]) / res)
                    yoff = int((gt_var[3] - gt_pa[3]) / res)
                    
                    var_bb0 = band_info['band'].ReadAsArray(xoff, yoff, par.XSize, par.YSize).astype(np.float32)
                    var_bb = var_bb0.flatten()
                    var_pa0 = var_bb[ind]
                    
                    var_pa = np.where(var_pa0 == nodata, float('NaN'), var_pa0)
                    mask2var = np.isnan(var_pa)
                    
                    if mask2var.all() == True:
                        dropcols[idx] = -idx
                        results[var] = {'mean': 'None', 'var': 'None', 'var2': 'None'}
                    else:
                        var_pa[mask2var] = np.interp(np.flatnonzero(mask2var), np.flatnonzero(~mask2var), var_pa[~mask2var])
                        var_pa = np.random.random_sample(len(var_pa),) / res + var_pa
                        
                        var_mean = round(aggfun(var_pa), 2)
                        var_var = round(np.var(var_pa), 2)
                        var_var2 = var_var * sum_pa_mask
                        
                        results[var] = {
                            'mean': str(var_mean),
                            'var': str(var_var),
                            'var2': str(var_var2)
                        }
                        print(f"pa {var} processado")

                # Generar líneas de texto dinámicamente
                means = ' '.join([results[v]['mean'] for v in ENV_VARS.keys()])
                vars1 = ' '.join([results[v]['var'] for v in ENV_VARS.keys()])
                vars2 = ' '.join([results[v]['var2'] for v in ENV_VARS.keys()])
                
                print("PA masked - results exported")
                with open(csvname, 'a') as wb:
                    var_line = f"{ecor} {pa} {means} {vars1} {sum_pa_mask} {vars2}\n"
                    wb.write(var_line)

    with open(csvname1, 'a') as wb:
        wb.write(f"{ecor}\n")


def run_batch(k, aggfun=np.mean):
    eco_list0 = np.genfromtxt(f'csv/ecoregs{k}.csv', dtype=int)
    eco_list = np.unique(eco_list0)
    for pm in range(len(eco_list)):
        ecor = eco_list[pm]
        print(ecor)
        ehabitat(ecor, '', '', k, aggfun)
    print(str(datetime.now()))
    print("BATCH END")


def run_batch_all(aggfun=np.mean):
    ecox_list0 = np.genfromtxt('csv/segm_done.csv', dtype=int)
    ecox_list = np.unique(ecox_list0)
    for pmx in range(len(ecox_list)):
        global park
        park = ecox_list[pmx]
        print(park)
        for k in range(1, 10):
            run_batch(k, aggfun)
    print(str(datetime.now()))
    print("BATCH END")


def run_park(paid, aggfun=np.mean):
    global park
    park = paid
    print(park)
    for k in range(1, 10):
        run_batch(k, aggfun)
    print(str(datetime.now()))
    print("BATCH END")
