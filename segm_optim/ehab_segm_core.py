#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / GDAL 3+ (2024) - REFACTORED FOR DYNAMIC VARIABLES
#### Statistics of all the segments at once from one segment raster per area and threshold (2026)
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
import hft_state
import sys
from osgeo import ogr, gdal

# Importar configuración global de variables
try:
    from config import ENV_VARS, RESOLUTION
except ImportError:
    print("ERROR: No se encuentra config.py. Asegúrate de crearlo en el mismo directorio.")
    sys.exit(1)

gdal.UseExceptions()
ogr.UseExceptions()

res = RESOLUTION # tamaño de píxel de las variables de ../inVars (rellenado por setup.py)
_thr_values = None  # {k: value} of the thresholds of the segmentation (hft_state.run_thresholds), read when needed
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


class SegmentPixels:
    """The pixels of the segments of a segment raster (LABELS: 2-D, one label per segment, 0 or less = no segment), in
    row-major order: their positions in the raster (idx), their labels (lab) and the number of pixels of every label
    (n, indexed by label). The statistics of all the segments are computed from them at once, variable by variable,
    with the same results as the per-segment computation of eHabpy up to 1.5 (one mask tiffs/pa_<segm_id>.tif per
    segment, whose pixels were taken in row-major order too)."""

    def __init__(self, labels):
        flat = np.asarray(labels).ravel()
        self.idx = np.flatnonzero(flat > 0)
        self.lab = flat[self.idx].astype(np.intp)  # bincount works with intp: converted once, not in every call
        self.size = int(self.lab.max()) + 1 if self.lab.size else 1
        self.n = np.bincount(self.lab, minlength=self.size)
        self._grouped = None

    def grouped(self):
        """(order, start): positions in lab of the pixels of every segment, one segment after the other and each in
        row-major order (a stable sort by label), and where each label starts in order. Computed once, when needed."""
        if self._grouped is None:
            order = np.argsort(self.lab, kind='stable')
            start = np.cumsum(self.n) - self.n
            self._grouped = (order, start)
        return self._grouped

    def stats(self, v, aggfun=np.mean):
        """Statistic (AGGFUN: np.mean, np.median or any function of a 1-D array), population variance and 'all NaN'
        of variable V in every segment, as arrays indexed by label (NaN for labels without pixels).
        V: float64 values of the pixels self.idx, NaN = no data. Its gaps are filled in place as before: by linear
        interpolation (np.interp) along the pixels of their segment in row-major order (only the segments that
        have gaps are looped over); a segment without any data in V is 'all NaN'."""
        lab, n = self.lab, self.n
        nan = np.isnan(v)
        allnan = np.zeros(self.size, dtype=bool)
        if nan.any():
            nnan = np.bincount(lab[nan], minlength=self.size)
            allnan = (nnan == n) & (n > 0)
            gaps = np.flatnonzero((nnan > 0) & ~allnan)
            if gaps.size:
                order, start = self.grouped()
                for s in gaps:
                    pos = order[start[s]:start[s] + n[s]]
                    seg = v[pos]
                    m = np.isnan(seg)
                    seg[m] = np.interp(np.flatnonzero(m), np.flatnonzero(~m), seg[~m])
                    v[pos] = seg
        del nan

        with np.errstate(invalid='ignore', divide='ignore'):
            mean = np.bincount(lab, weights=v, minlength=self.size) / n
            # population variance (np.var) in two passes, sum of the squared deviations from the mean of the segment
            # (not E[x²] - mean², which loses the precision of small variances of large values)
            dev = mean[lab]
            np.subtract(v, dev, out=dev)
            np.square(dev, out=dev)
            var = np.bincount(lab, weights=dev, minlength=self.size) / n
        del dev

        stat = mean if aggfun is np.mean else self._per_segment(v, aggfun, allnan)
        return stat, var, allnan

    def _per_segment(self, v, aggfun, allnan):
        """AGGFUN (np.median or another function) of the values of V of every segment with data, in row-major order
        as before: grouped() puts them one segment after the other, a contiguous slice each, so the median is exactly
        that of np.median (the middle value, or the mean of the two middle values when the number of pixels is even).
        A sort of all the pixels by (segment, value) (np.lexsort) gives the same medians but took 15-20 s per variable
        on a grid of 17 million pixels; this loop about a second, since the grouping is computed once for all the
        variables of the area and threshold."""
        order, start = self.grouped()
        sv = v[order]
        stat = np.full(self.size, np.nan)
        for s in np.flatnonzero((self.n > 0) & ~allnan):
            stat[s] = float(aggfun(sv[start[s]:start[s] + self.n[s]]))
        return stat


def segment_stats(labels, values, aggfun=np.mean, nodata=None):
    """Statistics of a variable in every segment of a segment raster, for a single variable (see SegmentPixels):
    LABELS and VALUES are 2-D arrays of the same shape; NoData = NaN or NODATA in VALUES.
    Returns (n, stat, var, allnan), arrays indexed by label."""
    seg = SegmentPixels(labels)
    v = np.asarray(values).ravel()[seg.idx].astype(np.float64, copy=False)
    if nodata is not None:
        v[v == nodata] = np.nan
    return (seg.n,) + seg.stats(v, aggfun)


def _read_pixels(band_info, gt_seg, shape, idx):
    """Values of a variable of global_maps (float64, NoData -> NaN) at the pixels IDX (row-major positions) of a raster
    with geotransform GT_SEG and SHAPE (rows, columns) on the grid of the variables. The offset of the window, a whole
    number of pixels, is rounded so that a floating-point error does not move it by one pixel. Pixels outside the
    variable are NaN (no data)."""
    gt_var, band = band_info['gt'], band_info['band']
    rows, cols = shape
    xoff = int(round((gt_seg[0] - gt_var[0]) / res))
    yoff = int(round((gt_var[3] - gt_seg[3]) / res))
    x0, y0 = max(xoff, 0), max(yoff, 0)
    x1, y1 = min(xoff + cols, band.XSize), min(yoff + rows, band.YSize)
    if (x0, y0, x1, y1) == (xoff, yoff, xoff + cols, yoff + rows):
        v = band.ReadAsArray(xoff, yoff, cols, rows).ravel()[idx].astype(np.float64, copy=False)
    else:
        win = np.full(shape, np.nan)
        if x1 > x0 and y1 > y0:
            win[y0 - yoff:y1 - yoff, x0 - xoff:x1 - xoff] = band.ReadAsArray(x0, y0, x1 - x0, y1 - y0)
        v = win.ravel()[idx]
        del win
    nodata = band_info['nodata']
    if nodata is not None:
        v[v == nodata] = np.nan
    return v


def threshold_label(k):
    """Value of the threshold k of the segmentation, for the messages (csv/segm_thresholds.csv: k is the index that
    names the files)."""
    global _thr_values
    if _thr_values is None:
        _thr_values = hft_state.run_thresholds()
    return hft_state.threshold_label(_thr_values[k]) if k in _thr_values else f'k = {k}'


def segment_raster_files(paid, k, base=''):
    """Segment raster of area PAID and threshold K written by segmentation_pca_par.py (labels = cat of the vector, 0
    outside the segments) and its table cat -> segm_id."""
    return (os.path.join(base, f'tiffs/segm_{paid}_{k}.tif'), os.path.join(base, f'csv/segm_{paid}_{k}_ids.csv'))


def read_segment_ids(path):
    """{segm_id: cat} of a table csv/segm_<id>_<k>_ids.csv ('cat segm_id')."""
    out = {}
    with open(path) as f:
        next(f, None)
        for line in f:
            p = line.split()
            if len(p) == 2:
                out[int(p[1])] = int(p[0])
    return out


def segment_rows(ecor, pa_list, aggfun=np.mean):
    """Rows of csv/park_<park>_hri_results<k>.csv for the segments PA_LIST (segm_id, sorted) of the segment raster of
    area park and threshold ECOR: 'ecoregion segm_id <var>pamean... <var>pavar... sumpamask <var>pavar2...', as
    before: in the order of PA_LIST, without the segments with no pixels (or that start left of or above the grid of
    the first variable), 'None' for a variable without any data in the segment, values written with repr()."""
    tif, ids_csv = segment_raster_files(park, ecor, nwpath)
    for f in (tif, ids_csv):
        if not os.path.isfile(f):
            sys.exit(f"ERROR: {os.path.relpath(f, nwpath)} not found: the segmentation was made by an older eHabpy "
                     "or did not finish. Run it again (python3 run_hft.py, or python3 segmentation_pca_par.py).")
    ds = gdal.Open(tif)
    gt_seg = ds.GetGeoTransform()
    labels = ds.GetRasterBand(1).ReadAsArray()
    shape = labels.shape
    ds = None
    seg = SegmentPixels(labels)
    del labels
    cat_of = read_segment_ids(ids_csv)
    print(f"Area {park}, threshold {threshold_label(int(ecor))}: {len(pa_list)} segments, {seg.idx.size} pixels")

    # As before, a segment that starts left of or above the grid of the first variable is skipped
    # (offset 0 = the segment starts in the first row or column of that raster)
    gt_ref = global_maps[next(iter(ENV_VARS))]['gt']
    xoff = int(round((gt_seg[0] - gt_ref[0]) / res))
    yoff = int(round((gt_ref[3] - gt_seg[3]) / res))
    skip = np.zeros(seg.size, dtype=bool)
    if xoff < 0 or yoff < 0:
        r, c = np.divmod(seg.idx, shape[1])
        first = np.iinfo(np.int64).max
        rmin, cmin = np.full(seg.size, first), np.full(seg.size, first)
        np.minimum.at(rmin, seg.lab, r)
        np.minimum.at(cmin, seg.lab, c)
        skip = (seg.n > 0) & ((cmin + xoff < 0) | (rmin + yoff < 0))
        del r, c

    # one pass per variable over all the segments (the progress: one unit per segment, as before)
    results = []
    for var in ENV_VARS:
        v = _read_pixels(global_maps[var], gt_seg, shape, seg.idx)
        results.append(seg.stats(v, aggfun))
        del v
        hft_state.progress('DONE', f'{len(pa_list) / len(ENV_VARS):.4f}')

    lines, missing = [], 0
    for sid in pa_list:
        cat = cat_of.get(int(sid))
        if cat is None or not 0 < cat < seg.size:
            missing += cat is None
            continue
        n = int(seg.n[cat])
        if n <= 0 or skip[cat]:
            continue
        means, vars1, vars2 = [], [], []
        for stat, var, allnan in results:
            if allnan[cat]:
                means.append('None')
                vars1.append('None')
                vars2.append('None')
            else:
                # without rounding: with 2 decimals the variances of indices between -1 and 1 (~0.01) lost precision
                vv = float(var[cat])
                means.append(repr(float(stat[cat])))
                vars1.append(repr(vv))
                vars2.append(repr(vv * n))
        lines.append(f"{ecor} {sid} {' '.join(means)} {' '.join(vars1)} {n} {' '.join(vars2)}\n")
    if missing:
        print(f"WARNING: {missing} segments of csv/park_{park}_{ecor}.csv are not in {os.path.relpath(ids_csv, nwpath)}")
    return lines


def ehabitat(ecor, nw, nwpathout, k, aggfun=np.mean):
    """Statistics of every segment of area park (global) and threshold ECOR (the segm_id of csv/park_<park>_<ecor>.csv)
    appended to csv/park_<park>_hri_results<k>.csv. They are computed from the segment raster of the area and threshold
    (tiffs/segm_<park>_<ecor>.tif) with NumPy, all the segments at once, variable by variable (segment_rows)."""
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
    pa_list = np.unique(pa_list0)  # sorted segm_id: the order of the rows

    if len(pa_list):
        lines = segment_rows(ecor, pa_list, aggfun)
        with open(csvname, 'a') as wb:
            wb.write(''.join(lines))
        print(f"{len(lines)} segments written to {os.path.relpath(csvname, outdir)}")

    with open(csvname1, 'a') as wb:
        wb.write(f"{ecor}\n")


def run_batch(k, aggfun=np.mean):
    eco = f'csv/ecoregs{k}.csv'
    if not os.path.isfile(eco):  # written by the segmentation for every threshold exported
        print(f"WARNING: {eco} not found: no segments of the threshold {threshold_label(k)} (k = {k}); is the segmentation up to date?")
        return
    eco_list0 = np.genfromtxt(eco, dtype=int)
    eco_list = np.unique(eco_list0)
    for pm in range(len(eco_list)):
        ecor = eco_list[pm]
        print(f"Threshold {threshold_label(int(ecor))} (k = {ecor})")
        ehabitat(ecor, '', '', k, aggfun)
    print(str(datetime.now()))
    print("BATCH END")


def levels():
	"""Threshold levels k to summarise: all those of the segmentation (csv/segm_thresholds.csv: the thresholds of
	SEGM_THRESHOLDS, 1-9 for the nine of 0.1-0.9), or only the fixed threshold of config.py."""
	only = hft_state.fixed_level()
	return [only] if only else list(hft_state.run_thresholds())


def remove_park_results(paid):
    """Borra las estadísticas de un área: ehabitat() añade filas a sus archivos, así que sin esto
    ejecutarlo otra vez duplicaría las filas."""
    for k in range(1, hft_state.MAX_THRESHOLDS + 1):
        f = f'csv/park_{paid}_hri_results{k}.csv'
        if os.path.isfile(f):
            os.remove(f)


def run_batch_all(aggfun=np.mean):
    for k in range(1, hft_state.MAX_THRESHOLDS + 1):
        if os.path.isfile(f'csv/ecoregs_done{k}.csv'):
            os.remove(f'csv/ecoregs_done{k}.csv')
    ecox_list0 = np.genfromtxt('csv/segm_done.csv', dtype=int)
    ecox_list = np.unique(ecox_list0)
    for pmx in range(len(ecox_list)):
        global park
        park = ecox_list[pmx]
        print(park)
        remove_park_results(park)
        for k in levels():
            f = f'csv/park_{park}_{k}.csv'
            hft_state.progress('UNITS', len(np.unique(np.atleast_1d(np.genfromtxt(f, dtype=str)))) if os.path.isfile(f) else 0)
        for k in levels():
            run_batch(k, aggfun)
    print(str(datetime.now()))
    print("BATCH END")


def run_park(paid, aggfun=np.mean):
    global park
    park = paid
    print(park)
    remove_park_results(park)
    for k in levels():
        run_batch(k, aggfun)
    print(str(datetime.now()))
    print("BATCH END")
