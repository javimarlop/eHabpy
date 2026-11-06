# eHabitat+ (eHabpy) — Modernization notes

This code base was written for **Python 2.7**, **GRASS GIS 7** and an old **R**
stack (rgdal, PySAL 1.x, monolithic `pysal`). It has been ported to run on
current versions:

| Component | Old | New target |
|-----------|-----|------------|
| Python    | 2.7 | **3.10 – 3.12** |
| GRASS GIS | 7.0.6 | **8.3 / 8.4** |
| GDAL      | 1.x / 2.x | **3.x** |
| NumPy     | 1.x | **1.24 – 2.2** |
| SciPy     | old | **1.10+** |
| PySAL     | `pysal` 1.x | **`libpysal` + `esda`** (PySAL 2.x) |
| R spatial | `rgdal` (retired 2023) | **`sf`** |
| R plots   | ggplot2 (old) | **ggplot2 3.4+/4.x** |

Nothing about the *algorithm* was changed — the aim was a faithful port that
produces the same results. Two exceptions are behaviour-preserving fixes for
constructs that are simply illegal in Python 3 (documented below).

See `README.md` for install/run instructions, `requirements.txt` /
`environment.yml` for dependencies.

---

## 1. Python 2 → Python 3 (all `.py`)

* `print x` → `print(x)` (every statement).
* `xrange` → `range`.
* `from __future__ import division` removed (Python 3 already true-divides;
  every `/` kept its floating-point meaning).
* `np.genfromtxt(..., dtype='string')` → `dtype=str` (NumPy dropped the
  `'string'` alias).
* `np.int(...)` → `int(...)`, `np.float` → `float` (removed in NumPy 1.24).
* `np.random.random_integers(N)` → `np.random.randint(1, N+1)`
  (`random_integers` was removed; note it was inclusive of the upper bound).

## 2. Removed / relocated library APIs

* **joblib**: `from sklearn.externals.joblib import Parallel, delayed`
  → `from joblib import Parallel, delayed` (scikit-learn deleted the vendored
  copy in 0.23). *scikit-learn is no longer a dependency at all.*
* **SciPy stats**: `from scipy.stats import chisqprob` (removed) →
  `from scipy.stats.distributions import chi2` and `chisqprob(x, k)` →
  `chi2.sf(x, k)` (identical values). `ehab_optim.py` already used `chi2.sf`.
* **GDAL**: `from gdalconst import GA_ReadOnly` →
  `from osgeo.gdalconst import GA_ReadOnly` (`lcgc.py`).
* **GDAL exceptions**: added `gdal.UseExceptions()` / `ogr.UseExceptions()`
  (GDAL ≥ 3.7 warns, GDAL 4.0 will error, if they are not enabled).
* **gdal_merge.py**: the `os.system('gdal_merge.py ...')` call in
  `ehab_optim*.py` is now a `merge_rasters()` helper that calls
  `osgeo_utils.gdal_merge` programmatically (same behaviour, no dependence on
  the script being on `$PATH`), with an `os.system` fallback.
* `os.system('rm ...')` for the temporary mask file → `os.remove(...)`.

## 3. Behaviour-preserving numeric fix (`ehab_optim.py`, `ehab_optim_median.py`)

The original computed the similarity raster with:

```python
pmhh = np.where(pmh <= 0.001, None, pmh)   # object array with None
pmhhmax = pmhh.max()
hr11 = np.where(pmhh > 0, 1, 0)
```

In Python 2, `None` sorted below every float, so `.max()` and `> 0` worked. In
Python 3 both raise `TypeError`. Replaced `None` with `0.0`:

```python
pmhh = np.where(pmh <= 0.001, 0.0, pmh)    # float array
```

Every downstream use of `pmhh` filters on `> 0`, so results are **identical**
(the low-similarity pixels are excluded either way). Verified numerically.

## 4. De-duplication of `getmeanvar.py` / `getmedianvar.py`

The originals were **3,395 lines each** — nine byte-for-byte copies of the same
routine (`ehabitat1` … `ehabitat9`) differing only by a numeric suffix in two
CSV names, and `getmedianvar.py` was a copy of `getmeanvar.py` with
`np.mean` → `np.median`.

That logic now lives once in **`segm_optim/ehab_segm_core.py`**, parameterized by:

* `k`      — the segmentation threshold level (1..9), used in the CSV names;
* `aggfun` — `np.mean` (getmeanvar) or `np.median` (getmedianvar).

`getmeanvar.py` and `getmedianvar.py` are now thin wrappers. **The public API is
unchanged** — `run_batch_all()`, `run_park(paid)`, and the historical
`run_batch1`…`run_batch9` / `ehabitat1`…`ehabitat9` names all still work:

```python
from getmeanvar import *      # or: from getmedianvar import *
run_batch_all()
```

(`ehab_optim.py` / `ehab_optim_median.py` were kept as two separate files, as in
the original, because a `from module import *` wrapper cannot override a global
that the imported function reads from its own module. There the only
computational difference is the per-HFT similarity statistic — `np.mean` vs
`np.median` for `AveHFTSim`/`MedianHFTSim` — which is selected by the module
constants `AGG` / `SIM_LABEL`. `ehab_optim_median.py` is generated from
`ehab_optim.py`.)

## 5. GRASS GIS 7 → 8 (`segmentation_pca_par.py`, `subpas_loop_segm_optim.py`,
   `ecoregs/ecoreg_buffers_loop.py`)

* **Session start-up**: the hard-coded `GISBASE`/paths are gone. A `find_gisbase()`
  helper reads `$GISBASE` or `grass --config path`, and `init_grass()` calls the
  GRASS 8 signature `gsetup.init(gisdbase, location, mapset)` (falling back to the
  GRASS 7 4-argument form). Database/location are configurable through the
  `GRASSDBASE`, `GRASSLOC`, `GRASSMAPSET`, `GISBASE` environment variables.
* `g.mapset ... gisdbase=` → `dbase=` (correct GRASS 8 option name).
* `i.segment ... seed=` → `seeds=`.
* `v.out.ogr ... ola=NAME dsn='.'` (GRASS 6/7) →
  `output=NAME.shp format='ESRI_Shapefile' output_layer=NAME` (GRASS 8; `dsn`
  was removed).
* `v.db.addcolumn col=` → `columns=`, `v.db.update col=/qcol=` →
  `column=/query_column=`, `v.extract out=`/`v.to.rast out=`/`r.to.vect out=`/
  `r.out.gdal out=` → `output=`, `g.region vect=` → `vector=`,
  `v.to.rast ... val=`/`labelcol=` → `value=`/`label_column=`,
  `r.mask ... maskc=` → `maskcats=`, `v.in.ogr out=` → `output=`.
* **Mask handling**: the old `g.rename rast=MASK,masc` trick (and the GRASS 6
  `g.remove rast=MASK`) used to drop a mask is replaced by `r.mask -r`
  (wrapped in a `drop_mask()` helper that ignores the "no mask present" error).
* `g.mremove type=... pattern=...` (removed in GRASS 7) →
  `g.remove type=... pattern=... -f`; `typ=`/`patt=`/`rast` → `type=`/
  `pattern=`/`raster`.
* `Pool(...)` is now guarded by `if __name__ == '__main__':` (required for
  `multiprocessing` on modern Python) and the worker count is set by
  `EHAB_NPROC`.

## 6. PySAL 1.x → `libpysal` + `esda`, and GeoPandas dissolve (`moranvar.py`)

* `import pysal` → `import libpysal, esda`.
* `pysal.rook_from_shapefile(shp)` → `libpysal.weights.Rook.from_shapefile(shp)`.
* `pysal.Moran(y, w)` → `esda.Moran(y, w)` (`.I` unchanged).
* The polygon **dissolve** (previously fiona + shapely `unary_union` +
  `itertools.groupby`, which is fragile across fiona versions) is now
  `geopandas`:
  `gpd.read_file(shp).dissolve(by='segm_id', as_index=False, aggfunc='first')`
  — same result (union per `segm_id`, first attributes kept), far more robust.
* `mors`/`varis` are converted to NumPy arrays before the vectorised
  normalisation (`(mors-min(mors))/(max(mors)-min(mors))`), which a plain
  Python list does not support.

## 7. R: rgdal → sf, and ggplot2 updates

* **`moranvar_plots.R`**
  * `library(rgdal)` removed.
  * `readOGR(dsn='shp', layer=...)` → `sf::st_read('shp', layer=..., quiet=TRUE)`.
  * `writeOGR(obj, dsn='results', layer=..., driver="ESRI Shapefile")` →
    `sf::st_write(obj, 'results/<layer>.shp', delete_layer=TRUE, quiet=TRUE)`.
  * `plot(spdf, col=...)` → `plot(st_geometry(sfobj), col=...)` (an `sf` object
    otherwise facet-plots every column).
  * `ggplot2:::rescale01` (internal, removed) → a local `rescale01()` helper.
* **`CreateRadialPlot.R`**: the deprecated/`defunct` line aesthetic
  `geom_path(..., size=)` → `linewidth=` (ggplot2 3.4+/4.x). `geom_point`/
  `geom_text` keep `size=`.
* **`thdi.R`** (already used `sf`): wrapped the ecoregion / class reads in
  `st_make_valid()` so modern GEOS does not reject the (invalid) WWF/MEOW
  ecoregion polygons during `st_intersection()`.

## 8. Housekeeping

* The vendored **`tqdm.py`** copies (Python 2, and shadowing the pip package)
  were deleted; scripts now do `from tqdm import tqdm`. Install with
  `pip install tqdm` (or via `environment.yml`).

---

## Known items to verify with real data

The raster inputs (`inVars/*.tif`), the GRASS database (`global_MW` location)
and the pre-computed intermediate CSV/shapefiles are **not** in this repository,
so the GRASS- and data-dependent steps could not be executed end-to-end here.
What *was* verified: every Python file byte-compiles under Python 3.12; the
scientific core (Mahalanobis distance via joblib, `chi2.sf`, the `None→0.0`
fix, `scipy.ndimage` labelling) runs correctly under NumPy 2.4 / SciPy 1.17;
all three R scripts parse under R 4.x.

Please check the following against your data when you run the full pipeline:

1. **Column name `hclust_mean` vs `hclst_m`.** `moranvar_plots.R` writes the HFT
   class column as `hclust_mean` (the ESRI Shapefile driver truncates it to
   `hclust_mea`), but `subpas_loop_segm_optim.py` and `lcgc.py` read a column
   named `hclst_m`. This mismatch pre-dates the port; align the name (or switch
   the class output to GeoPackage) for the segmentation→similarity hand-off.
2. **GRASS map/column names** (`gcmask`, `cspas`, `ecoregs_moll`, `pre`,
   `eprsqrt`, `ndvimax2`, the `ehabplus_cs,rasterized_parks,javier` mapsets,
   the `wdpaid`/`eco_id` columns) are specific to the original DOPA database and
   must exist in your GRASS location.
3. **`v.out.ogr` append** (`flags='a'`) and the per-worker mapset creation in
   `segmentation_pca_par.py` should be validated against your GRASS 8 build.

The experimental scratch scripts under `pas/tmp/` and `segm_optim/tmp/` are
superseded duplicates of the main scripts and were left untouched (not part of
the documented pipeline).
