import os
import scanpy as sc 
import anndata as ad
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
import geopandas as gpd
from shapely.geometry import Point, Polygon
import scipy.io
import scipy.sparse as sp
import gzip
from pathlib import Path


def construct_centers(cell_segmentation, a=100.0):
    # a = 100.0  # µm center-to-center pitch
    bounds = cell_segmentation.total_bounds  # [minx, miny, maxx, maxy] in µm
    minx, miny, maxx, maxy = bounds

    # Build hex centers covering the tissue bbox
    dx = a
    dy = np.sqrt(3)/2 * a

    xs = np.arange(minx - a, maxx + a, dx)
    ys = np.arange(miny - a, maxy + a, dy)

    centers = []
    rows = []
    cols = []

    r = 0
    for y in ys:
        offset = (r % 2) * (a/2)
        for c, x in enumerate(xs):
            centers.append((x + offset, y))
            rows.append(r)
            cols.append(c)
        r += 1

    centers_df = pd.DataFrame({"row": rows, "col": cols, "x": [p[0] for p in centers], "y": [p[1] for p in centers]})
    return centers_df, centers, rows, cols


def construct_spots(centers_df, centers, rows, cols, spot_r = 27.5):
    # spot_r = 27.5  # µm radius
    spot_geom = [Point(xy).buffer(spot_r, resolution=32) for xy in centers]
    spots = gpd.GeoDataFrame(
        centers_df.assign(spot_id=[f"r{r}c{c}" for r,c in zip(rows, cols)]),
        geometry=spot_geom,
        crs=None
    )
    return spots
    




def write_10x_mex_coordinate(vadata, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    umi = np.array(vadata.X.sum(axis=1)).ravel()

    keep = umi > 0
    adata_spots = vadata[keep].copy()
    X = adata_spots.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    # 10x expects genes x barcodes
    M = X.T.tocoo()

    # 10x matrices are integer counts
    M.data = M.data.astype(np.int64)

    # IMPORTANT: write COO coordinate format
    scipy.io.mmwrite(out_dir / "matrix.mtx", M)

    # gzip
    with open(out_dir / "matrix.mtx", "rb") as f_in, gzip.open(out_dir / "matrix.mtx.gz", "wb") as f_out:
        f_out.write(f_in.read())
    (out_dir / "matrix.mtx").unlink()

    # barcodes (spots)
    with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as f:
        for bc in adata_spots.obs_names.astype(str):
            f.write(bc + "\n")

    # features (genes): 3-column 10x format
    features = pd.DataFrame({
        0: adata_spots.var_names.astype(str),          # gene_id
        1: adata_spots.var_names.astype(str),          # gene_name
        2: ["Gene Expression"] * adata_spots.n_vars    # feature_type
    })
    with gzip.open(out_dir / "features.tsv.gz", "wt") as f:
        features.to_csv(f, sep="\t", header=False, index=False)

