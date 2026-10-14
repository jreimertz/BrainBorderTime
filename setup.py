#! usr/bin/python

# Load packages
import scanpy as sc
import anndata as ad
import polars as pl

# Set up variables
file_path = "/data"
color_list = ["#779CAF", "#B8A7D1", "#F1AB86", "#ED5551", "#A7CFD1", 
               "#B6BCA6", "#E8A0D6", "#E0E0E0", "#ED9F58", "#C57B57",
               "#54CCC9", "#735684", "#7c1d00", "#E5BD57", "#6B86DB", 
               "#D1BFC5", "#76B041", "#B74E8A", "#98DFAF", "#CEDB56", 
               "#53685F", "#AAAAAA", "#ff0000", "#1F1F1F", "#3300bf"]
time_cols = ["#DCDCDC", "#87CBDD", "#3193AD", "#1C6C82"]

# altas object
atlas = ad.io.read_h5ad(f"{file_path}/atlas.h5ad")

marker_genes = pl.read_parquet(f"{file_path}/top_markers_scanpy.parquet")

scaled_exp = pl.read_parquet(f"{file_path}/top_20_scaled_exp.parquet")
