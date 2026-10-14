#!/usr/bin/python

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ~                    all_functions.py                     ~
#' ~                                                         ~
#' ~            Author: Justin Colville-Reimertz             ~
#' ~       Description: Python script containing most        ~
#' ~              functions used to generate the             ~
#' ~                  Brain Border Time App                  ~
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###################
# Import packages #
###################

import enum
import os
import sys
import math
import scanpy as sc
import anndata as ad
import polars as pl
import pyarrow
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from matplotlib import rc_context
from typing import Optional, Union
from shiny import App, reactive, render, req, ui
from shiny.types import ImgData
from zarr import group

###################
# Manipulate Data #
###################

def apply_filters(adata: ad.AnnData, 
                  cells: list = [], 
                  timepoints: list = []
                  ) -> ad.AnnData:
    '''
    Subset AnnData object using defined filters

    Parameters
    ----------
    adata: AnnData object to apply filters to
    cells: list of cell types to filter. Optional
    timepoints: list of time points to filter. Optional

    Returns
    -------
    Subsetted AnnData object
    '''
    # Apply cell filters
    if len(cells) > 0:
        adata = adata[adata.obs['Broad_Cell_Annotation'].isin(cells)]

        # Sort broad and fine cell types
        adata.obs["broad_count"] = adata.obs.groupby(
            "Broad_Cell_Annotation", observed=True).cumcount()+1
        adata.obs["fine_count"] = adata.obs.groupby(
            "Fine_Cell_Annotation", observed=True).cumcount()+1
        adata.obs["broad_count"] = adata.obs.groupby(
            "Broad_Cell_Annotation", observed=True
            )['broad_count'].transform("max")
        adata.obs["fine_count"] = adata.obs.groupby(
            "Fine_Cell_Annotation", observed=True
            )['fine_count'].transform("max")
        adata.obs = adata.obs.sort_values(by = ["broad_count", "fine_count"], 
                                      ascending=False)
        
    # Apply time filters
    if len(timepoints) > 0:
        adata = adata[adata.obs['time'].isin(timepoints)]
    
    return adata

def get_group_col(col: str) -> str:
    '''
    Function to match nicely formatted group inputs with the appropriate meta 
    data column for use in various plotting functions

    Parameters
    ----------
    col: String with nicely formatted meta data column

    Returns
    -------
    AnnData meta data column
    '''
    if col != None:
        if col == "Cell Type":
            return "Broad_Cell_Annotation"
        elif col == "Time":
            return "time"
        elif col == "Mouse":
            return "mouse"
        else:
            return col
    else:
        return col

def build_scaled_exp_matrix(adata: ad.AnnData,
                            genes_list: list,
                            group: str,
                            groups_list: list,
                            #scale_mat: pl.DataFrame = pl.DataFrame() | None = None,
                            out_path: str = "scaled_exp.parquet"
                            ):
    '''
    Create a mean scaled expression matrix for a given list of genes

    Parameters
    ----------
    adata: AnnData object with scaled expression values
    genes_list: list of genes to include in the final matrix
    group: Meta data column to group scaled expression by
    groups_list: list of meta data values to include in the final matrix
    scale_mat: dataframe with the scaled matrix
    out_path: filepath to save the matrix to

    Returns
    -------
    Polars df with genes as row names and cell types as columns
    '''
    # Check if scaled expression matrix already exists at provided file path
    # If no file found, create a scaled expression matrix
    if not os.path.isfile(out_path):
        # Create scaled data matrix layer if doesn't already exist
        if "scaled" not in adata.layers:
            adata.layers["scaled"] = sc.pp.scale(adata, copy = True)

        # Subset the object for genes of interest
        adata = adata[:, adata.var_names.isin(genes_list)]
        # Fetch the scaled counts data
        scaled_df = adata.to_df(layer="scaled")
        # Confirm the indexes match up
        if all(scaled_df.index != adata.obs.index):
            sys.exit("Indexes don't match between adata and scaled data")
        # Add specified group column to scaled data
        scaled_df[group] = adata.obs[group]
        # Convert from Pandas format to Polars
        scaled_df = pl.from_pandas(scaled_df)
        # Transform from wide to long format
        merged_cells = scaled_df.unpivot(
            on=genes_list,
            value_name="scaled_exp", 
            variable_name="Gene",
            index = group)
        # Pivot on genes and get the mean scaled expression
        scaled_exp = merged_cells.pivot(
            on="Broad_Cell_Annotation", 
            index="Gene",
            aggregate_function="mean")
        
        groups_order = adata.obs[group].cat.categories.to_list()
        reordered_cols = groups_order.copy()
        reordered_cols.insert(0, "Gene")
        scaled_mat = scaled_exp.select(reordered_cols).sort(
            by=pl.col("Gene").sort_by(pl.Series(genes_list)))

        # Save the expression matrix to the provided file path
        scaled_mat.write_parquet(out_path)
        
    else:
        scaled_mat = pl.read_parquet(out_path)

    return scaled_mat

#########
# Plots #
#########

def get_style(mode):
    '''
    Function to set plot style depending on if the app is in light or dark mode

    Parameters
    ----------
    mode: dark mode switch id

    Returns
    -------
    Plotly theme
    '''
    if mode == "dark":
        plt.style.use("dark_background")
        return "plotly_dark"
    else:
        plt.style.use("default")
        return "plotly"

def build_umap(adata: ad.AnnData,
               mode: str,
               colors: list,
               group: str,
               split: str | None = None,
               ncol: int | None = None,
               static: bool = True,
               figsize: tuple = (15,15),
               save_fig: bool = False,
               plot_file: str = "umap.png",
               **kwargs
               ):
    '''
    Build a static or interactive UMAP plot

    Parameters
    ----------
    adata: AnnData object to plot
    mode: denotes whether app is in light or dark mode
    colors: list of colors to color points by
    group: meta data column to group the plot by
    split: meta data column to split the plot by
    ncol: number of columns to plot the output by if `split` provided. Optional
    static: boolean value denoting whether the plot should be static (`True`) 
    or interactive (`False`). If set to interactive, plot will be rendered 
    using Plotly to add hover information. Defaults to `True`
    coords_df: dataframe of UMAP coordinates and meta data information
    figsize: tuple to specify figure size. Defaults to `(15,15)`
    save_fig: boolean value denoting if the generated plot should be saved
    plot_file: string with filename to save fig to

    Returns
    -------
    Plot Figure
    '''
    # Set default colors using the provided group column
    group_vals = adata.obs[group].cat.categories
    colors_dict = {k:colors[i] for i,k in enumerate(group_vals)}

    # Set theme using the mode
    if mode == "dark":
        plt.style.use("dark_background")
        theme = "plotly_dark"
    else:
        plt.style.use("default")
        theme = "plotly_white"

    # Create UMAP plot
    if static:
        if split is None:
            with rc_context({"figure.figsize":figsize}):
                fig = sc.pl.umap(
                    adata,
                    color = group,
                    palette = colors_dict,
                    return_fig = True,
                    **kwargs
                )
                ax = fig.axes[0]
                ax.legend(ncol=1, bbox_to_anchor=(1,0.95),
                          fontsize = 25, frameon=False)
                if save_fig:
                    print("saving plot to file")
                    plt.savefig(fname=plot_file, bbox_inches="tight") 
        else:
            fig = split_umap(
                adata,
                group_by = group,
                split_by = split,
                cols = colors_dict,
                ncol = ncol,
                figsize=figsize,
                save_fig = save_fig,
                plot_file = plot_file,
                **kwargs
            )
    else:
        coords_df = get_UMAP(
            adata=adata,
            obs_cols=["Broad_Cell_Annotation", "Fine_Cell_Annotation", 
                      "time", "mouse"]
            )
        
        # Set point order for the legend
        cat_order = [c for c in adata.obs[group].cat.categories]

        # Create figure with new df
        fig = px.scatter(
            data_frame = coords_df,
            x = "UMAP1", 
            y = "UMAP2",
            color = group, 
            hover_data = ["Broad_Cell_Annotation", "Fine_Cell_Annotation",
                           "time", "mouse"],
            color_discrete_map = colors_dict,
            custom_data = ["barcode", "Broad_Cell_Annotation", 
                           "Fine_Cell_Annotation", "time", "mouse"],
            category_orders={group: cat_order},
            template = theme,
            #width=600,
            height=600,
            **kwargs
        )
        # Update the hover annotations
        fig.update_traces(
            hovertemplate = "<br>".join([
                "Barcode:   %{customdata[0]}",
                "Cell Type: %{customdata[1]}",
                "Subtype:   %{customdata[2]}",
                "Time:        %{customdata[3]}",
                "Mouse:     %{customdata[4]}"
            ]))
        
        fig.update_yaxes(
            scaleanchor="x",
            scaleratio=1,
        )

    return fig

def build_dot(adata: ad.AnnData,
              mode: str,
              colors: str,
              group_by: str,
              gene_df: pl.DataFrame | list,
              num_genes: int = 3,
              **kwargs):
    '''
    Function to build dot plot

    Parameters
    ----------
    adata: AnnData object to plot
    mode: denotes whether the app is in light or dark mode
    colors: continuous color palette
    group_by: obs column to group plot by
    gene_df: Dataframe of marker genes with assigned scores or list of genes
    num_genes: Number of marker genes to include in dot plot

    Returns
    -------
    Plot Figure
    '''
    # Set theme using the mode
    if mode == "dark":
        plt.style.use("dark_background")
    else:
        plt.style.use("default")

    # Determine if genes are marker genes or list of genes to be plotted
    if type(gene_df) is list:
        # If list of genes, just use that
        genes_list = gene_df.copy()
        num_genes = len(genes_list)
    elif type(gene_df) is pl.DataFrame:
        # Subset the marker genes using the provided number of genes
        subset_markers = gene_df.select(
            pl.all().top_k_by("scores", k=num_genes).over(
                "group", mapping_strategy="explode")
            )
        genes_list = subset_markers.get_column("names").unique(
            maintain_order=True).to_list()
    else:
        print(type(gene_df))
    
    # Set up plot axes
    fig,ax = plt.subplots(figsize=(10,6*num_genes))

    sc.pl.dotplot(adata, genes_list, group_by, 
                  swap_axes=True, cmap = colors, 
                  ax=ax, show=False)
    ax.tick_params("x", labelrotation=60)
    
    return ax

def build_feat(adata: ad.AnnData,
               mode: str,
               gene: str | list,
               colors: str,
               #na_color: str,
               split_by: str | None = None,
               **kwargs):
    '''
    Function to build feature plot

    Parameters
    ----------
    adata: AnnData object to plot
    mode: denotes whether the app is in light or dark mode
    gene: str or list of genes to plot
    colors: continuous color palette
    na_color: str denoting color for NA points
    split_by: str denoting obs column to facet plots by. Optional

    Returns
    -------
    Plot Figure
    '''
    # Check params
    if split_by is not None and split_by not in adata.obs.columns:
        msg = f"`split_by` col, {split_by} not found in `{adata}.obs`"
        sys.exit(msg)

    # Set theme using the mode
    if mode == "dark":
        theme = "plotly_dark"
    else:
        theme = "plotly_white"

    # If split by time, ensure timepoints in correct order
    if split_by == "time":
        cat_order = {"time": ["ZT0", "ZT6", "ZT12", "ZT18"]}
    else:
        cat_order = None

    # Pull UMAP coords as df and merge with relevant meta data cols
    umap_coords = get_UMAP(adata, 
                           obs_cols=["Broad_Cell_Annotation", 
                                     "Fine_Cell_Annotation", 
                                     "time", "mouse"])
    # Pull normalized expression data for specified genes
    gene_counts = sc.get.obs_df(adata, gene)
    gene_df = pl.from_pandas(gene_counts, include_index=True
                             ).rename({"None": "barcode"})

    meta = umap_coords.join(gene_df, on="barcode")

    # Generate feature plot(s)
    fig = px.scatter(
        data_frame = meta.sort(gene), 
        x = "UMAP1", y = "UMAP2",
        color = gene, 
        hover_data = [gene, "Broad_Cell_Annotation", 
                      "Fine_Cell_Annotation",
                      "time", "mouse"],
        color_continuous_scale=colors, 
        range_color=[0,math.ceil(meta[gene].max())],
        opacity=0.8,
        facet_col=split_by,
        category_orders=cat_order,
        facet_col_wrap=2,
        custom_data=[adata.obs.index, gene, "Broad_Cell_Annotation", 
                     "Fine_Cell_Annotation", "time", "mouse"],
        template=theme,
        height=600,
        **kwargs
    )
    
    # Update the hover annotations
    fig.update_traces(
        hovertemplate = "<br>".join([
            "Barcode:      %{customdata[0]}",
            "Expression: %{customdata[1]:.2f}",
            "Cell Type:    %{customdata[2]}",
            "Subtype:      %{customdata[3]}",
            "Time:           %{customdata[4]}",
            "Mouse:        %{customdata[5]}"
        ]))
    
    # Clean up subplot annotations
    if split_by is not None:
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

    fig.update_yaxes(
        scaleanchor="x",
        scaleratio=1,
    )

    return fig

def build_heatmap(scaled_exp: pl.DataFrame,
                  colors: str,
                  mode: str):
    '''
    Build a scaled expresison heatmap

    Parameters
    ----------
    scaled_exp: dataframe of scaled expression values
    colors: continous color palette
    mode: denotes whether the app is in light or dark mode

    Returns
    -------
    Plot figure
    '''
    # Set theme using the mode
    if mode == "dark":
        theme = "plotly_dark"
    else:
        theme = "plotly"

    #cell_types = scaled_exp.columns.isin(["Unnamed: 0", "Gene"]) == False
    # Create scaled exression heatmap
    fig = px.imshow(scaled_exp.drop("Gene"), 
                    x = scaled_exp.drop("Gene").columns, 
                    y= scaled_exp.get_column("Gene"),
                    labels={"x":"Cell Type", 
                            "y":"Gene", 
                            "color":"Scaled Expression"},
                    aspect="auto", 
                    #width = 1100, 
                    height = 600,
          color_continuous_scale=colors, 
          template=theme)
    fig.update_traces(
        hovertemplate = "%{y}</b>: %{z:.2f}",
        name = ""
    )
    fig.update_layout(hovermode = "x")
    fig.update_yaxes(showticklabels = False)
    fig.update_xaxes(tickangle = -45, dtick = {"size":15})

    return fig

def build_violin(adata: ad.AnnData,
                 mode: str,
                 gene: str,
                 colors: list,
                 group_by: str,
                 split_by: str | None = None,
                 **kwargs):
    '''
    Function to build violin plot

    Parameters
    ----------
    adata: AnnData object to plot
    mode: denotes whether the app is in light or dark mode
    colors: discrete color palette
    group_by: obs column to group plot by
    gene: Gene to plot

    Returns
    -------
    Plot Figure
    '''
    # Set theme using the mode
    if mode == "dark":
        plt.style.use("dark_background")
        linecolor = "lightgray"
    else:
        plt.style.use("default")
        linecolor = "black"

    # Set up df with selected gene
    df = pl.from_pandas(sc.get.obs_df(adata,
                   [gene, "time", group_by]))
    
    # Set cells order
    cells_order = [c for c in adata.obs[group_by].cat.categories]
    
    if split_by is not None:
        p = sns.catplot(df, x=group_by, y=gene, hue=group_by, kind="violin", 
                        palette=colors, inner="quart", linecolor=linecolor,
                        col="time", col_wrap=2, sharex=True, sharey=True,
                        col_order=["ZT0", "ZT6", "ZT12", "ZT18"],
                        hue_order = cells_order, order=cells_order, height=3, aspect=1)
        p.set_titles("{col_name}", y=0.9)

    else:
        p = sns.catplot(df, x=group_by, y=gene, hue=group_by, kind="violin", 
                        palette=colors, inner="quart", linecolor=linecolor,
                        hue_order = cells_order, order=cells_order, height=3, aspect=1)
        
    p.tick_params(axis="x", labelrotation=60, pad=1)
    p.set_axis_labels("", "Normalized Expression")
    p.tight_layout()

    return p
    

####################
# Helper Functions #
####################

def split_umap(adata: ad.AnnData, 
               group_by: str,
               split_by: str, 
               cols: dict,
               ncol: int | None = 2, 
               nrow: int | None = None, 
               figsize: tuple = (15,15),
               save_fig: bool = False,
               plot_file: str = "split_umap.png",
               **kwargs):
    '''
    Helper function to split umap plots

    Parmeters
    ---------
    adata: AnnData object to plot
    group_by: meta data column to group plots by
    split_by: meta data column to split plots by
    cols: dict with category values as keys and colors as values
    ncol: number of columns to plot the output by. Optional
    nrow: number of rows to plot the output by. Optional
    save_fig: boolean value denoting if the generated plot should be saved
    plot_file: string with filename to save fig to

    Return
    ------
    Plot figure
    '''
    categories = adata.obs[split_by].cat.categories
    if ncol is None and nrow is None:
        ncol = int(np.ceil(math.sqrt(len(categories))))
        nrow = int(np.ceil(math.sqrt(len(categories))))
    elif nrow is None and ncol is not None:
        nrow = int(np.ceil(len(categories) / ncol))
    elif ncol is None and nrow is not None:
        ncol = int(np.ceil(len(categories) / nrow))
    fig, axs = plt.subplots(nrow, ncol, 
                            figsize=(4*ncol, 4*nrow))
    axs = axs.flatten()
    for i, cat in enumerate(categories):
        ax = axs[i]
        adata_sub = adata[adata.obs[split_by] == cat]
        split_col = {cat:cols[cat]}
        with rc_context({"figure.figsize":figsize}):
            sc.pl.umap(adata_sub, ax=ax, color = group_by,
                    palette = split_col, show=False, 
                    title=cat, **kwargs)
            ax.set_title(cat, fontstyle = "normal", fontsize = "20",
                        fontweight = "bold")
            ax.get_legend().remove()

            #fig = plt.tight_layout()
            if save_fig:
                print("saving plot to file")
                plt.savefig(fname=plot_file, bbox_inches="tight") 
    
    return plt.tight_layout()

def get_UMAP(adata: ad.AnnData,
             obs_cols: list = []):
    '''
    Pull UMAP_X coordinates from adata object and merge with any desired obs 
    cols

    Parameters
    ----------
    adata: AnnData object
    obs_cols: list of columns found in `adata.obs`

    Returns
    -------
    Polars DataFrame
    '''
    # Confirm provided cols are in obs
    cleaned_cols = [col for col in obs_cols if col in adata.obs.columns]
    # Pull UMAP coords as a dataframe
    umap_coords = pl.DataFrame(adata.obsm["X_umap"], schema=["UMAP1","UMAP2"])
    # Add UMAP coords to selected meta data cols
    meta = pl.from_pandas(
        data=adata.obs[cleaned_cols], include_index=True
        ).rename({"None": "barcode"})

    return meta.hstack(umap_coords)

def get_image(img_file: str,
              width: str = "100%",
              height: str = "100%"):
    """ 
    load image from provided file path

    Parameters
    ----------
    img_file: file path for image to load

    Returns
    -------
    rendered image
    """
    img: ImgData = {"src": f"/data/figures/{img_file}",
                    "width": width,
                    "height": height}
    return img
