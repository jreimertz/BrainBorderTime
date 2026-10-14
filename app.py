#! usr/bin/python


# import packages
from networkx import is_empty
from shiny import App, reactive, render, render_ui, req, ui, module
from shinywidgets import output_widget, render_widget
import pyarrow

from all_functions import *
#from setup import file_path, color_list, time_cols, marker_genes, scaled_exp
from setup import *

import scanpy as sc
import plotly.express as px

# Set scanpy to not autoshow plots
sc.settings.autoshow = False

# Define modules
@module.ui
def nav_panel_ui(genePanel: str="", type: str=""):
    if type == "Feat":
        return ui.nav_panel(
            genePanel,
            output_widget(
                f"gene{type}",
                fill=True,
                fillable=True,
                height="600px"
            ),
            value=f"{genePanel}{type}Panel"
        )
    else:
        return ui.nav_panel(
            genePanel,
            ui.output_plot(
                f"gene{type}",
                fill=True,
                height="600px"
            ),
            value=f"{genePanel}{type}Panel"
        )

@module.server
def nav_panel_server(input, output, session, 
                     adata: ad.AnnData,
                     mode: str, 
                     colors: str, 
                     time_split: bool = True, 
                     groupby: str | None = None,
                     genePanel: str=""):
    
    @render_widget
    def geneFeat():
        if time_split:
            splitby = "time"
        else:
            splitby = None

        if type(colors) is str:
            return build_feat(
                adata=adata,
                mode=mode,
                gene=genePanel,
                colors=colors,
                split_by=splitby
            )
    
    @render.plot(alt="Selected Cells and Genes Violin Plots")
    def geneVln():
        if time_split:
            splitby = "time"
        else:
            splitby = None
        if type(colors) is list and groupby is not None:
            return build_violin(
                adata=adata,
                mode=mode,
                gene=genePanel,
                colors=colors,
                group_by=groupby,
                split_by=splitby
            )

# Define UI
app_ui = ui.page_navbar(
    # Information and landing page
    ui.nav_panel("Home", 
        ui.h1("Brain Border Time: A single-cell atlas of brain immune niches around the day-night cycle",
            align = "center"),
            ui.hr(),
            ui.layout_columns(
                ui.card(
                    ui.card(
                        ui.output_image("statUMAP"),
                        height="50vh",
                        fill=True
                    ),
                    ui.card(
                        ui.output_image("timeUMAP"),
                        fill=False,
                        height="20vh"
                    ),
                    ui.card_footer("Brain Border Time UMAP"),
                    #height="100%",
                    #width="100%",
                    full_screen=False
                ),
                ui.layout_columns(
                    ui.card(
                        ui.card_header(
                            "Background Information",
                            align = "center"),
                        ui.p("""Life on Earth is governed by the rhythmic 
                             rotation of the planet, dividing time into 
                             periods of light and dark to which most life 
                             forms are entrained. In mammals, time-keeping 
                             is genetically encoded by the circadian clock: a 
                             transcription-translation feedback loop lasting 
                             twenty-four hours that synchronizes cellular 
                             physiology to time-of-day. Core organismal 
                             functions are carried out under circadian 
                             regulation including even innate immunity, the 
                             body’s first line of defense against external 
                             pathogens. We set out to characterize how the 
                             brain's unique immune compartment changes around 
                             the day-night cycle via single-cell RNA sequencing 
                             (scRNAseq) of brain immune niches every six hours 
                             around a 24h clock. This web platform is intended 
                             to allow for exploration of the dataset.""",
                             style = "20px", align = "center"
                            )
                    ),
                    ui.card(
                        ui.card_header("Methods", align = "center"),
                        ui.p("""To comprehensively characterize how brain 
                             immune niches change around the clock, we 
                             performed single-cell RNA sequencing (scRNAseq) 
                             of the brain and leptomeninges of healthy adult 
                             mice at Zeitgeber time points ZT0, ZT6, ZT12, and 
                             ZT18. Zeitgeber Time (ZT) refers to the number of 
                             hours elapsed from the beginning of the light 
                             cycle; this approach thus sampled around a single 
                             24-hour light-dark cycle. Because mice are 
                             nocturnal, ZT0 and ZT6 – light cycle time points 
                             – represent the rest phase, and ZT12 and ZT18 
                             represent the active phase.""", 
                             style = "20px", align = "center"),
                        ui.p("""To sensitively sample rarer immune populations 
                             including BAMs and infiltrating leukocytes, we 
                             used fluorescence-activated cell sorting (FACS) in 
                             combination with acute labeling of circulating 
                             immune cells to de-enrich microglia while still 
                             excluding peripheral contamination. To arrest 
                             cellular transcription at the time of sample 
                             collection, prevent dissociation-induced gene 
                             signatures, and ensure transcriptional profiles 
                             were faithful to collection time, samples were 
                             treated with validated transcriptional and 
                             translational inhibitors throughout the collection 
                             pipeline.""", style = "20px", align = "center")
                    ),
                    row_heights=[3,4,3],
                    col_widths=12
                )
            ),
            # Citation information
            ui.card(
                ui.card_header("Citation", align = "center"),
                ui.p("""If you find this tool or data helpful please 
                        use the following citation:""", align = "center"),
                ui.p("[Citation here]", align = "center")
            )
    ),
    # Visualize the full dataset
    ui.nav_panel("Atlas",
        ui.h4("Visual Representation of the Full Atlas", 
                align = "center"),
        ui.hr(),
        ui.layout_columns(
            ui.accordion(
                ui.accordion_panel(
                    ui.h5("Explore specific cells and/or genes",
                          align = "center"),
                    ui.layout_columns(
                        ui.card(
                            ui.card_body(
                                ui.input_selectize(
                                    "selectCellViz",
                                    "Selct one or more cell types",
                                    choices={},
                                    multiple=True,
                                    #options=[maxItems=5]
                                ),
                                height="300px"
                            )
                        ),
                        ui.card(
                            ui.card_body(
                                ui.input_selectize(
                                    "selectGeneViz",
                                    "Select one or more genes",
                                    choices={},
                                    multiple=True,
                                    #options=[maxItems=5]
                                ),
                                height="300px"
                            )
                        ),
                        ui.card(
                            ui.card_body(
                                ui.input_action_button(
                                    "searchViz",
                                    "Explore Selected",
                                    width="100%"
                                ),
                                height="100px"
                            )
                        ),
                        #row_heights=4,
                        #col_widths=[4,4]
                    ),
                    width="100%",
                    value = "searchBarViz"
                ),
                open=False
            )
        ),
        ui.navset_pill(
            ui.nav_panel(
                "UMAP",
                ui.layout_columns(
                    ui.card(
                        ui.layout_sidebar(
                            ui.sidebar(
                                "Adjust Parameters",
                                ui.accordion(
                                    ui.accordion_panel(
                                        "Group Selection",
                                        ui.input_selectize(
                                            "selectGroupViz",
                                            "Select Grouping Category",
                                            choices={
                                                "Broad_Cell_Annotation":"Cell Type", 
                                                "time":"Time", 
                                                "mouse":"Mouse"
                                                }
                                        )
                                    ),
                                    #ui.accordion_panel(
                                    #    "Color Selection",
                                    #    ui.output_ui("colors")
                                    #),
                                    open=True,
                                    multiple=True
                                ),
                                position="right",
                                fill = True
                            ),
                            ui.card_body(
                                output_widget(
                                    "hoverUMAP",
                                    fill=True,
                                    height="600px"
                                ),
                            ),
                        ),
                        ui.card_footer("Brain Border Time UMAP"),
                        full_screen=True,
                        fill=True
                    ),
                    row_heights=12
                )
            ),
            ui.nav_panel(
                "Dot Plot",
                ui.layout_columns(
                    ui.card(
                        ui.layout_sidebar(
                            ui.sidebar(
                                "Adjust Parameters",
                                ui.accordion(
                                    ui.accordion_panel(
                                        "Color Palette",
                                        ui.input_selectize(
                                            "selectPaletteDot",
                                            "Select Color Palette",
                                            choices={
                                                "magma_r":"Dark Magma",
                                                "magma":"Light Magma",
                                                "plasma_r":"Dark Plasma",
                                                "plasma":"Light Plasma",
                                                "inferno_r":"Dark Inferno",
                                                "inferno":"Light Inferno",
                                                "viridis_r":"Dark Viridis",
                                                "viridis":"Light Viridis"
                                                }
                                        )
                                    ),
                                    ui.accordion_panel(
                                        "Marker Gene Count",
                                        ui.input_numeric(
                                            "numMarkers",
                                            "Select Number of Marker Genes",
                                            value=2,
                                            min=1,
                                            max=3
                                        )
                                    ),
                                    open=True,
                                    multiple=True
                                ),
                                position="right"
                            ),
                            ui.card_body(
                                ui.output_plot("fullDot",
                                                height="150vh")
                            ),
                            ui.card_footer("Brain Border Time Top Marker Genes")
                        ),
                        full_screen=True,
                    ),
                    row_heights=12
                )
            ),
            ui.nav_panel(
                "Heatmap",
                ui.layout_columns(
                    ui.card(
                        ui.layout_sidebar(
                            ui.sidebar(
                                "Adjust Parameters",
                                ui.accordion(
                                    ui.accordion_panel(
                                        "Color Palette",
                                        ui.input_selectize(
                                            "selectPaletteHeat",
                                            "Select Color Palette",
                                            choices=px.colors.named_colorscales()
                                        )
                                    ),
                                    open=True
                                ),
                                position="right",
                                fill = True
                            ),
                            ui.card_body(
                                output_widget(
                                    "fullHeatmap",
                                    fill=True),
                                height="600px",
                                #width="100%",
                            )
                        ),
                        ui.card_footer("Top 20 Marker Genes Heatmap"),
                        full_screen=True,
                        fill=True
                    ),
                    row_heights=12
                )
            )
        )

    ),
    # Explore specific cell and gene selections
    ui.nav_panel("Explore",
        ui.h1("Visually Explore the Atlas",
            align = "center"),
        ui.hr(),
        ui.layout_columns(
            ui.accordion(
                ui.accordion_panel(
                    ui.h5("Explore Specific cells and/or genes",
                          align="center"),
                    ui.layout_columns(
                        ui.card(
                            ui.card_body(
                                ui.input_selectize(
                                    "selectCellExp",
                                    "Select one or more cell types",
                                    choices=atlas.obs["Broad_Cell_Annotation"
                                                      ].cat.categories.to_list(),
                                    #selected=None,
                                    multiple=True,
                                    #options=[maxItems=5]
                                ),
                                height="300px"
                            )
                        ),
                        ui.card(
                            ui.card_body(
                                ui.input_selectize(
                                    "selectGeneExp",
                                    "Select one or more genes",
                                    choices=atlas.var_names.to_list(),
                                    multiple=True,
                                    #options=[maxItems=5]
                                ),
                                height="300px"
                            )
                        ),
                        ui.card(
                            ui.card_body(
                                ui.input_action_button(
                                    "searchExp",
                                    "Explore Selected",
                                    width="100%"
                                ),
                                height="100px"
                            )
                        ),
                        #col_widths=12
                    ),
                    width="100%",
                    value="searchBarExp"
                ),
                open=False
            ),
        ),
        ui.navset_pill(
            #ui.panel_conditional(
                ui.nav_panel(
                    "Feature Plot",
                    ui.layout_columns(
                        ui.card(
                            ui.layout_sidebar(
                                ui.sidebar(
                                    "Adjust Parameters",
                                    ui.accordion(
                                        ui.accordion_panel(
                                            "Time Selection",
                                            ui.input_checkbox(
                                                "splitTimeFeat",
                                                "Split by Time",
                                                value=True
                                            )
                                        ),
                                        ui.accordion_panel(
                                            "Color Selection",
                                            ui.input_selectize(
                                                "selectPaletteFeat",
                                                "Select Color Palette",
                                                choices={
                                                    "magma_r":"Dark Magma",
                                                    "magma":"Light Magma",
                                                    "plasma_r":"Dark Plasma",
                                                    "plasma":"Light Plasma",
                                                    "inferno_r":"Dark Inferno",
                                                    "inferno":"Light Inferno",
                                                    "viridis_r":"Dark Viridis",
                                                    "viridis":"Light Viridis"
                                                }
                                            )
                                        ),
                                        open=True,
                                        multiple=True
                                    ),
                                    position="right",
                                    fill=True
                                ),
                                ui.card_body(
                                    ui.output_ui(
                                        "ui_feat_tabs",
                                        fillable=True,
                                        #fill=True
                                    ),
                                    fill=True,
                                    fillable=True
                                )
                            ),
                            #ui.card_footer("selectGeneExp"),
                            ui.card_footer("Selected Genes Feature Plots"),
                            value="featGeneTab",
                            full_screen=True,
                            fill=True,
                            #height="100%"
                        ),
                        row_heights=12
                    ),
                    value="featTab",
                ),
                ui.nav_panel(
                    "Violin Plot",
                    ui.layout_columns(
                        ui.card(
                            ui.layout_sidebar(
                                ui.sidebar(
                                    "Adjust Parameters",
                                    ui.accordion(
                                        ui.accordion_panel(
                                            "Time Selection",
                                            ui.input_checkbox(
                                                "splitTimeVln",
                                                "Split By Time",
                                                value=True
                                            )
                                        ),
                                        #ui.accordion_panel(
                                        #    "Color Selection",
                                        #    ui.output_ui("colorsVln")
                                        #),
                                        open=True,
                                        multiple=True
                                    ),
                                    position="right",
                                    fillable=False
                                ),
                                ui.card_body(
                                    ui.output_ui(
                                        "ui_vln_tabs",
                                        #fillable=True,
                                        fill=True
                                    ),
                                    fill=True,
                                    fillable=True
                                )
                            ),
                            ui.card_footer("Selected Genes Violin Plot"),
                            value="vlnGeneTab",
                            full_screen=True,
                            fill=True
                            #height="100%"
                        ),
                        #row_heights=12
                    ),
                    value="vlnTab"
                ),
                ui.nav_panel(
                    "UMAP",
                    ui.layout_columns(
                        ui.card(
                            ui.layout_sidebar(
                                ui.sidebar(
                                    "Adjust Parameters",
                                    ui.accordion(
                                        ui.accordion_panel(
                                            "Color Selection",
                                            ui.p("Coming soon"),
                                            #ui.output_ui("colorsCell")
                                        ),
                                        open=True,
                                        multiple=True
                                    ),
                                    position="right",
                                    fill=True
                                ),
                                ui.card_body(
                                    output_widget(
                                        "cellUMAP",
                                        fill=True
                                    ),
                                ),
                            ),
                            ui.card_footer("Selected Cells Subtype UMAP"),
                            full_screen=True,
                            fill=True
                        )
                    )
                ),
                ui.nav_panel(
                    "Dot Plot",
                    ui.layout_columns(
                        ui.card(
                            ui.layout_sidebar(
                                ui.sidebar(
                                    "Adjust Parameters",
                                    ui.accordion(
                                        ui.accordion_panel(
                                            "Color Palette",
                                            ui.input_selectize(
                                                "selectCellPaletteDot",
                                                "Select Color Palette",
                                                choices={
                                                    "magma_r":"Dark Magma",
                                                    "magma":"Light Magma",
                                                    "plasma_r":"Dark Plasma",
                                                    "plasma":"Light Plasma",
                                                    "inferno_r":"Dark Inferno",
                                                    "inferno":"Light Inferno",
                                                    "viridis_r":"Dark Viridis",
                                                    "viridis":"Light Viridis"
                                                }
                                            )
                                        )
                                    ),
                                    position="right"
                                ),
                                ui.card_body(
                                    ui.output_plot("cellDot",
                                                   fill=True),
                                    height="600px",
                                    #width="100%"
                                ),
                                ui.card_footer("Selected Genes Dot Plot")
                            )
                        )
                    )
                )
        ),
        value="exploreTab"
    ),
    # Download the full dataset or a subset
    ui.nav_panel("Download",
        ui.h1("Download the dataset",
            align = "center"),
            ui.hr(),
            ui.p("Coming Soon")
    ),
    ui.nav_spacer(),
    # Dark mode toggle button
    ui.nav_control(
        ui.input_dark_mode(id="mode")
    ),
title = "Brain Border Time",
id = "page"
)

# Define server
def server(input, output, session):

    ## set up
    @reactive.effect
    def _():
        ui.update_selectize(
            "selectCellViz",
            choices=atlas.obs["Broad_Cell_Annotation"].cat.categories.to_list(),
            selected=None,
            server=True
        )

    @reactive.effect
    def _():
        ui.update_selectize(
            "selectGeneViz",
            choices=atlas.var_names.to_list(),
            selected=None,
            server=True
        )

    @reactive.effect
    def load_atlas():
        from setup import atlas

    @reactive.effect
    @reactive.event(input.make_light)
    def _():
        ui.update_dark_mode("light")

    @reactive.effect
    @reactive.event(input.make_dark)
    def _():
        ui.update_dark_mode("dark")

    ## Home
    @render.image
    def statUMAP():
        if input.mode() == "dark":
            return get_image("atlas_dark.png",
                             height="120%")
        else:
            return get_image("atlas_light.png",
                             height="120%")

    @render.image
    def timeUMAP():
        if input.mode() == "dark":
            return get_image("atlas_time_split_dark.png",
                             height="40%")
        else:
            return get_image("atlas_time_split_light.png",
                             height="40%")
    
    ## Atlas
    @reactive.effect
    @reactive.event(input.searchViz)
    def _():
        # Switch the active tab to the explore tab
        ui.update_nav_panel("page", target="exploreTab", method="show")
        ui.update_navset("page", selected="exploreTab")
        # Update selected values in explore search bar
        ui.update_selectize(
            "selectCellExp",
            selected=input.selectCellViz()
        )
        ui.update_selectize(
            "selectGeneExp",
            selected=input.selectGeneViz()
        )

    @render_widget
    def hoverUMAP():
        return build_umap(
            adata=atlas,
            mode=input.mode(),
            colors=color_list,
            group=input.selectGroupViz(),
            static=False
        )

    @output
    @render.plot(alt="Top Marker Genes Dot Plot")
    def fullDot():
        #marker_genes = read_csv("../data/top_markers_scanpy.csv")
        build_dot(
            adata=atlas, 
            mode=input.mode(),
            colors=input.selectPaletteDot(),
            group_by="Broad_Cell_Annotation",
            gene_df=marker_genes,
            num_genes=input.numMarkers()
        )

    @render_widget
    def fullHeatmap():
        return build_heatmap(
            scaled_exp=scaled_exp,
            colors=input.selectPaletteHeat(),
            mode=input.mode()
        )
    
    ## Explore
    # list of selected genes and cells
    selectedGenes = reactive.value([])
    selectedCells = reactive.value([])
    # Reactive subsetted atlas object
    atlasCells = reactive.value(atlas)
    # Reactive cell grouping variable
    groupElementsCell = reactive.value("Broad_Cell_Annotation")

    @reactive.effect
    @reactive.event(input.searchExp)
    def _():
        # Update selected values in explore search bar
        ui.update_selectize(
            "selectCellViz",
            selected=input.selectCellExp()
        )
        ui.update_selectize(
            "selectGeneViz",
            selected=input.selectGeneExp()
        )

    @reactive.effect
    @reactive.event(input.searchExp)
    def _():
        # If no genes selected previously just add all selected genes
        if not selectedGenes.get():
            selectedGenes.set(input.selectGeneExp())
        else:
            # Add new gene tabs and remove any required gene tabs
            selectedGenes.set([g for g in input.selectGeneExp() 
                               if g in selectedGenes.get()] + 
                               [g for g in input.selectGeneExp()
                                if g not in selectedGenes.get()])
            
        # If no cells selected previously just add all selected cells    
        if not selectedCells():
            selectedCells.set(input.selectCellExp())
        else:
            # Add new cells and remove any required cells
            selectedCells.set([c for c in input.selectCellExp()
                               if c in selectedCells()] +
                               [c for c in input.selectCellExp()
                                if c not in selectedCells()])
            
        # If any cells selected, filter the atlas to only include those cells
        if selectedCells():
            atlasCells.set(atlas[
                atlas.obs["Broad_Cell_Annotation"].isin(selectedCells())
                ])
        # If no cells selected include the full atlas
        else:
            atlasCells.set(atlas)

        # Determine whether to use broad cell annotations or fine cell 
        # annotations based on the number of fine cell annotations
        if len(atlasCells.get().obs["Fine_Cell_Annotation"].cat.categories
               ) <= len(color_list):
            groupElementsCell.set("Fine_Cell_Annotation")
        else:
            groupElementsCell.set("Broad_Cell_Annotation")
        
        @render.ui
        def ui_feat_tabs():
            for g in selectedGenes():
                nav_panel_server(
                    f"{g}FeatPanel",
                    adata=atlasCells(),
                    mode=input.mode(),
                    colors=input.selectPaletteFeat(),
                    time_split=input.splitTimeFeat(),
                    genePanel=g)
            return ui.navset_tab(
                *[nav_panel_ui(f"{g}FeatPanel", genePanel=g, type="Feat"
                            ) for g in selectedGenes()],
                id="featTabs"
            )

        @render.ui
        def ui_vln_tabs():
            for g in selectedGenes():
                nav_panel_server(
                    f"{g}VlnPanel",
                    adata=atlasCells(),
                    mode=input.mode(),
                    colors=color_list,
                    groupby=groupElementsCell(),
                    time_split=input.splitTimeVln(),
                    genePanel=g)
            return ui.navset_tab(
                *[nav_panel_ui(f"{g}VlnPanel", genePanel=g, type="Vln"
                            ) for g in selectedGenes()],
                id="vlnTabs"
            )

    @reactive.effect
    @reactive.event(input.searchViz)
    def _():
        # If no genes selected previously just add all selected genes
        if not selectedGenes.get():
            selectedGenes.set(input.selectGeneViz())
        else:
            # Add new gene tabs and remove any required gene tabs
            selectedGenes.set([g for g in input.selectGeneViz() 
                               if g in selectedGenes()] + 
                               [g for g in input.selectGeneViz()
                                if g not in selectedGenes()])
            
        # If no cells selected previously just add all selected cells
        if not selectedCells():
            selectedCells.set(input.selectCellViz())
        else:
            # Add new cells and remove and required cells
            selectedCells.set([c for c in input.selectCellViz()
                               if c in selectedCells()] +
                               [c for c in input.selectCellViz()
                                if c not in selectedCells()])
            
        # If any cells selected, filter the atlas to only include those cells
        if selectedCells():
            atlasCells.set(atlas[
                atlas.obs["Broad_Cell_Annotation"].isin(selectedCells())
                ])
        # If no cells selected include the full atlas
        else:
            atlasCells.set(atlas)

        # Determine whether to use broad cell annotations or fine cell 
        # annotations based on the number of fine cell annotations
        if len(atlasCells.get().obs["Fine_Cell_Annotation"].cat.categories
               ) <= len(color_list):
            groupElementsCell.set("Fine_Cell_Annotation")
        else:
            groupElementsCell.set("Broad_Cell_Annotation")

        @render.ui
        def ui_feat_tabs():
            for g in selectedGenes():
                nav_panel_server(
                    f"{g}FeatPanel",
                    adata=atlasCells(),
                    mode=input.mode(),
                    colors=input.selectPaletteFeat(),
                    time_split=input.splitTimeFeat(),
                    genePanel=g)
            return ui.navset_tab(
                *[nav_panel_ui(f"{g}FeatPanel", genePanel=g, type="Feat"
                            ) for g in selectedGenes()],
                id="featTabs"
            )

        @render.ui
        def ui_vln_tabs():
            for g in selectedGenes():
                nav_panel_server(
                    f"{g}VlnPanel",
                    adata=atlasCells(),
                    mode=input.mode(),
                    colors=color_list,
                    groupby=groupElementsCell(),
                    time_split=input.splitTimeVln(),
                    genePanel=g)
            return ui.navset_tab(
                *[nav_panel_ui(f"{g}VlnPanel", genePanel=g, type="Vln"
                            ) for g in selectedGenes()],
                id="vlnTabs"
            )
    
    @render_widget
    def cellUMAP():
        return build_umap(
            adata=atlasCells(),
            mode=input.mode(),
            colors=color_list,
            group=groupElementsCell(),
            static=False
        )
    
    @output
    @render.plot(alt="Selected Cells and Genes Dot Plot")
    def cellDot():
        #marker_genes = read_csv("../data/top_markers_scanpy.csv")
        build_dot(
            adata=atlasCells(), 
            mode=input.mode(),
            colors=input.selectCellPaletteDot(),
            group_by=groupElementsCell(),
            gene_df=list(selectedGenes())
        )

app = App(app_ui, server)
