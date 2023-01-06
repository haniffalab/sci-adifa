import json

from flask import current_app, flash
from sqlalchemy import exc
import scanpy as sc
import pandas as pd
import base64
from io import BytesIO
from matplotlib.figure import Figure

import numpy as np
import matplotlib as mpl
import cv2

from adifa import models
from adifa.resources.errors import (
    InvalidDatasetIdError,
    DatabaseOperationError,
    DatasetNotExistsError,
)


def get_matrixplot(
    datasetId,
    var_names,
    groupby,
    use_raw=None,
    log=False,
    num_categories=7,
    figsize=None,
    dendrogram=False,
    title=None,
    cmap="viridis",
    colorbar_title="Mean expression\\nin group",
    gene_symbols=None,
    var_group_positions=None,
    var_group_labels=None,
    var_group_rotation=None,
    layer=None,
    standard_scale=None,
    values_df=None,
    swap_axes=False,
    show=None,
    save=None,
    ax=None,
    return_fig=True,
    vmin=None,
    vmax=None,
    vcenter=None,
    norm=None,
    **kwds
):
    if not datasetId > 0:
        raise InvalidDatasetIdError

    try:
        dataset = models.Dataset.query.get(datasetId)
    except exc.SQLAlchemyError as e:
        raise DatabaseOperationError

    try:
        adata = current_app.adata[dataset.filename]
    except (ValueError, AttributeError) as e:
        raise DatasetNotExistsError

    var_intersection = list(set(adata.var.index) & set(var_names))

    if adata.obs[groupby].dtype == "bool":
        adata.obs[groupby] = adata.obs[groupby].astype("str").astype("category")

    plot = sc.pl.matrixplot(
        adata,
        var_intersection,
        groupby,
        use_raw,
        log,
        num_categories,
        figsize,
        dendrogram,
        title,
        cmap,
        colorbar_title,
        gene_symbols,
        var_group_positions,
        var_group_labels,
        var_group_rotation,
        layer,
        standard_scale,
        values_df,
        swap_axes,
        show,
        save,
        ax,
        return_fig,
        vmin,
        vmax,
        vcenter,
        norm,
    )

    if isinstance(plot.categories, pd.IntervalIndex):
        categories = [
            "({}, {}]".format(interval.left, interval.right)
            for interval in plot.categories.values
        ]
    else:
        categories = list(plot.categories)

    output = {
        "categories": categories,
        "var_names": plot.var_names,
        "values_df": json.loads(plot.values_df.to_json()),
        "min_value": str(plot.values_df.min().min()),
        "max_value": str(plot.values_df.max().max()),
        "excluded": list(set(var_names).difference(adata.var.index)),
    }

    return output

def get_spatial_plot(datasetId, plot_value="MACROPHAGE",cat="cell_labels_lvl2", Mode='celltype_percentage_across_sections', color='plasma', scale_mode='auto', scale_max=15, scale_min=0):
    # MACROPHAGE WT1
    if not datasetId > 0:
        raise InvalidDatasetIdError

    try:
        dataset = models.Dataset.query.get(datasetId)
    except exc.SQLAlchemyError as e:
        raise DatabaseOperationError

    try:
        adata = current_app.adata[dataset.filename]
    except (ValueError, AttributeError) as e:
        raise DatasetNotExistsError

    # Calculate tables to hold potential data to plot

    cat1 = 'spatial_location'    # make this the column which the masks relate to i.e. 12 sections
    cat2 = cat  # change to annotations of interest

    adata.obs[cat1] = adata.obs[cat1].astype('category')
    adata.obs[cat2] = adata.obs[cat2].astype('category')
    
    mode = Mode   # celltype_counts, celltype_percentage_within_sections, celltype_percentage_across_sections, gene_expression
    Cmap = mpl.colormaps[color]     # using premade colormaps e.g. viridis, plasma, inferno, magma, cividis, Reds
    scale = scale_mode             # for the color bar: auto, manual
    scale_lower_value = scale_min
    scale_upper_value = scale_max
    

    ###########################################

    # Begin making plot

    if mode == 'celltype_counts':
        # counts table
        cat1_labels = adata.obs[cat1]
        cat2_labels = adata.obs[cat2]
        counts_table = pd.crosstab(cat1_labels,cat2_labels)

        df_of_values = counts_table[[plot_value]]
        values = []
        for col in df_of_values:
            value = list(df_of_values[col].values)
            values.extend(value)
        
    elif mode == 'celltype_percentage_within_sections':
        # counts table
        cat1_labels = adata.obs[cat1]
        cat2_labels = adata.obs[cat2]
        counts_table = pd.crosstab(cat1_labels,cat2_labels)

        # percentage_table_row table:
        percentage_table_row = round((counts_table.T/counts_table.sum(axis=1)).T*100,2)

        df_of_values = percentage_table_row[[plot_value]]
        values = []
        for col in df_of_values:
            value = list(df_of_values[col].values)
            values.extend(value)

    elif mode == 'celltype_percentage_across_sections':
        # counts table
        cat1_labels = adata.obs[cat1]
        cat2_labels = adata.obs[cat2]
        counts_table = pd.crosstab(cat1_labels,cat2_labels)

        # percentage_table_row table:
        percentage_table_column = round((counts_table/counts_table.sum())*100,2)

        df_of_values = percentage_table_column[[plot_value]]
        values = []
        for col in df_of_values:
            value = list(df_of_values[col].values)
            values.extend(value)
            
    elif mode == 'gene_expression':
        pd.set_option('mode.chained_assignment', None)
        gene_expression_table = pd.DataFrame(columns=adata.var_names, index=adata.obs[cat1].cat.categories)
        for clust in gene_expression_table.index.to_list():
            gene_expression_table.loc[clust, :] = adata[adata.obs[cat1].isin([clust]),:].X.mean(0)

        df_of_values = gene_expression_table[[plot_value]]
        df_of_values = df_of_values.T
        values = []
        for col in df_of_values:
            value = list(df_of_values[col].values)
            values.extend(value)
            
    else:
        raise Exception('Mode option not correct. Please use one of the following: manual, celltype_counts, celltype_percentage or gene_expression')
        
    # Load in all masks and orignal base image
    mask1 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_01.jpg')
    mask2 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_02.jpg')
    mask3 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_03.jpg')
    mask4 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_04.jpg')
    mask5 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_05.jpg')
    mask6 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_06.jpg')
    mask7 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_07.jpg')
    mask8 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_08.jpg')
    mask9 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_09.jpg')
    mask10 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_10.jpg')
    mask11 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_11.jpg')
    mask12 = cv2.imread(f'{current_app.root_path}/static/images/Mask_section_12.jpg')
    embryo = cv2.imread(f'{current_app.root_path}/static/images/Embryo_mask_sections_outlined_reverse.jpg')
    embryo_outline = cv2.imread(f'{current_app.root_path}/static/images/Embryo_mask_sections_outlined.jpg')

    # Detect borders for masks
    embryo_grey = cv2.cvtColor(embryo_outline, cv2.COLOR_BGR2GRAY)
    ret, thresh = cv2.threshold(embryo_grey, 50, 255, cv2.THRESH_BINARY)
    contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    cv2.drawContours(embryo, contours, -1, (0,0,0), 3) # black

    # set threshold for every mask
    list_masks = [mask1,mask2,mask3,mask4,mask5,mask6,mask7,mask8,mask9,mask10,mask11,mask12]

    count = 0
    for mask in list_masks:
        _, mask_keep = cv2.threshold(mask, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
        list_masks[count] = mask_keep.copy()
        count+=1

    # create a color scale on the range of values inputted

    if scale == 'auto':
        norm = mpl.colors.Normalize( vmin=min(values) , vmax=max(values) )
        sm = mpl.cm.ScalarMappable(cmap=Cmap, norm=norm)
        
    elif scale == 'manual':
        norm = mpl.colors.Normalize( vmin=scale_lower_value , vmax=scale_upper_value )
        sm = mpl.cm.ScalarMappable(cmap=Cmap, norm=norm)

    else:
        raise Exception('Scale option not correct. Please use either auto or manual')
        
    # Set the color for each mask 
    mask_color = np.copy(embryo)
    count = 0
    for mask in list_masks:
        mask_color[(list_masks[count]==255).all(-1)] = list(int((255*x)) for x in list(sm.to_rgba(values[count]))[0:3])
        count+=1


    # plot the final mask which holds all the other masks and their corresponding colors as well
    fig = Figure(figsize=(4,5))
    ax1, ax2 = fig.subplots(nrows=2, gridspec_kw={"height_ratios":[1, 0.05]})

    im = ax1.imshow(mask_color, interpolation='nearest')
    ax1.set_axis_off()

    cb = fig.colorbar(sm, cax=ax2, orientation='horizontal', pad=0.2)
    cb.ax.tick_params(labelsize=10)
                
    if mode == 'celltype_counts':
        fig.suptitle(f'Number of counts for celltype {plot_value}', fontsize=10, y=0.98, wrap=True)
        cb.set_label("No. cells", fontsize=10)
        
    elif mode == 'celltype_percentage_within_sections':
        fig.suptitle(f'Percentage of celltype {plot_value} compared to other celltypes within section', fontsize=10, y=0.98, wrap=True)
        cb.set_label("Percentage %", fontsize=10)
        
    elif mode == 'celltype_percentage_across_sections':
        fig.suptitle(f'Percentage of  celltype {plot_value} compared across sections', fontsize=10, y=0.98, wrap=True)
        cb.set_label("Percentage %", fontsize=10)
            
    elif mode == 'gene_expression':
        fig.suptitle(f'Mean gene expression of {plot_value} for each section',fontsize=10, y=0.98, wrap=True)
        cb.set_label("Expression", fontsize=10)

    buf = BytesIO()

    fig.savefig(buf, format="png")
    image = base64.b64encode(buf.getbuffer()).decode("ascii")

    return image
