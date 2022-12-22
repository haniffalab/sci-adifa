import json

from flask import current_app, flash
from sqlalchemy import exc
import scanpy as sc
import pandas as pd
import base64
from io import BytesIO
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
import cv2

from flask import url_for, Response

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

def get_spatial_plot(datasetId):

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
    cat2 = 'cell_labels_lvl2'  # change to annotations of interest

    adata.obs[cat1] = adata.obs[cat1].astype('category')
    adata.obs[cat2] = adata.obs[cat2].astype('category')

    # counts table
    cat1_labels = adata.obs[cat1]
    cat2_labels = adata.obs[cat2]
    counts_table = pd.crosstab(cat1_labels,cat2_labels)

    # percentage_table_row table:
    percentage_table_row = counts_table.T
    percentage_table_row = (round((percentage_table_row/percentage_table_row.sum())*100,2)).T

    # percentage_table_row table:
    percentage_table_column = counts_table
    percentage_table_column = (round((percentage_table_column/percentage_table_column.sum())*100,2))

    gene_expression_table = pd.DataFrame(columns=adata.var_names, index=adata.obs[cat1].cat.categories)                                                                                                 
    for clust in adata.obs[cat1].cat.categories:
        gene_expression_table.loc[clust] = adata[adata.obs[cat1].isin([clust]),:].X.mean(0)

    mode = 'celltype_percentage_across_sections'   # celltype_counts, celltype_percentage_within_sections, celltype_percentage_across_sections, gene_expression
    cmap = plt.cm.viridis      # using premade colormaps e.g. viridis, plasma, inferno, magma, cividis, Reds
    scale = 'auto'             # for the color bar: auto, manual

    #################################################################################
    # Optional depending on choise of mode

    # if chose celltype_counts or celltype_percentage mode
    celltype_to_plot = 'MACROPHAGE'

    # if chose gene_expression mode
    gene_to_plot = 'WT1'   # OR4F5, VCAM1, P2RY12, MLANA, PEML, TYRP1, SHH

    # if chose scale manual
    scale_lower_value = 0
    scale_upper_value = 15

    #################################################################################
        
    if mode == 'celltype_counts':
        df_of_values = counts_table[[celltype_to_plot]]
        values = []
        for col in df_of_values:
            value = list(df_of_values[col].values)
            values.extend(value)
        
    elif mode == 'celltype_percentage_within_sections':
        df_of_values = percentage_table_row[[celltype_to_plot]]
        values = []
        for col in df_of_values:
            value = list(df_of_values[col].values)
            values.extend(value)

    elif mode == 'celltype_percentage_across_sections':
        df_of_values = percentage_table_column[[celltype_to_plot]]
        values = []
        for col in df_of_values:
            value = list(df_of_values[col].values)
            values.extend(value)
            
    elif mode == 'gene_expression':
        df_of_values = gene_expression_table[[gene_to_plot]]
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
    #cv2.drawContours(embryo, contours, -1, (255,255,255), 3) # black

    # set threshold for every mask
    _, mask1 = cv2.threshold(mask1, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask2 = cv2.threshold(mask2, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask3 = cv2.threshold(mask3, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask4 = cv2.threshold(mask4, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask5 = cv2.threshold(mask5, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask6 = cv2.threshold(mask6, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask7 = cv2.threshold(mask7, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask8 = cv2.threshold(mask8, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask9 = cv2.threshold(mask9, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask10 = cv2.threshold(mask10, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask11 = cv2.threshold(mask11, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    _, mask12 = cv2.threshold(mask12, thresh=180, maxval=255, type=cv2.THRESH_BINARY)

    # create a color scale on the range of values inputted

    if scale == 'auto':
        norm = mpl.colors.Normalize( vmin=min(values) , vmax=max(values) )
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        
    elif scale == 'manual':
        norm = mpl.colors.Normalize( vmin=scale_lower_value , vmax=scale_upper_value )
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    else:
        raise Exception('Scale option not correct. Please use either auto or manual')
        
    # Set the color for each mask 
    mask1_color = np.copy(embryo)
    print(mask1_color.shape)
    mask1_color[np.all(mask1==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[0]))[0:3])
    mask2_color = np.copy(mask1_color)
    mask2_color[np.all(mask2==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[1]))[0:3])
    mask3_color = np.copy(mask2_color)
    mask3_color[np.all(mask3==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[2]))[0:3])
    mask4_color = np.copy(mask3_color)
    mask4_color[np.all(mask4==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[3]))[0:3])
    mask5_color = np.copy(mask4_color)
    mask5_color[np.all(mask5==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[4]))[0:3])
    mask6_color = np.copy(mask5_color)
    mask6_color[np.all(mask6==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[5]))[0:3])
    mask7_color = np.copy(mask6_color)
    mask7_color[np.all(mask7==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[6]))[0:3])
    mask8_color = np.copy(mask7_color)
    mask8_color[np.all(mask8==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[7]))[0:3])
    mask9_color = np.copy(mask8_color)
    mask9_color[np.all(mask9==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[8]))[0:3])
    mask10_color = np.copy(mask9_color)
    mask10_color[np.all(mask10==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[9]))[0:3])
    mask11_color = np.copy(mask10_color)
    mask11_color[np.all(mask11==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[10]))[0:3])
    mask12_color = np.copy(mask11_color)
    mask12_color[np.all(mask12==255, -1)] = list(int((255*x)) for x in list(sm.to_rgba(values[11]))[0:3])

    # plot the final mask which holds all the other masks and their corresponding colors as well
    fig = Figure(figsize=(5,3))
    #ax = fig.subplots()
    ax = fig
    
    # plt.figure()

    plt.imshow(mask12_color)

    # ax.axis('off')

    cb = plt.colorbar(ax)
    cb.ax.tick_params(labelsize=20)
        
        
    if mode == 'celltype_counts':
        #plt.title(f'Number of counts for celltype {celltype_to_plot}', fontsize=40, y=1.05)

        if scale == 'manual':
            cb.set_label("No. cells" + " (scale values representitive to manual set upper and lower thresholds)", fontsize=20, rotation=270, labelpad=70)
        elif scale == 'auto':
            cb.set_label("No. cells", fontsize=30, rotation=270, labelpad=70)
        
    elif mode == 'celltype_percentage_within_sections':
        ax.title(f'Percentage of celltype {celltype_to_plot} compared to other celltypes within section', fontsize=40, y=1.05)
        
        if scale == 'manual':
            cb.set_label("Percentage %" + " (scale values representitive to manual set upper and lower thresholds)", fontsize=20, rotation=270, labelpad=70)
        elif scale == 'auto':
            cb.set_label("Percentage %", fontsize=30, rotation=270, labelpad=70)
        
    elif mode == 'celltype_percentage_across_sections':
        ax.title(f'Percentage of  celltype {celltype_to_plot} compared across sections', fontsize=40, y=1.05)
        
        if scale == 'manual':
            cb.set_label("Percentage %" + " (scale values representitive to manual set upper and lower thresholds)", fontsize=20, rotation=270, labelpad=70)
        elif scale == 'auto':
            cb.set_label("Percentage %", fontsize=30, rotation=270, labelpad=70)
            
    elif mode == 'gene_expression':
        ax.title(f'Mean gene expression of {gene_to_plot} for each section',fontsize=40, y=1.05)
        
        if scale == 'manual':
            cb.set_label("Expression" + " (scale values representitive to manual set upper and lower thresholds)", fontsize=20, rotation=270, labelpad=70)
        elif scale == 'auto':
            cb.set_label("Expression", fontsize=30, rotation=270, labelpad=70)

    plt.savefig(f'{current_app.root_path}/static/images/Plotting_elmer.png', dpi=100)

    #fig = Figure()
    # ax = fig.subplots()
    # ax.plot([1,2])

    #img = mpimg.imread(f'{current_app.root_path}/static/images/Plotting_elmer.png')

    buf = BytesIO()

    #FigureCanvas(plt.gcf()).print_png(buf)
    #return Response(buf.getvalue(), mimetype='image/png')
    fig.savefig(buf, format="png")
    image = base64.b64encode(buf.getbuffer()).decode("ascii")
    print(image)

    return image
    # return url_for('static', filename='/images/Plotting_elmer.png')