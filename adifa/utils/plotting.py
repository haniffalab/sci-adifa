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
    **kwds,
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


def get_spatial_plot(
    datasetId,
    cat="cell_labels_lvl2",
    plot_value=["MACROPHAGE", "IMMUNE"],
    mode="percentage_within_sections",
    color="viridis",
    scale_mode="auto",
    scale_max=15,
    scale_min=0,
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

    # Calculate tables to hold potential data to plot

    cat1 = "spatial_location"  # make this the column which the masks relate to i.e. 12 sections
    cat2 = cat  # change to annotations of interest

    Cmap = mpl.colormaps[
        color
    ]  # using premade colormaps e.g. viridis, plasma, inferno, magma, cividis, Reds
    scale = scale_mode  # for the color bar: auto, manual
    scale_lower_value = scale_min
    scale_upper_value = scale_max

    ###########################################

    # Begin making plot

    if not mode or not plot_value or not len(plot_value):
        values = [0]*len(adata.obs[cat1])

    elif mode == "gene_expression":
        plot_value = plot_value[0]
        df_of_values = (adata.varm["Sectional_gene_expression"].T)[plot_value]
        values = list(df_of_values.values)

    elif mode in [
        "counts",
        "percentage_within_sections",
        "percentage_across_sections",
    ]:

        if len(plot_value) > 1:
            adata.obs["combined_annotation"] = adata.obs[cat2].copy().astype(str)
            for value in plot_value:
                adata.obs.loc[
                    adata.obs[cat2].isin([value]), "combined_annotation"
                ] = "combined_annotation"
            cat2 = "combined_annotation"
            plot_value = ["combined_annotation"]

        adata.obs[cat1] = adata.obs[cat1].astype("category")
        adata.obs[cat2] = adata.obs[cat2].astype("category")

        counts_table = pd.crosstab(adata.obs[cat1], adata.obs[cat2])

        if mode == "counts":
            df_of_values = counts_table[plot_value]
            values = []
            for col in df_of_values:
                value = list(df_of_values[col].values)
                values.extend(value)

        elif mode == "percentage_across_sections":
            percentage_table_column = round(
                (counts_table / counts_table.sum()) * 100, 2
            )
            df_of_values = percentage_table_column[plot_value]
            values = []
            for col in df_of_values:
                value = list(df_of_values[col].values)
                values.extend(value)

        elif mode == "percentage_within_sections":
            percentage_table_row = round(
                (counts_table.T / counts_table.sum(axis=1)).T * 100, 2
            )
            df_of_values = percentage_table_row[plot_value]
            values = []
            for col in df_of_values:
                value = list(df_of_values[col].values)
                values.extend(value)

    else:
        raise Exception(
            "Mode option not correct. Please use one of the following: manual, counts, percentage or gene_expression"
        )

    #################################################################################
    # create a color scale on the range of values inputted

    if scale == "auto":
        norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values))
        sm = mpl.cm.ScalarMappable(cmap=Cmap, norm=norm)

    elif scale == "manual":
        norm = mpl.colors.Normalize(vmin=scale_lower_value, vmax=scale_upper_value)
        sm = mpl.cm.ScalarMappable(cmap=Cmap, norm=norm)

    else:
        raise Exception("Scale option not correct. Please use either auto or manual")

    #################################################################################

    base_img = np.full(adata.uns["shape"], 255, dtype=np.uint8)

    count = 0
    for key in adata.uns["polygons"].keys():
        cv2.fillPoly(
            base_img,
            pts=tuple([adata.uns["polygons"][key][0]]),
            color=tuple(
                list(int((255 * x)) for x in list(sm.to_rgba(values[count]))[0:3])
            ),
        )
        count += 1

    #################################################################################

    # plot the final mask which holds all the other masks and their corresponding colors as well
    fig = Figure(figsize=(4, 5))
    ax1, ax2 = fig.subplots(nrows=2, gridspec_kw={"height_ratios": [1, 0.05]})

    im = ax1.imshow(base_img, interpolation="nearest")
    ax1.set_axis_off()

    cb = fig.colorbar(sm, cax=ax2, orientation="horizontal", pad=0.2)
    cb.ax.tick_params(labelsize=10)

    #################################################################################

    if mode and plot_value and len(plot_value):

        if mode == "gene_expression":
            fig.suptitle(
                f"Mean gene expression of {plot_value} for each section",
                fontsize=10,
                y=0.98,
                wrap=True,
            )
            cb.set_label("Expression", fontsize=10)

        elif mode == "counts":
                fig.suptitle(
                    f"Number of counts for {plot_value[0]}", fontsize=10, y=0.98, wrap=True
                )
                cb.set_label("No. cells", fontsize=10)

        elif mode == "percentage_within_sections":
            fig.suptitle(
                f"Percentage of {plot_value[0]} within section",
                fontsize=10,
                y=0.98,
                wrap=True,
            )
            cb.set_label("Percentage %", fontsize=10)

        elif mode == "percentage_across_sections":
            fig.suptitle(
                f"Percentage of {plot_value[0]} across sections",
                fontsize=10,
                y=0.98,
                wrap=True,
            )
            cb.set_label("Percentage %", fontsize=10)

    #################################################################################
    buf = BytesIO()

    fig.savefig(buf, format="png")
    image = base64.b64encode(buf.getbuffer()).decode("ascii")

    return image
