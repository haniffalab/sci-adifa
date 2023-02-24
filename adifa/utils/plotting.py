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
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
import cv2
import joypy
import math
import scipy
from scipy import stats
import re
import datetime
from dateutil.relativedelta import relativedelta

mpl.use("agg")

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
    tick_no=8,
    scale_log=False,
    plot_covid=False,
    use_premade_info=False,
    datetime_add_info_col="haniffa_ID",
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

    # Check settings defined in function input

    cat1 = "spatial_location"  # make this the column which the masks relate to i.e. 12 sections
    cat2 = cat  # change to annotations of interest
    cat3 = datetime_add_info_col
    Cmap = mpl.colormaps[
        color
    ]  # using premade colormaps e.g. viridis, plasma, inferno, magma, cividis, Reds
    scale = scale_mode  # for the color bar: auto, manual
    scale_lower_value = scale_min
    scale_upper_value = scale_max

    # No input
    if (
        (not mode
        or (mode != "gene_expression" and not cat2)
        or (
            cat2
            and adata.obs[cat2].dtype == "category"
            and (not plot_value or not len(plot_value))
        ))
        and mode != "proportion"
    ):
        values = [0] * len(adata.obs[cat1])

    elif mode == "gene_expression":
        plot_value = plot_value[0]
        df_of_values = (adata.varm["Sectional_gene_expression"].T)[plot_value]
        values = list(df_of_values.values)

    else:
        # check if string / category / object column in YYYY-MM-DD format and reassign dtype
        if adata.obs[cat2].dtype == "category" and pd.api.types.is_string_dtype(
            adata.obs[cat2].cat.categories.dtype
        ):
            if all(
                adata.obs[cat2].str.match(
                    "^(\d{4})-(0[1-9]|1[0-2]|[1-9])-([1-9]|0[1-9]|[1-2]\d|3[0-1])$"
                )
            ):
                adata.obs[cat2] = adata.obs[cat2].astype("datetime64[ns]")

        # Input is datetime
        if adata.obs[cat2].dtype == "datetime64[ns]":
            if use_premade_info == True:
                Dates = adata.uns[cat2]["dates"]
                dates = []
                for d in Dates:
                    dates.append(datetime.date(*list(d)))
                labels = [
                    "{0:%d %b %Y}:\n{1}".format(d, l)
                    for l, d in zip(adata.uns[cat2]["labels"], dates)
                ]
            else:
                df_new = (
                    (adata.obs[[cat3, cat2]])
                    .reset_index(drop=True)
                    .drop_duplicates()
                    .reset_index(drop=True)
                ).set_index(cat3)

                date_dict = {}
                for i in df_new.index:
                    date_dict[i] = df_new.loc[i][0]

                labels = [
                    "{0:%d %b %Y}:\n{1}".format(d, l)
                    for l, d in zip(date_dict.keys(), date_dict.values())
                ]
                dates = date_dict.values()

            buf = plot_datetime(cat2, dates, labels, plot_covid)

            return base64.b64encode(buf.getvalue()).decode("ascii")

        # input is int or float
        elif adata.obs[cat2].dtype in ["float64", "int32", "int64"]:

            cat_use = "log_scaled_col" if scale_log else cat

            x_lim_range = [
                adata.obs[[cat_use, cat1]][cat_use].min(),
                adata.obs[[cat_use, cat1]][cat_use].max(),
            ]

            ticks_array = []
            ticks = np.linspace(
                float(x_lim_range[0]), float(x_lim_range[1]), num=4, dtype=float
            )
            ticks_array.append(ticks[0])

            if scale_log == True:
                adata.obs["log_scaled_col"] = np.log10(adata.obs[cat2])
                for i in ticks[1:]:
                    ticks_array.append(math.ceil(i))
            else:
                for i in ticks[1:]:
                    ticks_array.append(math.ceil(i / 100) * 100)

            fig, axes = joypy.joyplot(
                data=adata.obs[[cat_use, cat1]],
                by=cat1,
                colormap=Cmap,
                fade=True,
                range_style="group",
                x_range=x_lim_range,
                tails=0,
                xlim="max",
                overlap=0,
                figsize=(4, 5),
                ylabelsize=10,
                xlabelsize=10,
            )

            fig.suptitle(
                f"{cat2} across {cat1}",
                fontsize=10,
                y=1.0
            )

            axes[-1].set_xticks(ticks_array, labels=["{:.1f}".format(x) for x in ticks_array])
            if scale_log == True:
                axes[-1].set_xlabel(
                    f"log10 of {cat2}", fontsize=10, color="black", alpha=1
                )
            else:
                axes[-1].set_xlabel(cat2, fontsize=10, color="black", alpha=1)
            axes[-1].xaxis.set_label_coords(0.5, -0.07)

            patches = [[]] * len(adata.obs[cat1].unique())
            counter = 0
            for i in axes:
                if counter > (len(patches) - 1):
                    break
                current_handles, current_labels = axes[
                    counter
                ].get_legend_handles_labels()
                patches[counter] = (current_handles[0].get_edgecolor()).tolist()[0]
                counter += 1

            counter = 0
            labels = list(adata.obs[cat1].unique())
            legend_patches = []

            gmeans = []
            for i in patches:
                std_ = float(
                    adata.obs[adata.obs[cat1].isin([labels[counter]])][cat2].std()
                )
                kur_ = round(
                    stats.kurtosis(
                        adata.obs[adata.obs[cat1].isin([labels[counter]])][cat2]
                    ),
                    2,
                )
                gmean_ = float(
                    stats.gmean(
                        adata.obs[adata.obs[cat1].isin([labels[counter]])][cat2]
                    )
                )
                gmeans.append(gmean_)
                legend_patches.append(
                    mpatches.Patch(
                        color=i,
                        label="{}:     {:.2f},     {:.2f},     {:.2f}".format(
                            labels[counter], std_, kur_, gmean_
                        ),
                    )
                )
                counter += 1

            fig.legend(
                handles=legend_patches,
                title="Section:         std,        K,        GM",
                loc="lower center",
                bbox_to_anchor=(0.5, -0.59),
                fontsize=10,
                title_fontsize=10,
            )

            for i in range(len(patches)):
                if scale_log == True:
                    X = np.log10(gmeans[i])
                else:
                    X = gmeans[i]
                axes[i].axvline(X, color="red", lw=2, alpha=1, ymax=0.1)

            line1 = Line2D([], [], color='red', marker='|', linestyle='None',
                              markersize=10, markeredgewidth=1.5, label='Geometric mean')
    
            fig.legend(handles=[line1], loc="lower center", bbox_to_anchor=(0.5, -0.65), fontsize=10)
    
            buf = BytesIO()
            fig.savefig(buf, format="png", bbox_inches="tight")
            return base64.b64encode(buf.getvalue()).decode("ascii")
        
        # dtype is string/category/object/bool
        else:

            if mode != "proportion" and len(plot_value) > 1:
                adata.obs["combined_annotation"] = (
                    adata.obs[cat2].copy().astype(str)
                )
                for value in plot_value:
                    adata.obs.loc[
                        adata.obs[cat2].isin([value]), "combined_annotation"
                    ] = "combined_annotation"
                cat2 = "combined_annotation"
                plot_value = ["combined_annotation"]

            adata.obs[cat1] = adata.obs[cat1].astype("category")
            adata.obs[cat2] = adata.obs[cat2].astype(str)
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
            
            elif mode == "proportion":
                values = (counts_table["True"]/counts_table.sum(axis=1)*100).values

    #################################################################################
    # create a color scale on the range of values inputted

    if scale == "auto":
        norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values))
        sm = mpl.cm.ScalarMappable(cmap=Cmap, norm=norm)

    elif scale == "manual":
        norm = mpl.colors.Normalize(vmin=scale_lower_value, vmax=scale_upper_value)
        sm = mpl.cm.ScalarMappable(cmap=Cmap, norm=norm)

    # print(min(values), max(values))
    # print(norm)
    # print(sm.get_clim())

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

    # if (
    #     mode != "gene_expression"
    #     and cat2
    #     and len(plot_value) == len(adata.obs[cat2].unique())
    # ):
    #     fig.suptitle(
    #         f"All categories of {cat2} selected!",
    #         fontsize=10,
    #         y=0.98,
    #         wrap=True,
    #     )
    #     cb.remove()

        # incase we want a legend in future
        #Label = "Number of elements in {c} are {n}".format(c=cat2, n=adata.obs[cat2].count())
        #legend_info = Line2D([], [], color='blue', marker="$\u24D8$", linestyle='None', markersize=10, label=Label)
        #ax1.legend(
        #    title=f"{cat2} information",
        #    handles = [legend_info],
        #    bbox_to_anchor=(1.1, -0.04), # 0.01
        #    labelspacing=2,
        #    fontsize=4,
        #    title_fontsize=8,
        #    prop={"size": 8},
        #)  # 1.4, 1


    # Option 1 - make standard image 
    if mode not in ["gene_expression", "proportion"] and len(plot_value) == 0:
        fig.suptitle(
            "Nothing selected to visualise",
            fontsize=10,
            y=0.98,
            wrap=True,
        )
        cb.remove()

    # Option 2 - load in premade image
    # elif mode!="gene_expression" and cat2 and len(plot_value) == 0:
    #    plt.close()

    #    fig = Figure(figsize=(6, 8))
    #    ax1= fig.subplots(nrows=1)
    #    im = ax1.imshow(adata.uns['Start_img'], interpolation="nearest", aspect='auto')
    #    ax1.set_axis_off()
    #    cb.remove()

    else:
        if mode in ["gene_expression", "proportion"] or (mode and plot_value and len(plot_value)):

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
                    f"Number of counts for {plot_value[0]}",
                    fontsize=10,
                    y=0.98,
                    wrap=True,
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
            
            elif mode == "proportion":
                fig.suptitle(
                    f"Percentage of truthful {cat2} values within section",
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


def plot_datetime(cat, dates, labels, plot_covid=False):
    if plot_covid == True:
        covid_start = datetime.date(2020, 3, 23)
        covid_end = datetime.date(2022, 2, 24)

        if covid_start < min(dates):
            start_date = (covid_start - relativedelta(months=1)).replace(day=1)
        else:
            start_date = (min(dates) - relativedelta(months=1)).replace(day=1)

        if covid_end > max(dates):
            end_date = (covid_end + relativedelta(months=1)).replace(day=1)
        else:
            end_date = (max(dates) + relativedelta(months=1)).replace(day=1)

    else:
        start_date = (min(dates) - relativedelta(months=1)).replace(day=1)
        end_date = (max(dates) + relativedelta(months=1)).replace(day=1)

    fig, ax = plt.subplots(figsize=(4, 5), constrained_layout=True)
    _ = ax.set_ylim(-2, 1.75)
    _ = ax.set_xlim(start_date, end_date)
    _ = ax.axhline(0, xmin=0, xmax=1, c="red", zorder=1)

    _ = ax.get_xaxis().set_major_locator(mdates.MonthLocator(interval=1))
    _ = ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%b %Y"))

    _ = ax.scatter(dates, np.zeros(len(dates)), s=120, c="green", zorder=2)
    _ = ax.scatter(dates, np.zeros(len(dates)), s=30, c="darkgreen", zorder=3)

    label_offsets = np.zeros(len(dates))
    label_offsets[::2] = 1.4  # 0.35
    label_offsets[1::2] = -1.8  # -0.7
    for i, (l, d) in enumerate(zip(labels, dates)):
        _ = ax.text(
            d,
            label_offsets[i],
            l,
            ha="center",
            color="royalblue",
            fontsize=10,
        )

    stems = np.zeros(len(dates))
    stems[::2] = 1  # 0.3
    stems[1::2] = -1  # -0.3
    markerline, stemline, baseline = ax.stem(dates, stems, use_line_collection=True)
    _ = plt.setp(markerline, marker=",", color="green", markersize=5)
    _ = plt.setp(stemline, color="green", linewidth=1.25)

    # hide lines around chart
    for spine in ["left", "top", "right", "bottom"]:
        _ = ax.spines[spine].set_visible(False)

    _ = ax.set_title(
        f"Timeline for {cat}",
        fontsize=10,
        y=1.1,
    )

    # add month ticks and labels
    for i in list(
        (pd.date_range(start=start_date, end=end_date, freq="1MS")).map(
            lambda d: str(d.date())
        )
    ):
        x_axis_pos = datetime.date(*list(map(int, i.replace("-", " ").split(" "))))
        _ = ax.axvline(x_axis_pos, ymin=0.5, ymax=0.57, c="black", zorder=1)
        if (x_axis_pos.strftime("%B")[:3]) == "Jan":
            _ = ax.text(x_axis_pos, -0.3, x_axis_pos.strftime("%B")[:3])
            _ = ax.text(x_axis_pos, 0.5, x_axis_pos.strftime("%Y"))
        else:
            _ = ax.text(x_axis_pos, -0.3, x_axis_pos.strftime("%B")[:3])

    if plot_covid == True:
        # add covid line segment
        x_min, x_max = ax.get_xlim()
        ticks = [(tick - x_min) / (x_max - x_min) for tick in ax.get_xticks()]
        tick_labels = ax.get_xticklabels()
        _ = ax.axvline(covid_start, ymin=0.6, ymax=0.7, c="purple", zorder=1)
        _ = ax.axvline(covid_end, ymin=0.6, ymax=0.7, c="purple", zorder=1)
        _ = ax.plot([covid_start, covid_end], [0.4, 0.4], linestyle="-", color="purple")
        _ = ax.text(covid_start, 0.7, "Covid restrictions", color="purple", fontsize=10)

    # hide tick labels
    _ = ax.set_xticks([])
    _ = ax.set_yticks([])

    buf = BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight")
    buf.seek(0)
    plt.close()

    return buf
