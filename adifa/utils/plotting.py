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
import plotly.graph_objs as go
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
    colormap="viridis",
    scale_mode="auto",
    scale_max=15,
    scale_min=0,
    tick_no=8,
    scale_log=False,
    plot_covid=True,
    use_premade_info=True,
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
        colormap
    ]  # using premade colormaps e.g. viridis, plasma, inferno, magma, cividis, Reds
    scale = scale_mode  # for the color bar: auto, manual
    scale_lower_value = scale_min
    scale_upper_value = scale_max

    # No input
    if (
        not mode
        or (mode != "gene_expression" and not cat2)
        or (
            cat2
            and adata.obs[cat2].dtype == "category"
            and (not plot_value or not len(plot_value))
        )
    ) and mode != "proportion":
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
                dates = [i.date() for i in dates]

            fig = plot_datetime(cat2, dates, labels, plot_covid)

            return fig.to_json()

        # input is int or float
        elif adata.obs[cat2].dtype in ["float64", "int32", "int64"]:

            fig = go.Figure()

            for index, section in enumerate(adata.obs[cat1].unique()):
                values = (adata.obs[cat2][adata.obs[cat1].isin([section])]).values

                if scale_log == True:
                    values = np.log(values)

                l = len(adata.obs[cat1].unique())
                c = cm.get_cmap(Cmap, l)
                fig.add_trace(
                    go.Violin(x=values, line_color=f"rgb{c(index/l)[:3]}", name=section)
                )

            fig.update_traces(orientation="h", side="positive", width=2, points=False)
            fig.update_layout(
                xaxis_showgrid=False,
                xaxis_zeroline=False,
                title_text=f"Ridgeplot of the continual variable {cat2} across {cat1}",
                xaxis_title=f"{cat2}",
                yaxis_title=f"{cat1}",
                showlegend=False,
            )
            fig.update_layout(
                {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"}
            )

            return fig.to_json()

        # dtype is string/category/object/bool
        else:

            if mode != "proportion" and len(plot_value) > 1:
                adata.obs["combined_annotation"] = adata.obs[cat2].copy().astype(str)
                for value in plot_value:
                    adata.obs.loc[
                        adata.obs[cat2].isin([value]), "combined_annotation"
                    ] = "combined_annotation"
                cat2 = "combined_annotation"
                plot_value = ["combined_annotation"]

            adata.obs[cat1] = adata.obs[cat1].astype("category")
            adata.obs[cat2] = adata.obs[cat2].astype(str).astype("category")

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
                percentage_table_row = (
                    counts_table[plot_value].div(counts_table.sum(axis=1), axis=0) * 100
                )
                values = percentage_table_row[plot_value].values.reshape(-1)

            elif mode == "proportion":
                values = (counts_table["True"] / counts_table.sum(axis=1) * 100).values

    #################################################################################
    # create a color scale on the range of values inputted

    if scale == "auto":
        vmax = max(values) if max(values) > 0 else 1
        vmin = min(values) if min(values) != vmax else 0
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        sm = mpl.cm.ScalarMappable(cmap=Cmap, norm=norm)

    elif scale == "manual":
        norm = mpl.colors.Normalize(vmin=scale_lower_value, vmax=scale_upper_value)
        sm = mpl.cm.ScalarMappable(cmap=Cmap, norm=norm)

    #################################################################################

    fig = go.Figure()

    for i, key in enumerate(adata.uns["polygons"].keys()):
        polygon0 = go.Scatter(
            x=list(*adata.uns["polygons"][key][:, :, 0, 0]),
            y=list(*adata.uns["polygons"][key][:, :, 0, 1]),
            showlegend=False,
            mode="lines",
            fill="toself",
            line=dict(
                color=(
                    "#%02x%02x%02x"
                    % tuple(
                        list(int((255 * x)) for x in list(sm.to_rgba(values[i]))[0:3])
                    )
                ),
                width=1,
            ),
            fillcolor=(
                "#%02x%02x%02x"
                % tuple(list(int((255 * x)) for x in list(sm.to_rgba(values[i]))[0:3]))
            ),
            hoveron="fills",
            # text=text,
            # hoverinfo="text+x+y",
            hoverinfo="x+y",
            # title = title
        )
        fig.add_trace(polygon0)

    if mode not in ["gene_expression", "proportion"] and len(plot_value) == 0:
        title = "Nothing selected to visualise"

    else:
        if mode in ["gene_expression", "proportion"] or (
            mode and plot_value and len(plot_value)
        ):

            if mode == "gene_expression":
                title = f"Mean gene expression of {plot_value[0]} for each section"
                # text = "<br>".join(
                #     [
                #         f"<b>{key}</b>",
                #         "",
                #         f"Section mean gene expression value: {plot_value[0]}",
                #     ]
                # )

            elif mode == "counts":
                title = f"Number of counts for {plot_value[0]}"
                # text = "<br>".join(
                #     [
                #         f"<b>{key}</b>",
                #         "",
                #         f"Number of {plot_value[0]} cells in {key}:    {values[count]} cells",
                #         f"Total number of {plot_value[0]} cells in data:    {sum(values)} cells",
                #     ]
                # )

            elif mode == "percentage_within_sections":
                title = f"Percentage of {plot_value[0]} compared within section"
                # text = "<br>".join(
                #     [
                #         f"<b>{key}</b>",
                #         "",
                #         f"{plot_value[0]} represents {values[count]}% of the cells within {key}",
                #     ]
                # )

            elif mode == "percentage_across_sections":
                title = f"Percentage of {plot_value[0]} across sections"
                # text = "<br>".join(
                #     [
                #         f"<b>{key}</b>",
                #         "",
                #         f"{key} contains {values[count]}% of the cells for {plot_value[0]}",
                #     ]
                # )

            elif mode == "proportion":
                title = f"Percentage of truthful {cat2} values within section"
                # text = "<br>".join(
                #     [
                #         f"<b>{key}</b>",
                #         "",
                #         f"{round(values[count],2)}% (2dp) of the cells within {key} are true for {cat2}",
                #     ]
                # )

    colorbar_trace = go.Scatter(
        x=[None],
        y=[None],
        mode="markers",
        marker=dict(
            colorscale=colormap,
            showscale=True,
            cmin=min(values),
            cmax=max(values),
            colorbar=dict(
                thickness=10, tickmode="auto"
            ),  # tickvals=[min(values), max(values)],
        ),
        hoverinfo="none",
        showlegend=False,
    )
    fig.add_trace(colorbar_trace)
    fig.update_layout(title_text=title)
    fig.update_xaxes(visible=False, fixedrange=True)
    fig.update_yaxes(visible=False, autorange="reversed", fixedrange=True)
    fig.update_layout(
        {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"}
    )

    return fig.to_json()


def plot_datetime(cat2, dates, labels, plot_covid=False):
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

    stems = np.zeros(len(dates))
    stems[::2] = 1  # 0.3
    stems[1::2] = -1  # -0.3

    data = [
        go.Scatter(
            x=dates,
            # y = np.zeros(len(dates)),
            y=stems,
            mode="markers",
            marker=dict(color="red"),
            text=labels,
            hoverinfo="text+x+y",
        )
    ]

    # Use the 'shapes' attribute from the layout to draw the vertical lines
    layout = go.Layout(
        shapes=[
            dict(
                type="line",
                xref="x",
                yref="y",
                x0=i,
                y0=0,
                x1=i,
                y1=stems[dates.index(i)],
                line=dict(color="black", width=1),
            )
            for i in dates
        ],
        title=f"Timeline for {cat2}",
    )

    # Plot the chart
    fig = go.Figure(data, layout)

    fig.add_hline(y=0)

    # add month ticks and labels
    for i in list(
        (pd.date_range(start=start_date, end=end_date, freq="1MS")).map(
            lambda d: str(d.date())
        )
    ):
        x_axis_pos = datetime.date(*list(map(int, i.replace("-", " ").split(" "))))
        fig.add_vline(x=x_axis_pos, y0=0.45, y1=0.55, line_width=1, line_color="black")
        if (x_axis_pos.strftime("%B")[:3]) == "Jan":
            fig.add_annotation(
                dict(
                    font=dict(color="black", size=7),
                    x=x_axis_pos,
                    y=-0.3,
                    showarrow=False,
                    text=x_axis_pos.strftime("%B")[:3],
                    textangle=0,
                    xanchor="left",
                    xref="x",
                    yref="y",
                )
            )
            fig.add_annotation(
                dict(
                    font=dict(color="black", size=10),
                    x=x_axis_pos,
                    y=0.5,
                    showarrow=False,
                    text=x_axis_pos.strftime("%Y"),
                    textangle=0,
                    xanchor="left",
                    xref="x",
                    yref="y",
                )
            )
        else:
            fig.add_annotation(
                dict(
                    font=dict(color="black", size=7),
                    x=x_axis_pos,
                    y=-0.3,
                    showarrow=False,
                    text=x_axis_pos.strftime("%B")[:3],
                    textangle=0,
                    xanchor="left",
                    xref="x",
                    yref="y",
                )
            )

    if plot_covid == True:
        fig.add_vline(
            x=covid_start, y0=0.55, y1=0.65, line_width=1, line_color="purple"
        )
        fig.add_vline(x=covid_end, y0=0.55, y1=0.65, line_width=1, line_color="purple")
        fig.add_shape(
            type="line",
            x0=covid_start,
            y0=0.25,
            x1=covid_end,
            y1=0.25,
            line=dict(color="purple", width=1),
            xref="x",
            yref="y",
        )
        fig.add_annotation(
            dict(
                font=dict(color="purple", size=10),
                x=covid_start,
                y=0.5,
                showarrow=False,
                text="Covid",
                textangle=0,
                xanchor="left",
                xref="x",
                yref="y",
            )
        )

    fig.update_xaxes(visible=False, fixedrange=True)
    fig.update_yaxes(visible=False, fixedrange=True)
    fig.update_layout(
        {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"}
    )
    fig.update_xaxes(type="date", dtick="M1")

    return fig
