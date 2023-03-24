import json
import textwrap
from itertools import chain
from functools import partial

from flask import current_app, flash
from sqlalchemy import exc
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import cm
import plotly.graph_objs as go
import datetime
from dateutil.relativedelta import relativedelta

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
    #tick_no=8,
    scale_log=False,
    plot_covid=True,
    use_premade_info=True,
    datetime_add_info_col="haniffa_ID",
    mask_set = '12_sections'   #'12_sections' (original),'head_body' (example 2),'multi_head_body' (example with multiple polygons having the same "name" e.g. head_0 and head_1 are 2 polygons but get coloured and have the same values )
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

    cat1 = adata.uns['Mask_selector'][mask_set]  # informs which column in .obs to use relative to the masks set selected
    cat2 = cat  # change to annotations of interest
    cmap = mpl.colormaps[
        colormap
    ]  # using premade colormaps e.g. viridis, plasma, inferno, magma, cividis, Reds
    scale = scale_mode  # for the color bar: auto, manual
    scale_lower_value = scale_min
    scale_upper_value = scale_max

    adata.obs[cat1] = adata.obs[cat1].astype("category")

    if not mode:
        return plot_categorical(adata, mode, cat1, cat2, plot_value, cmap, colormap, mask_set)
    if mode in ["counts", "percentage_within_sections", "percentage_across_sections"]:
        if len(plot_value) > 1:
            adata.obs["combined_annotation"] = adata.obs[cat2].copy().astype(str)
            for value in plot_value:
                adata.obs.loc[
                    adata.obs[cat2].isin([value]), "combined_annotation"
                ] = "combined_annotation"
            cat2 = "combined_annotation"
            plot_value = "combined_annotation"

            adata.obs[cat2] = adata.obs[cat2].astype(str).astype("category")

        return plot_categorical(adata, mode, cat1, cat2, plot_value, cmap, colormap, mask_set)
    elif mode == "gene_expression":
        return plot_gene_expression(adata, cat1, mask_set, plot_value[0], cmap, colormap)
    elif mode == "distribution":
        return plot_distribution(adata, cat1, cat2, cmap, scale_log)
    elif mode == "proportion":
        return plot_proportion(adata, cat1, cat2, cmap, colormap, mask_set, plot_value)
    elif mode == "date":
        return plot_date(
            adata, cat2, use_premade_info, plot_covid, datetime_add_info_col
        )
    else:
        raise


def plot_gene_expression(adata, cat1, mask_set, gene, cmap, colormap):
    df_of_values = (adata.varm[mask_set + "_Sectional_gene_expression"].T)[gene]
    values = list(df_of_values.values)

    title = f"Mean gene expression of <br> {gene} <br> for each section"
    text_template = partial(
        "<br>".join(
            [
                "<b>{key}</b>",
                "",
                "Section mean gene expression value: {value:.3f}",
            ]
        ).format,
    )

    return plot_polygons(adata, cat1, mask_set, values, title, cmap, colormap, text_template)


def plot_proportion(adata, cat1, cat2, cmap, colormap, mask_set, plot_value):

    if not plot_value or not len(plot_value):
        values = [0] * len(adata.obs[cat1])
        title = "Nothing selected to visualise"

        return plot_polygons(adata, cat1, mask_set, values, title, cmap, colormap)
    
    else:
        adata.obs[cat1] = adata.obs[cat1].astype("category")
        adata.obs[cat2] = adata.obs[cat2].astype(str).astype("category")

        counts_table = pd.crosstab(adata.obs[cat1], adata.obs[cat2])

        if len(plot_value) > 1:
            values = [1] * len(adata.obs[cat1])
            title = "Selected to visualise both true and false counts"
            return plot_polygons(adata, cat1, mask_set, values, title, cmap, colormap)

        else:
            values = (counts_table[plot_value[0]] / counts_table.sum(axis=1) * 100).values

            title = f"Percentage of {plot_value[0]} <br> {cat2} <br> values within section"
            text_template = partial(
                "<br>".join(
                    [
                        "<b>{key}</b>",
                        "",
                        "{value:.2f}% of the cells within {key} are true for {cat2}",
                    ]
                ).format,
                cat2=cat2,
            )

            return plot_polygons(adata, cat1, mask_set, values, title, cmap, colormap, text_template)


def plot_categorical(adata, mode, cat1, cat2, plot_value, cmap, colormap, mask_set):

    if not mode or not plot_value or not len(plot_value):
        values = [0] * len(adata.obs[cat1])
        title = "Nothing selected to visualise"

        return plot_polygons(adata, cat1, mask_set, values, title, cmap, colormap)
    elif isinstance(plot_value, list):
        plot_value = plot_value[0]

    assert adata.obs[cat2].dtype == "category"

    counts_table = pd.crosstab(adata.obs[cat1], adata.obs[cat2])

    if mode == "counts":
        df_of_values = counts_table[plot_value]
        values = df_of_values.values.tolist()

        title = f"Number of counts for <br> {plot_value}"
        text_template = partial(
            "<br>".join(
                [
                    "<b>{key}</b>",
                    "",
                    "Number of {plot_value} cells in {key}:\t{value} cells",
                    "Total number of {plot_value} cells in data:\t{sum_values} cells",
                ]
            ).format,
            plot_value=plot_value,
        )

    elif mode == "percentage_across_sections":
        percentage_table_column = round((counts_table / counts_table.sum()) * 100, 2)
        values = percentage_table_column[plot_value].values.tolist()

        title = f"Percentage of <br> {plot_value} <br> across sections"
        text_template = partial(
            "<br>".join(
                [
                    "<b>{key}</b>",
                    "",
                    "{key} contains {value:.2f}% of the cells for {plot_value}",
                ]
            ).format,
            plot_value=plot_value,
        )

    elif mode == "percentage_within_sections":
        percentage_table_row = (
            counts_table[plot_value].div(counts_table.sum(axis=1), axis=0) * 100
        )
        values = percentage_table_row.values.reshape(-1).tolist()

        title = f"Percentage of <br> {plot_value} <br> compared within section"
        text_template = partial(
            "<br>".join(
                [
                    "<b>{key}</b>",
                    "",
                    "{plot_value} represents {value:.2f}% of the cells within {key}",
                ]
            ).format,
            plot_value=plot_value,
        )

    return plot_polygons(adata, cat1, mask_set, values, title, cmap, colormap, text_template)


def plot_polygons(
    adata,
    cat1,
    mask_set,
    values,
    title,
    cmap,
    colormap="viridis",
    text_template=None,
    scale="auto",
    scale_lower_value=0,
    scale_upper_value=100,
):

    print(cat1)
    print(values)


    values_dict = dict(zip(adata.obs[cat1].unique(),values))

    if scale == "auto":
        vmax = max(values) if max(values) > 0 else 1
        vmin = min(values) if min(values) != vmax else 0
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)

    elif scale == "manual":
        norm = mpl.colors.Normalize(vmin=scale_lower_value, vmax=scale_upper_value)
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)

    fig = go.Figure()

    for i, key in enumerate(adata.uns[mask_set + "_polygons"].keys()):
        
        
        polygon0 = go.Scatter(
            x=list(*adata.uns[mask_set + "_polygons"][key][:, :, 0, 0]),
            y=list(*adata.uns[mask_set + "_polygons"][key][:, :, 0, 1]),
            showlegend=False,
            mode="lines",
            fill="toself",
            line=dict(
                color=(
                    "#%02x%02x%02x"
                    % tuple(
                        list(int((255 * x)) for x in list(sm.to_rgba([val for k, val in values_dict.items() if k in key][0]))[0:3])
                    )
                ),
                width=1,
            ),
            fillcolor=(
                "#%02x%02x%02x"
                % tuple(list(int((255 * x)) for x in list(sm.to_rgba([val for k, val in values_dict.items() if k in key][0]))[0:3]))
            ),
            hoveron="fills",
            text=text_template(key=key, value=[val for k, val in values_dict.items() if k in key][0], sum_values=sum(values))
            if text_template
            else None,
            hoverinfo="text+x+y",
            # hoverinfo="x+y",
        )
        fig.add_trace(polygon0)

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
    fig.update_layout(title_text=wrap_text(title), title_x=0.5)
    fig.update_xaxes(visible=False, fixedrange=True)
    fig.update_yaxes(
        visible=False, autorange="reversed", fixedrange=True, scaleanchor="x"
    )
    fig.update_layout(
        {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"}
    )

    return fig.to_json()


def plot_distribution(adata, cat1, cat2, cmap, scale_log=False):
    assert adata.obs[cat2].dtype in ["float64", "int32", "int64"]

    fig = go.Figure()

    for index, section in enumerate(adata.obs[cat1].unique()):
        values = (adata.obs[cat2][adata.obs[cat1].isin([section])]).values

        if scale_log == True:
            values = np.log(values)

        l = len(adata.obs[cat1].unique())
        c = cm.get_cmap(cmap, l)
        fig.add_trace(
            go.Violin(x=values, line_color=f"rgb{c(index/l)[:3]}", name=section)
        )

    fig.update_traces(orientation="h", side="positive", width=2, points=False)
    fig.update_layout(
        xaxis_showgrid=False,
        xaxis_zeroline=False,
        title_text=wrap_text(
            f"Ridgeplot of continual variable <br> {cat2} <br> across {cat1}"
        ),
        title_x=0.5,
        xaxis_title=f"{cat2}",
        yaxis_title=f"{cat1}",
        showlegend=False,
    )
    fig.update_layout(
        {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"}
    )

    return fig.to_json()


def plot_date(
    adata,
    cat2,
    use_premade_info=True,
    plot_covid=False,
    datetime_add_info_col="haniffa_ID",
):

    adata.obs[cat2] = adata.obs[cat2].astype("datetime64[ns]")

    if use_premade_info == True:
        Dates = adata.uns['pre_made_date_col_' + cat2]["dates"]
        dates = []
        for d in Dates:
            dates.append(datetime.date(*list(d)))
        labels = [
            "{0:%d %b %Y}:\n{1}".format(d, l)
            for l, d in zip(adata.uns['pre_made_date_col_' + cat2]["labels"], dates)
        ]
    else:
        df_new = (
            (adata.obs[[datetime_add_info_col, cat2]])
            .reset_index(drop=True)
            .drop_duplicates()
            .reset_index(drop=True)
        ).set_index(datetime_add_info_col)

        date_dict = {}
        for i in df_new.index:
            date_dict[i] = df_new.loc[i][0]

        labels = [
            "{0:%d %b %Y}:\n{1}".format(d, l)
            for l, d in zip(date_dict.keys(), date_dict.values())
        ]
        dates = date_dict.values()
        dates = [i.date() for i in dates]

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

    return fig.to_json()


def wrap_text(text, width=45):
    return "<br>".join(
        list(chain(*[textwrap.wrap(x, width=width) for x in text.split("<br>")]))
    )
