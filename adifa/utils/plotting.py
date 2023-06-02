from __future__ import annotations
from typing import Union
import json
import textwrap
from itertools import chain
from functools import partial

from flask import current_app
from sqlalchemy import exc
import zarr
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import cm
import plotly.graph_objs as go
import datetime
from dateutil.relativedelta import relativedelta

from adifa import models
from adifa.utils.adata_utils import (
    parse_data,
    parse_array,
    parse_group,
    get_group_index,
)
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

    vars = get_group_index(adata["var"])[:]
    var_intersection = [
        x for x in var_names if x in set(vars)
    ]  # preserve var_names order

    sorter = np.argsort(vars)
    var_indx = sorter[np.searchsorted(vars, var_intersection, sorter=sorter)]

    if type(adata["obs"][groupby]).__name__ == "Group":
        obs_df = pd.DataFrame(
            [str(x) for x in list(parse_group(adata["obs"][groupby]))],
            index=get_group_index(adata["obs"])[:],
            columns=[groupby],
        )
    else:
        obs_arr = parse_array(adata["obs"][groupby], adata)
        obs_df = pd.DataFrame(
            [int(x) for x in list(obs_arr)]
            if obs_arr.dtype in ["int"]
            else [float(x) for x in list(obs_arr)]
            if obs_arr.dtype in ["float"]
            else [str(x) for x in list(obs_arr)],
            index=get_group_index(adata["obs"])[:],
            columns=[groupby],
        )
    tempdata = sc.AnnData(adata["X"].oindex[:, var_indx])

    tempdata.var_names = var_intersection
    tempdata.obs = obs_df

    plot = sc.pl.matrixplot(
        tempdata,
        var_intersection,  # multiple var
        groupby,  # single obs
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
        "excluded": list(set(var_names).difference(vars)),
    }

    return output


def get_spatial_plot(
    datasetId: int,
    mask: str,
    mode: str = None,
    cat: str = None,
    plot_value: list[str] = [],
    colormap: str = "viridis",
    # scale_mode: str = "auto",
    # scale_max: int = 15,
    # scale_min: int = 0,
    # tick_no: int = 8,
    scale_log: bool = False,
    plot_covid: bool = False,
    use_premade_info: bool = True,
    datetime_add_info_col: str = "haniffa_ID",
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

    cat1 = adata.uns["masks"][mask]["obs"][()]
    cat2 = cat  # change to annotations of interest
    cmap = mpl.colormaps[
        colormap
    ]  # using premade colormaps e.g. viridis, plasma, inferno, magma, cividis, Reds
    # scale = scale_mode  # for the color bar: auto, manual
    # scale_lower_value = scale_min
    # scale_upper_value = scale_max

    obs_cat1 = pd.Categorical(parse_data(adata.obs[cat1]))
    obs_cat2 = parse_data(adata.obs[cat2], adata)

    if not mode:
        return plot_categorical(
            adata=adata,
            mode=mode,
            obs_cat1=obs_cat1,
            obs_cat2=obs_cat2,
            plot_value=plot_value,
            cmap=cmap,
            colormap=colormap,
            mask=mask,
        )
    if mode in ["counts", "percentage_within_sections", "percentage_across_sections"]:
        # cat2 is categorical
        if len(plot_value) > 1:
            obs_cat2_df = pd.DataFrame(
                obs_cat2.astype(str),
                index=get_group_index(adata.obs),
                columns=["combined_annotation"],
            )
            for value in plot_value:
                obs_cat2_df.loc[
                    obs_cat2.isin([value]), "combined_annotation"
                ] = "combined_annotation"
            obs_cat2 = pd.Categorical(obs_cat2_df["combined_annotation"])
            plot_value = "combined_annotation"

        return plot_categorical(
            adata=adata,
            mode=mode,
            obs_cat1=obs_cat1,
            obs_cat2=obs_cat2,
            plot_value=plot_value,
            cmap=cmap,
            colormap=colormap,
            mask=mask,
        )
    elif mode == "gene_expression":
        return plot_gene_expression(
            adata=adata,
            obs_cat1=obs_cat1,
            mask=mask,
            gene=plot_value[0],
            cmap=cmap,
            colormap=colormap,
        )
    elif mode == "distribution":
        return plot_distribution(
            cat1=cat1,
            cat2=cat2,
            obs_cat1=obs_cat1,
            obs_cat2=obs_cat2,
            cmap=cmap,
            scale_log=scale_log,
        )
    elif mode in ["proportion_within_sections", "proportion_across_sections"]:
        return plot_proportion(
            adata=adata,
            mode=mode,
            cat2=cat2,
            obs_cat1=obs_cat1,
            obs_cat2=obs_cat2,
            cmap=cmap,
            colormap=colormap,
            mask=mask,
            plot_value=plot_value[0],
        )
    elif mode == "date":
        return plot_date(
            adata=adata,
            cat2=cat2,
            obs_cat2=obs_cat2,
            use_premade_info=use_premade_info,
            plot_covid=plot_covid,
            datetime_add_info_col=datetime_add_info_col,
        )
    else:
        raise


def plot_gene_expression(
    adata: sc.AnnData,
    obs_cat1: pd.Categorical,
    mask: str,
    gene: str,
    cmap,
    colormap: str,
):
    group_df = parse_group(adata.varm[adata.uns["masks"][mask]["varm"][()]])
    df_of_values = (group_df.T)[gene]
    values = list(df_of_values.values)

    # from sklearn.preprocessing import MinMaxScaler
    # scaler = MinMaxScaler()
    # values = np.concatenate(scaler.fit_transform(np.array(values).reshape(-1, 1))).tolist()

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

    return plot_polygons(
        adata, obs_cat1, mask, values, title, cmap, colormap, text_template
    )


def plot_proportion(
    adata: sc.AnnData,
    mode: str,
    cat2: str,
    obs_cat1: pd.Categorical,
    obs_cat2: pd.Categorical,
    cmap,
    colormap: str,
    mask: str,
    plot_value: str,
):

    if not plot_value:
        values = [0] * len(obs_cat1)
        title = "Nothing selected to visualise"

        return plot_polygons(adata, obs_cat1, mask, values, title, cmap, colormap)

    else:
        obs_cat2 = pd.Categorical(obs_cat2.astype("str"))

        counts_table = pd.crosstab(obs_cat1, obs_cat2)

        if mode == "proportion_within_sections":

            values = (counts_table[plot_value] / counts_table.sum(axis=1) * 100).values

            title = f"Percentage of {plot_value} <br> {cat2} <br> values within section"
            text_template = partial(
                "<br>".join(
                    [
                        "<b>{key}</b>",
                        "",
                        "{value:.2f}% of the cells within {key} are{b} {cat2}",
                    ]
                ).format,
                cat2=cat2,
                b="" if plot_value == "True" else " not",
            )

        elif mode == "proportion_across_sections":
            values = (counts_table / counts_table.sum() * 100)[plot_value].values

            title = (
                f"Percentage of {plot_value} <br> {cat2} <br> values across sections"
            )
            text_template = partial(
                "<br>".join(
                    [
                        "<b>{key}</b>",
                        "",
                        "{key} contains {value:.2f}% of the cells which are{b} {cat2}",
                    ]
                ).format,
                cat2=cat2,
                b="" if plot_value == "True" else " not",
            )

        return plot_polygons(
            adata, obs_cat1, mask, list(values), title, cmap, colormap, text_template
        )


def plot_categorical(
    adata: sc.AnnData,
    mode: str,
    obs_cat1: pd.Categorical,
    obs_cat2: pd.Categorical,
    plot_value: Union[list[str], str],
    cmap,
    colormap: str,
    mask: str,
):

    if not mode or not plot_value or not len(plot_value):
        values = [0] * len(obs_cat2)
        title = "Nothing selected to visualise"

        return plot_polygons(adata, obs_cat1, mask, values, title, cmap, colormap)
    elif isinstance(plot_value, list):
        plot_value = plot_value[0]

    assert obs_cat2.dtype == "category"

    counts_table = pd.crosstab(obs_cat1, obs_cat2)

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

    return plot_polygons(
        adata, obs_cat1, mask, values, title, cmap, colormap, text_template
    )


def plot_polygons(
    adata: sc.AnnData,
    obs_cat1: pd.Categorical,
    mask: str,
    values: list[float],
    title: str,
    cmap,
    colormap: str = "viridis",
    text_template=None,
    scale: str = "auto",
    scale_lower_value: int = 0,
    scale_upper_value: int = 100,
):

    values_dict = dict(zip(obs_cat1.unique(), values))

    if scale == "auto":
        vmax = max(values) if max(values) > 0 else 1
        vmin = min(values) if min(values) != vmax else 0
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)

    elif scale == "manual":
        norm = mpl.colors.Normalize(vmin=scale_lower_value, vmax=scale_upper_value)
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)

    fig = go.Figure()

    for i, key in enumerate(adata.uns["masks"][mask]["polygons"].keys()):

        polygon0 = go.Scatter(
            x=list(*adata.uns["masks"][mask]["polygons"][key][:, :, 0, 0]),
            y=list(*adata.uns["masks"][mask]["polygons"][key][:, :, 0, 1]),
            showlegend=False,
            mode="lines",
            fill="toself",
            line=dict(
                color=(
                    "#%02x%02x%02x"
                    % tuple(
                        list(
                            int((255 * x))
                            for x in list(
                                sm.to_rgba(
                                    [val for k, val in values_dict.items() if k in key][
                                        0
                                    ]
                                )
                            )[0:3]
                        )
                    )
                ),
                width=1,
            ),
            fillcolor=(
                "#%02x%02x%02x"
                % tuple(
                    list(
                        int((255 * x))
                        for x in list(
                            sm.to_rgba(
                                [val for k, val in values_dict.items() if k in key][0]
                            )
                        )[0:3]
                    )
                )
            ),
            hoveron="fills",
            text=wrap_text(
                text_template(
                    key=key,
                    value=[val for k, val in values_dict.items() if k in key][0],
                    sum_values=sum(values),
                )
            )
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


def plot_distribution(
    cat1: str,
    cat2: str,
    obs_cat1: pd.Categorical,
    obs_cat2: np.ndarray,
    cmap,
    scale_log: bool = False,
):
    assert obs_cat2.dtype in ["float64", "int32", "int64"]

    fig = go.Figure()

    for index, section in enumerate(obs_cat1.unique()):
        values = obs_cat2[obs_cat1.isin([section])]

        # hover_info_text = "<br>".join(["<b>{section}</b>",
        #    "",
        #    f"Min value: {min(values)}",
        #    f"Max value: {max(values)}",
        #    f"Mean value: {np.mean(values)}",
        #    f"Median value: {np.median(values)}"])

        # hover_info_text = partial(
        #    "<br>".join(
        #        [
        #            "<b>{section}</b>",
        #            "",
        #            f"Min value: {min(values)}",
        #            f"Max value: {max(values)}",
        #            f"Mean value: {np.mean(values)}",
        #            f"Median value: {np.median(values)}",
        #        ]
        #    )
        # )

        if scale_log == True:
            values = np.log(np.square(values))
            x_title = f"{cat2} (log(x^2))"
        else:
            x_title = f"{cat2}"

        l = len(obs_cat1.unique())
        c = cm.get_cmap(cmap, l)
        fig.add_trace(
            go.Violin(
                x=values,
                line_color=f"rgb{c(index/l)[:3]}",
                name=section,
                hoverinfo="none",  # seems to be an issue with plotly violin that it shows all stats if use x or same number boxes as y and won't update text properly
                # hovertemplate=hover_info_text,
                # hovertext = hover_info_text
            )
        )

    fig.update_traces(orientation="h", side="positive", width=2, points=False)
    fig.update_layout(
        xaxis_showgrid=False,
        xaxis_zeroline=False,
        title_text=wrap_text(
            f"Ridgeplot of continual variable <br> {cat2} <br> across {cat1}"
        ),
        title_x=0.5,
        xaxis_title=x_title,
        yaxis_title=f"{cat1}",
        showlegend=True,
        legend=dict(font=dict(size=8, color="black")),
    )
    fig.update_layout(
        {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"}
    )

    return fig.to_json()


def plot_date(
    adata: sc.AnnData,
    cat2: str,
    obs_cat2: pd.Categorical,
    use_premade_info: bool = True,
    plot_covid: bool = False,
    datetime_add_info_col: str = "haniffa_ID",
):

    obs_cat2 = obs_cat2.astype("datetime64[ns]")

    if use_premade_info == True:
        Dates = list(adata.uns["premade_date_information"][cat2]["dates"][:])
        dates = []
        for d in Dates:
            dates.append(datetime.date(*list(d)))
        labels = [
            "{0:%d %b %Y}:\n{1}".format(d, l)
            for l, d in zip(
                list(adata.uns["premade_date_information"][cat2]["labels"][:]), dates
            )
        ]
    else:
        df_new = (
            pd.DataFrame(
                {
                    datetime_add_info_col: parse_data(adata.obs[datetime_add_info_col]),
                    cat2: obs_cat2,
                }
            )
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

    data = [
        go.Scatter(
            x=dates,
            y=np.zeros(len(dates)),
            # y=stems,
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
                y1=-1,
                line=dict(color="black", width=1),
            )
            for i in dates
        ],
    )

    # Plot the chart
    fig = go.Figure(data, layout)

    fig.update_layout(
        title={
            "text": f"Timeline for {cat2}",
            "y": 0.9,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        }
    )

    if plot_covid:
        covid_start = datetime.date(2020, 3, 23)
        covid_end = datetime.date(2022, 2, 24)

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

    fig.update_xaxes(
        showline=True,
        linewidth=2,
        linecolor="black",
        visible=True,
        fixedrange=False,
        autorange=True,
        rangeslider=dict(autorange=True, thickness=0.3, bgcolor="#e4f7fe"),
    )
    fig.update_yaxes(visible=False, fixedrange=True)
    fig.update_layout(
        {"plot_bgcolor": "rgba(0,0,0,0)", "paper_bgcolor": "rgba(0,0,0,0)"},
        autosize=True,
    )
    fig.update_xaxes(type="date")

    return fig.to_json()


def wrap_text(text: str, width: int = 45):
    return "<br>".join(
        list(chain(*[textwrap.wrap(x, width=width) for x in text.split("<br>")]))
    )
