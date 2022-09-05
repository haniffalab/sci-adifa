import json

from flask import current_app, flash
from sqlalchemy import exc
import scanpy as sc
import pandas as pd
import numpy as np

from adifa import models
from adifa.utils.adata_utils import parse_array, parse_group, get_group_index
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
        obs_arr = parse_array(adata, adata["obs"][groupby])
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
