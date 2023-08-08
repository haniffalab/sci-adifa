import os
import re
import hashlib

from flask import current_app
from scipy.sparse import spmatrix
from sqlalchemy import exc
import muon as mu
import scanpy as sc
import numpy as np

from adifa import models
from adifa.resources.errors import (
    InvalidDatasetIdError,
    DatabaseOperationError,
    DatasetNotExistsError,
)


@current_app.template_filter("modality")
def mod_name(mod):
    if mod == "rna":
        return "RNA"
    if mod == "prot":
        return "Protein"
    if mod == "muon":
        return "Muon"
    # atac-seq
    else:
        return mod


def get_annotations(adata):
    annotations = {"obs": {}, "obsm": [], "var": {}}

    switcher = {
        "category": type_category,
        "bool": type_bool,
        "int": type_numeric,
        "float": type_numeric,
        "complex": type_numeric,
    }

    for name in adata.obs:
        # Map numpy dtype to a simple type for switching
        dtype = re.sub(r"[^a-zA-Z]", "", adata.obs[name].dtype.name)
        # Get the function from switcher dictionary
        func = switcher.get(dtype, type_discrete)
        # Define a safe key
        key = hashlib.md5(name.encode("utf-8")).hexdigest()
        # Define obs
        obs = func(adata.obs[name])
        obs["name"] = name
        obs["id"] = key

        if isinstance(adata, mu.MuData) and len(name.split(":")) > 1:
            group, sufix = name.split(":")[0], ":".join(name.split(":")[1:])
            obs["group"], obs["display_name"] = (
                (group, sufix) if group in adata.mod.keys() else ("default", name)
            )
        else:
            obs["group"], obs["display_name"] = ("default", name)

        annotations["obs"][key] = obs

    # remove unwanted obsm arrays
    if isinstance(adata, mu.MuData):
        annotations["obsm"] = [
            value for value in adata.obsm if value not in adata.mod.keys()
        ]
        for mod in adata.mod.keys():
            annotations["obsm"].extend([mod + ":" + value for value in adata[mod].obsm])
    else:
        annotations["obsm"] = adata.obsm_keys()

    annotations["var"] = adata.var_names.tolist()

    return annotations


def get_degs(adata):
    try:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata.var.sort_values(by=["means"], ascending=False)
        df = (
            adata.var[adata.var["highly_variable"] == True]
            .sort_values(by=["means"], ascending=False)
            .head(10)
        )
        return df.index.tolist()
    except Exception:
        return False


def get_bounds(datasetId, obsm):
    if not datasetId > 0:
        raise InvalidDatasetIdError

    try:
        dataset = models.Dataset.query.get(datasetId)
    except exc.SQLAlchemyError:
        raise DatabaseOperationError

    parts = obsm.split(":")
    if len(parts) > 1:
        modality = parts[0]
        obsm = parts[1]
    else:
        modality = False

    try:
        if dataset.filename.endswith(".h5ad"):
            adata = current_app.adata[dataset.filename]
        elif dataset.filename.endswith(".h5mu") and dataset.modality != "muon":
            adata = current_app.adata[dataset.filename][dataset.modality]
        elif dataset.filename.endswith(".h5mu") and modality:
            adata = current_app.adata[dataset.filename][modality]
        elif dataset.filename.endswith(".h5mu"):
            adata = current_app.adata[dataset.filename]
    except (ValueError, AttributeError):
        raise DatasetNotExistsError

    # Normalised [-1,1] @TODO
    adata.obsm[obsm] = (
        2.0 * (adata.obsm[obsm] - np.min(adata.obsm[obsm])) / np.ptp(adata.obsm[obsm])
        - 1
    )

    # Embedded coordinate bounds
    output = {
        "x": {
            "min": adata.obsm[obsm][:, 0].min().item(),
            "max": adata.obsm[obsm][:, 0].max().item(),
        },
        "y": {
            "min": adata.obsm[obsm][:, 1].min().item(),
            "max": adata.obsm[obsm][:, 1].max().item(),
        },
    }

    return output


def get_coordinates(datasetId, obsm):
    if not datasetId > 0:
        raise InvalidDatasetIdError

    try:
        dataset = models.Dataset.query.get(datasetId)
    except exc.SQLAlchemyError:
        raise DatabaseOperationError

    parts = obsm.split(":")
    if len(parts) > 1:
        modality = parts[0]
        obsm = parts[1]
    else:
        modality = False

    try:
        if dataset.filename.endswith(".h5ad"):
            adata = current_app.adata[dataset.filename]
        elif dataset.filename.endswith(".h5mu") and dataset.modality != "muon":
            adata = current_app.adata[dataset.filename][dataset.modality]
        elif dataset.filename.endswith(".h5mu") and modality:
            adata = current_app.adata[dataset.filename][modality]
        elif dataset.filename.endswith(".h5mu"):
            adata = current_app.adata[dataset.filename]
    except (ValueError, AttributeError):
        raise DatasetNotExistsError

    # Normalised [-1,1] @TODO
    adata.obsm[obsm] = (
        2.0 * (adata.obsm[obsm] - np.min(adata.obsm[obsm])) / np.ptp(adata.obsm[obsm])
        - 1
    )

    # True resolution sample generation
    output = []
    for x in adata.obsm[obsm]:
        # output.append(x[:2].tolist())
        output.append([round(num, 4) for num in x[:2].tolist()])

    return output


def get_labels(datasetId, feature="", obs="", modality=""):
    dataset = models.Dataset.query.get(datasetId)
    if dataset.filename.endswith(".h5ad"):
        adata = current_app.adata[dataset.filename]
    elif dataset.filename.endswith(".h5mu") and dataset.modality != "muon":
        adata = current_app.adata[dataset.filename][dataset.modality]
    elif dataset.filename.endswith(".h5mu") and feature:
        adata = current_app.adata[dataset.filename][modality]
    elif dataset.filename.endswith(".h5mu") and obs:
        adata = current_app.adata[dataset.filename]

    if feature:
        try:
            feature_idx = adata.var_names.get_loc(feature)
            output = [
                str(round(float(x), 4))
                for x in (
                    adata.X[:, feature_idx].toarray().reshape(-1)
                    if isinstance(adata.X, spmatrix)
                    else adata.X[:, feature_idx]
                )
            ]
        except KeyError:
            # @todo HANDLE ERROR
            output = [0] * len(adata.obs.index)
        except IndexError:
            # @todo HANDLE ERROR
            output = [0] * len(adata.obs.index)
    elif obs:
        try:
            output = adata.obs[obs].fillna(np.nan).astype(str).tolist()
        except KeyError:
            # @todo HANDLE ERROR
            output = [0] * len(adata.obs.index)
        except IndexError:
            # @todo HANDLE ERROR
            output = [0] * len(adata.obs.index)

    return output


def search_features(datasetId, searchterm, modality):
    dataset = models.Dataset.query.get(datasetId)
    if dataset.filename.endswith(".h5ad"):
        adata = current_app.adata[dataset.filename]
    if dataset.filename.endswith(".h5mu") and dataset.modality == "muon":
        adata = current_app.adata[dataset.filename][modality]
    elif dataset.filename.endswith(".h5mu"):
        adata = current_app.adata[dataset.filename][dataset.modality]

    output = [g for g in adata.var_names if searchterm.lower() in g.lower()]

    return output


def gene_search(datasetId, searchterm):
    dataset = models.Dataset.query.get(datasetId)
    if dataset.filename.endswith(".h5ad") or dataset.modality == "muon":
        adata = current_app.adata[dataset.filename]
    elif dataset.filename.endswith(".h5mu"):
        adata = current_app.adata[dataset.filename][
            dataset.modality
        ]  # adata = current_app.adata
    genes = [g for g in adata.var_names if searchterm in g]

    output = []
    for gene in genes:
        sample = {"name": gene}
        output.append(sample)

    return output


def cat_expr_w_counts(datasetId, cat, gene, func="mean"):
    from numpy import NaN

    dataset = models.Dataset.query.get(datasetId)
    if dataset.filename.endswith(".h5ad") or dataset.modality == "muon":
        adata = current_app.adata[dataset.filename]
    elif dataset.filename.endswith(".h5mu"):
        adata = current_app.adata[dataset.filename][dataset.modality]
    groupall = adata[:, [gene]].to_df().join(adata.obs[cat]).groupby(cat)
    groupexpr = (
        adata[:, [gene]]
        .to_df()
        .replace(float(adata[:, [gene]].X.min()), NaN)
        .join(adata.obs[cat])
        .groupby(cat)
    )

    if func == "mean":
        expr = groupexpr.mean()
    elif func == "median":
        expr = groupexpr.median()

    countpc = (groupexpr.count() * 100 / groupall.count()).astype(int)

    output = [
        {
            "gene": gene,
            "cat": group,
            "expr": float(expr.loc[group, gene]),
            "count": int(countpc.loc[group, gene]),
        }
        for group in groupall.groups.keys()
    ]

    return output


def mode(d):
    from scipy import stats as s

    mode = s.mode(d)
    return int(mode[0])


def type_category(obs):
    categories = [str(i) for i in obs.cat.categories.values.flatten()]

    if len(categories) > 100:
        return {
            "type": "categorical",
            "group": obs.name.split(":")[0]
            if len(obs.name.split(":")) > 1
            else "default",
            "is_truncated": True,
            "values": dict(enumerate(categories[:99], 1)),
        }

    return {
        "type": "categorical",
        "group": obs.name.split(":")[0] if len(obs.name.split(":")) > 1 else "default",
        "is_truncated": False,
        "values": dict(enumerate(categories, 1)),
    }


def type_bool(obs):
    return {
        "type": "categorical",
        "group": obs.name.split(":")[0] if len(obs.name.split(":")) > 1 else "default",
        "values": {0: "True", 1: "False"},
    }


def type_numeric(obs):
    accuracy = 4
    return {
        "type": "continuous",
        "group": obs.name.split(":")[0] if len(obs.name.split(":")) > 1 else "default",
        "min": round(series_min(obs), accuracy),
        "max": round(series_max(obs), accuracy),
        "mean": round(series_mean(obs), accuracy),
        "median": round(series_median(obs), accuracy),
    }


def type_discrete(obs):
    return {"type": "discrete"}


def series_max(s):
    if s.isna().all():
        return 0
    else:
        return s.max().item()


def series_min(s):
    if s.isna().all():
        return 0
    else:
        return s.min().item()


def series_mean(s):
    if s.isna().all():
        return 0
    else:
        return s.mean().item()


def series_median(s):
    if s.isna().all():
        return 0
    else:
        return s.median().item()


def disease_filename():
    return (
        os.path.join(current_app.config.get("DATA_PATH"), "disease.csv")
        if os.path.isfile(
            os.path.join(current_app.config.get("DATA_PATH"), "disease.csv")
        )
        else os.path.join(current_app.root_path, "data", "disease.csv")
    )
