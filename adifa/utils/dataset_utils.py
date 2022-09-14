import os
import hashlib
import gc

from flask import current_app
import scanpy as sc
import muon as mu

from adifa import db
from adifa import models
from adifa.utils import adata_utils


def generate_hash(file):
    current_app.logger.info("Hashing " + os.path.basename(file))
    md5_object = hashlib.md5()
    block_size = 128 * md5_object.block_size
    with open(file, "rb") as a_file:
        chunk = a_file.read(block_size)
        while chunk:
            md5_object.update(chunk)
            chunk = a_file.read(block_size)

    hash = md5_object.hexdigest()
    return hash


def process_anndata(adata, filename, hash, modality="rna"):
    annotations = adata_utils.get_annotations(adata)

    # check if exists
    try:
        record = models.Dataset.query.filter_by(
            filename=filename, modality=modality
        ).first()
    except Exception as e:
        current_app.logger.error("Error: " + str(e))
        return

    exists = bool(record)
    if not exists:
        new = models.Dataset()
        new.published = 1
        new.filename = filename
        new.modality = modality
        new.hash = hash
        new.data_obs = annotations.get("obs")
        new.data_obsm = annotations.get("obsm")
        new.data_var = annotations.get("var")
        new.genes_deg = adata_utils.get_degs(adata)
        new.title = filename
        try:
            db.session.add(new)
            db.session.commit()
        except Exception as e:
            current_app.logger.error("Error: " + str(e))
        else:
            current_app.logger.info("Adding " + filename)
    else:
        record.hash = hash
        record.published = 1
        record.data_obs = annotations.get("obs")
        record.data_obsm = annotations.get("obsm")
        record.data_var = annotations.get("var")
        record.genes_deg = adata_utils.get_degs(adata)
        try:
            db.session.commit()
        except Exception as e:
            current_app.logger.error("Error: " + str(e))
        else:
            current_app.logger.info("Reprocessed " + filename)


def auto_discover():
    for filename in os.listdir(current_app.config.get("DATA_PATH")):
        file = os.path.join(current_app.config.get("DATA_PATH"), filename)
        if filename.endswith(".h5mu"):
            current_app.logger.info("Inspecting " + filename)
            mudata = mu.read(file)
            hash = generate_hash(file)

            for modality in list(mudata.mod.keys()):
                adata = mudata[modality]
                process_anndata(adata, filename, hash, modality)
            process_anndata(mudata, filename, hash, "muon")
            del mudata

        if filename.endswith(".h5ad"):
            current_app.logger.info("Inspecting " + filename)
            adata = sc.read(file)
            hash = generate_hash(file)
            process_anndata(adata, filename, hash)
            del adata

        else:
            continue
        gc.collect()

    return 5


def load_files():
    # load data files and populate database
    current_app.adata = dict()
    try:
        datasets = models.Dataset.query.all()
    except Exception as e:
        current_app.logger.error(e)
        return

    for dataset in datasets:
        file = os.path.join(current_app.config.get("DATA_PATH"), dataset.filename)
        if os.path.isfile(file):
            if dataset.filename.endswith(".h5ad"):
                current_app.adata[dataset.filename] = sc.read(file)
            if dataset.filename.endswith(".h5mu"):
                if dataset.filename not in current_app.adata:
                    current_app.adata[dataset.filename] = mu.read(file)
            dataset.published = 1
            db.session.commit()
            current_app.logger.info("Loading " + dataset.title)
        else:
            dataset.published = 0
            db.session.commit()
            current_app.logger.info("Missing " + dataset.title)
