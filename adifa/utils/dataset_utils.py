import os
import hashlib

from flask import current_app
import scanpy as sc

from adifa import db
from adifa import models
from adifa.utils import adata_utils


def auto_discover():
    for filename in os.listdir(current_app.config.get("DATA_PATH")):
        if filename.endswith(".h5ad"):
            # process file
            current_app.logger.info("Inspecting " + filename)
            adata = sc.read(os.path.join(current_app.config.get("DATA_PATH"), filename))
            annotations = adata_utils.get_annotations(adata)
            # generate hash
            current_app.logger.info("Hashing " + filename)
            md5_object = hashlib.md5()
            block_size = 128 * md5_object.block_size
            with open(
                os.path.join(current_app.config.get("DATA_PATH"), filename), "rb"
            ) as a_file:
                chunk = a_file.read(block_size)
                while chunk:
                    md5_object.update(chunk)
                    chunk = a_file.read(block_size)

            hash = md5_object.hexdigest()

            # check if exists
            try:
                record = models.Dataset.query.filter_by(hash=hash).first()
            except Exception as e:
                current_app.logger.error("Error: " + str(e))
                return

            exists = bool(record)
            if not exists:
                new = models.Dataset()
                new.published = 1
                new.filename = filename
                new.hash = hash
                new.data_obs = annotations.get("obs")
                new.data_obsm = annotations.get("obsm")
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
                record.filename = filename
                record.published = 1
                record.data_obs = annotations.get("obs")
                record.data_obsm = annotations.get("obsm")
                record.genes_deg = adata_utils.get_degs(adata)
                try:
                    db.session.commit()
                except Exception as e:
                    current_app.logger.error("Error: " + str(e))
                else:
                    current_app.logger.info("Reprocessed " + filename)
        else:
            continue


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
            current_app.adata[dataset.filename] = sc.read(file)
            dataset.published = 1
            db.session.commit()
            current_app.logger.info("Loading " + dataset.title)
        else:
            dataset.published = 0
            db.session.commit()
            current_app.logger.info("Missing " + dataset.title)
