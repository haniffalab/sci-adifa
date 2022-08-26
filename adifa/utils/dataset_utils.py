import os
import hashlib
import gc

from flask import current_app
import zarr

from adifa import db
from adifa import models
from adifa.utils import adata_utils


class AnndataGroup(zarr.Group):
    def __getitem__(self, item):
        if item == "_index":
            current_app.logger.info("_index")
            try:
                return super().__getitem__(self.attrs["_index"])
            except KeyError:
                return super().__getitem__(item)
        else:
            return super().__getitem__(item)


def auto_discover():
    for zarr_dir in [
        x
        for x in os.listdir(current_app.config.get("DATA_PATH"))
        if os.path.isdir(os.path.join(current_app.config.get("DATA_PATH"), x))
    ]:
        if zarr_dir.endswith(".zarr"):  #
            # process file
            current_app.logger.info("Inspecting " + zarr_dir)
            adata = AnndataGroup(
                zarr.open(
                    os.path.join(current_app.config.get("DATA_PATH"), zarr_dir), "r"
                ).store
            )
            annotations = adata_utils.get_annotations(adata)
            # generate hash
            current_app.logger.info("Hashing " + zarr_dir)
            hash = adata["X"].hexdigest()  #

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
                new.filename = zarr_dir
                new.hash = hash
                new.data_obs = annotations.get("obs")
                new.data_obsm = annotations.get("obsm")
                # new.genes_deg = adata_utils.get_degs(adata)
                new.title = zarr_dir
                try:
                    db.session.add(new)
                    db.session.commit()
                except Exception as e:
                    current_app.logger.error("Error: " + str(e))
                else:
                    current_app.logger.info("Adding " + zarr_dir)
            else:
                record.filename = zarr_dir
                record.published = 1
                record.data_obs = annotations.get("obs")
                record.data_obsm = annotations.get("obsm")
                # record.genes_deg = adata_utils.get_degs(adata)
                try:
                    db.session.commit()
                except Exception as e:
                    current_app.logger.error("Error: " + str(e))
                else:
                    current_app.logger.info("Reprocessed " + zarr_dir)
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
        zarr_dir = os.path.join(current_app.config.get("DATA_PATH"), dataset.filename)
        if os.path.isdir(zarr_dir):
            current_app.adata[dataset.filename] = AnndataGroup(
                zarr.open(zarr_dir, "r").store
            )
            dataset.published = 1
            db.session.commit()
            current_app.logger.info("Loading " + dataset.title)
        else:
            dataset.published = 0
            db.session.commit()
            current_app.logger.info("Missing " + dataset.title)
