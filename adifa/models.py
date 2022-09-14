from email.policy import default
import os
from datetime import datetime

from flask import current_app

from adifa import db


class Dataset(db.Model):
    __tablename__ = "datasets"
    id = db.Column(db.Integer, primary_key=True)
    published = db.Column(db.Boolean, default=False, nullable=False)
    filename = db.Column(db.String(120), nullable=False)
    hash = db.Column(db.String(120), nullable=False)
    title = db.Column(db.String(120), nullable=False)
    desc = db.Column(db.String(500), nullable=True)
    modality = db.Column(db.String(120), nullable=False)
    date_created = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    date_modified = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    data_obs = db.Column(db.JSON, default={})
    data_var = db.Column(db.JSON, default={})
    data_uns = db.Column(db.JSON, default={})
    data_obsm = db.Column(db.JSON, default={})
    data_varm = db.Column(db.JSON, default={})
    pub_doi = db.Column(db.String(120), nullable=True)
    pub_link = db.Column(db.String(512), nullable=True)
    pub_author = db.Column(db.String(120), nullable=True)
    pub_group = db.Column(db.String(120), nullable=True)
    pub_date = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
    download_link = db.Column(db.String(120), nullable=True)
    password = db.Column(db.String(120), nullable=True)
    genes_deg = db.Column(db.JSON, default={})
    genes_curated = db.Column(db.JSON, default={})

    def __repr__(self):
        return f"Dataset('{self.filename}')"

    @property
    def serialize(self):
        """Return object data in easily serializable format"""
        return {
            "id": self.id,
            "filename": self.filename,
            "size": os.path.getsize(
                os.path.join(current_app.config.get("DATA_PATH"), self.filename)
            ),
            "modality": self.modality,
            "data_obs": self.data_obs,
            "data_obsm": self.data_obsm,
            "data_var": self.data_var,
            "title": self.title,
        }
