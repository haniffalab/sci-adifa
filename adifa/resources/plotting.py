import json

from flask import request
from flask_restful import Resource

from adifa.utils.plotting import get_matrixplot


class Matrixplot(Resource):
    def get(self, id):
        modality = request.args.get("modality", "rna", type=str)
        groupby = request.args.get("groupby", "", type=str)
        var_names = request.args.getlist("var_names")

        return get_matrixplot(id, modality, var_names, groupby)
