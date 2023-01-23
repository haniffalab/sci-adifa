import json

from flask import request
from flask_restful import Resource

from adifa.utils.plotting import get_matrixplot, get_spatial_plot


class Matrixplot(Resource):
    def get(self, id):
        groupby = request.args.get("groupby", "", type=str)
        var_names = request.args.getlist("var_names")

        return get_matrixplot(id, var_names, groupby)


class Spatial(Resource):
    def get(self, id):
        cat = request.args.get("cat", None, type=str)
        plot_value = request.args.getlist("plot_value[]", None)
        isGene = "gene" in request.args

        if isGene and cat:
            return get_spatial_plot(id, plot_value=cat, mode="gene_expression")
        elif cat:
            return get_spatial_plot(id, cat=cat, plot_value=plot_value)            
        else:
            return get_spatial_plot(id)
