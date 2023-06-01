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
        mask = request.args.get("mask", type=str)
        mode = request.args.get("mode", None, type=str)
        cat = request.args.get("cat", None, type=str)
        plot_value = request.args.getlist("plot_value[]", None)
        colormap = request.args.get("colormap", type=str)
        scale_log = request.args.get(
            "scale_log", False, type=lambda v: v.lower() == "true"
        )

        return get_spatial_plot(
            id,
            mask=mask,
            cat=cat,
            plot_value=plot_value,
            mode=mode,
            colormap=colormap,
            scale_log=scale_log,
        )
