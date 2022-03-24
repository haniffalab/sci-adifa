import json

from flask import request
from flask_restful import Resource

from adifa.utils.plotting import get_matrixplot


class Matrixplot(Resource):
    def get(self, id):
        groupby = request.args.get('groupby', '', type=str)
        var_names = json.loads(request.args.get('var_names', '', type=str))

        return get_matrixplot(id, var_names, groupby)

