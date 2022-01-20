from flask import Blueprint
from flask_restful import Api

from adifa.resources.about import About
from adifa.resources.datasets import Bounds, CellByGeneAggregates, Coordinates, Dataset, Datasets, DiseaseGeneList, Genelist, Genesearch, Labels
from adifa.resources.errors import errors


bp = Blueprint('api_v1', __name__)
api = Api(bp, errors=errors)
   
api.add_resource(About, '/about')
api.add_resource(Bounds, '/bounds')
api.add_resource(CellByGeneAggregates, '/datasets/<id>/cxg')
api.add_resource(Coordinates, '/coordinates')
api.add_resource(Dataset, '/datasets/<id>')
api.add_resource(Datasets, '/datasets')
api.add_resource(DiseaseGeneList, '/datasets/<id>/diseases')
api.add_resource(Genelist, '/datasets/<id>/genes')
api.add_resource(Genesearch, '/datasets/<id>/genesearch')
api.add_resource(Labels, '/labels')
