import orjson

from flask import Blueprint, make_response
from flask_restful import Api

from adifa.resources.about import About
from adifa.resources.datasets import Bounds, CellByGeneAggregates, Coordinates, Dataset, Datasets, DiseaseGeneList, Genelist, Genesearch, Labels


class AdifaAPI(Api):
    # This class overrides 'handle_error' method of Api class, to extend
    # global exception handing functionality of flask-restful
    def handle_error(self, e):
        code = getattr(e, 'code', 500)
        message = getattr(e, 'description', 'Internal Server Error')
        errors = getattr(e, 'errors', None)
        to_dict = getattr(e, 'to_dict', None)

        if code == 500:
            logger.exception(e)

        if to_dict:
            data = to_dict()
        else:
            data = {'code': code, 'message': message}

        if isinstance(errors, list):
            data['errors'] = []
            for error in errors:
                data['errors'].append(repr(error))
        elif errors:
            data['errors'] = repr(errors)

        return self.make_response(data, code)


bp = Blueprint('api_v1', __name__)
api = AdifaAPI(bp)
   
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
