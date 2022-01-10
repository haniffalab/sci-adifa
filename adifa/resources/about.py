from flask_restful import Resource


class About(Resource):
    def get(self):
        return {'name': 'AdifaAPI', 'version': '1'}
