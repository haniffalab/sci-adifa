import os

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'best-password-here'
    SQLALCHEMY_DATABASE_URI = os.environ.get('SQLALCHEMY_DATABASE_URI') or 'sqlite:///../instance/adifa.sqlite'
    SQLALCHEMY_TRACK_MODIFICATIONS = os.environ.get('SQLALCHEMY_TRACK_MODIFICATIONS') or False
    API_VERSION = os.environ.get('API_VERSION') or 1
    API_PREFIX = os.environ.get('API_PREFIX') or '/api'
    API_SERVER = os.environ.get('API_SERVER') or '/'
    DATA_PATH = os.environ.get('DATA_PATH') or './instance/'
    JSONIFY_PRETTYPRINT_REGULAR = False
    JSON_SORT_KEYS = False
