import os

from dotenv import load_dotenv

basedir = os.path.abspath(os.path.dirname(__file__))
load_dotenv(os.path.join(basedir, '.env'))

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY')
    SQLALCHEMY_DATABASE_URI = os.environ.get('SQLALCHEMY_DATABASE_URI')
    SQLALCHEMY_TRACK_MODIFICATIONS = os.environ.get('SQLALCHEMY_TRACK_MODIFICATIONS')
    API_VERSION = os.environ.get('API_VERSION')
    API_PREFIX = os.environ.get('API_PREFIX')
    API_SERVER = os.environ.get('API_SERVER')
    DATA_PATH = os.environ.get('DATA_PATH')
    JSONIFY_PRETTYPRINT_REGULAR = False
    JSON_SORT_KEYS = False
