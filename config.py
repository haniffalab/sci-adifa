import os


class Config(object):
    SECRET_KEY = os.environ.get("SECRET_KEY") or "best-password-here"
    SQLALCHEMY_DATABASE_URI = (
        os.environ.get("SQLALCHEMY_DATABASE_URI")
        or "sqlite:///../instance/adifa.sqlite"
    )
    SQLALCHEMY_TRACK_MODIFICATIONS = (
        os.environ.get("SQLALCHEMY_TRACK_MODIFICATIONS") or False
    )
    API_VERSION = os.environ.get("API_VERSION") or 1
    API_PREFIX = os.environ.get("API_PREFIX") or "/api"
    API_SERVER = os.environ.get("API_SERVER") or "/"
    DATA_PATH = os.environ.get("DATA_PATH") or "./instance"
    HOME_URL = os.environ.get("HOME_URL") or False
    KEEP_OBS_ORDER = os.environ.get("KEEP_OBS_ORDER", "False").lower() in ("1", "true")
    JSONIFY_PRETTYPRINT_REGULAR = False
    JSON_SORT_KEYS = False
