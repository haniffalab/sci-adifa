import os

import pytest

from sqlalchemy import text

from adifa import create_app
from adifa import db as _db


@pytest.fixture(scope="session")
def app():
    """Create and configure a new app instance for each test."""
    # create a temporary file to isolate the database for each test
    basedir = os.path.abspath(os.path.dirname(__file__))
    db_path = os.path.join(basedir, "test.sqlite")
    sql_path = os.path.join(basedir, "data.sql")

    # create the app with common test config
    app = create_app(
        {
            "TESTING": True,
            "SQLALCHEMY_DATABASE_URI": "sqlite:///" + db_path,
            "DATA_PATH": "./test/data",
        }
    )

    # create the database and load test data
    with app.app_context():
        _db.init_app(app)
        _db.create_all()
        with open(sql_path) as f:
            _db.engine.execute(text(f.read()))
        yield app
        _db.session.remove()  # looks like db.session.close() would work as well
        _db.drop_all()

    # remove the temporary database
    os.unlink(db_path)


@pytest.fixture(scope="session")
def client(app):
    """A test client for the app."""
    return app.test_client()


@pytest.fixture(scope="session")
def runner(app):
    """A test runner for the app's Click commands."""
    return app.test_cli_runner()


@pytest.fixture(scope="session")
def db(app):
    _db.app = app
    return _db


@pytest.fixture(scope="function")
def session(db, request):
    """Creates a database session for a test."""
    session = db.create_scoped_session()
    return session
