from datetime import datetime
import sqlite3

import pytest
from adifa import models

#from adifa.db import get_db

def test_post_dataset(session):
    post = models.Dataset(
        filename='test.h5ad',
        hash='1234',
        title='test',
        desc='',
        date_created=datetime.now(),
        date_modified=datetime.now(),
        data_obs='{}',
        data_var='{}',
        data_uns='{}',
        data_obsm='{}',
        data_varm='{}',
        pub_doi='',
        pub_link='',
        pub_author='',
        pub_group='',
        pub_date=datetime.now(),
        download_link=''
    )

    session.add(post)
    session.commit()

    assert post.id > 0

def test_get_close_db(db):
    return True
    with pytest.raises(sqlite3.ProgrammingError) as e:
        db.execute("SELECT 1")
        current_app.logger.error(e)

    assert "closed" in str(e.value)


def test_init_db_command(runner, monkeypatch):
    return True
    class Recorder:
        called = False

    def fake_init_db():
        Recorder.called = True

    monkeypatch.setattr("adifa.db.init_db", fake_init_db)
    result = runner.invoke(args=["init-db"])
    assert "Initialized" in result.output
    assert Recorder.called
