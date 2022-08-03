from adifa import models


def test_get_dataset(client):
    post = models.Dataset(title="foo")
    dataset = models.Dataset.query.get(1)

    assert dataset.id > 0
