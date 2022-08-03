from adifa import create_app


def test_config():
    """Test create_app without passing test config."""
    assert not create_app().testing
    assert create_app({"TESTING": True}).testing


def test_hello(client):
    response = client.get("/hello")
    assert response.data == b"Hello, World!"


def test_privacy(client):
    response = client.get("/privacy")
    assert b"Privacy" in response.data


def test_scatterplot(client):
    response = client.get("/dataset/1/scatterplot")
    assert b"Scatterplot in embedding basis" in response.data


def test_matrixplot(client):
    response = client.get("/dataset/1/matrixplot")
    assert b"Heatmap of the mean expression values" in response.data


def test_password(client):
    response = client.get("/dataset/1/password")
    assert b"This dataset is protected" in response.data
