def test_about(client):
    response = client.get("/api/v1/about")
    assert b"\"name\": \"AdifaAPI\"" in response.data
    assert b"\"version\": \"1\"" in response.data

# def test_bounds(client):
#     response = client.get("/api/v1/bounds?datasetId=1")
#     assert response.status_code == 400

# def test_datasets(client):
#     response = client.get("/api/v1/datasets")
#     assert b"\"datasets\": [" in response.data