def test_index(client):
    response = client.get("/api/v1/about")
    assert response.data == b"{\"name\":\"AdifaAPI\",\"version\":\"1\"}"
