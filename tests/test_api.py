def test_index(client):
    response = client.get("/api/v1/about")
    assert b"\"name\": \"AdifaAPI\"" in response.data
    assert b"\"version\": \"1\"" in response.data
