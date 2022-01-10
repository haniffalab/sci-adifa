def test_privacy(client):
    response = client.get("/privacy")
    assert b"Privacy" in response.data
