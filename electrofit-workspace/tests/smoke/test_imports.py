import importlib.util as u

def test_imports():
    assert u.find_spec("electrofit")
    assert u.find_spec("electrofit.cli.app")
    # legacy workflows package removed; no further assertions