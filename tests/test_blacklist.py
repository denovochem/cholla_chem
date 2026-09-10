from cholla_chem.utils import blacklist as blacklist_module
from cholla_chem.utils.blacklist import filter_blacklisted, get_blacklist_set

# ---------------------------------------------------------------------------
# get_blacklist_set
# ---------------------------------------------------------------------------


def test_get_blacklist_set_returns_frozenset():
    """get_blacklist_set should return a frozenset."""
    result = get_blacklist_set()
    assert isinstance(result, frozenset)


def test_get_blacklist_set_is_nonempty():
    """The blacklist should contain at least some known bad names."""
    result = get_blacklist_set()
    assert len(result) > 0


def test_get_blacklist_set_is_lowercase():
    """All entries in the blacklist should be lowercase for case-insensitive matching."""
    result = get_blacklist_set()
    for name in result:
        assert name == name.lower(), f"Blacklist entry '{name}' is not lowercase"


def test_get_blacklist_set_contains_known_entries():
    """Known blacklisted names should be present (lowercase)."""
    result = get_blacklist_set()
    assert "example 2" in result
    assert "compound 9" in result
    assert "resin" in result


def test_get_blacklist_set_is_cached():
    """Repeated calls should return the same cached object."""
    first = get_blacklist_set()
    second = get_blacklist_set()
    assert first is second


def test_get_blacklist_set_graceful_fallback(monkeypatch):
    """If the file cannot be loaded, an empty frozenset should be returned."""

    class FakeFiles:
        def joinpath(self, *args):
            raise FileNotFoundError("Simulated missing file")

    monkeypatch.setattr(
        blacklist_module.resources, "files", lambda pkg: FakeFiles(), raising=True
    )
    monkeypatch.setattr(blacklist_module, "_BLACKLIST_CACHE", None)

    result = get_blacklist_set()
    assert result == frozenset()


# ---------------------------------------------------------------------------
# filter_blacklisted
# ---------------------------------------------------------------------------


def test_filter_blacklisted_removes_matching():
    """Known blacklisted names should be removed."""
    names = ["ethanol", "Example 2", "water", "Compound 9"]
    result = filter_blacklisted(names)
    assert "ethanol" in result
    assert "water" in result
    assert "Example 2" not in result
    assert "Compound 9" not in result


def test_filter_blacklisted_case_insensitive():
    """Case variants of blacklisted names should be filtered."""
    names = ["COMPOUND 9", "compound 9", "Compound 9", "ethanol"]
    result = filter_blacklisted(names)
    assert result == ["ethanol"]


def test_filter_blacklisted_preserves_order():
    """Non-blacklisted names should preserve their original order."""
    names = ["ethanol", "Example 2", "water", "Compound 9", "methanol"]
    result = filter_blacklisted(names)
    assert result == ["ethanol", "water", "methanol"]


def test_filter_blacklisted_preserves_casing():
    """Original casing of non-blacklisted names should be preserved."""
    names = ["Ethanol", "Water"]
    result = filter_blacklisted(names)
    assert result == ["Ethanol", "Water"]


def test_filter_blacklisted_empty_input():
    """Empty input should return empty list."""
    assert filter_blacklisted([]) == []


def test_filter_blacklisted_all_blacklisted():
    """If all names are blacklisted, result should be empty."""
    names = ["Example 2", "Compound 9", "resin"]
    result = filter_blacklisted(names)
    assert result == []


def test_filter_blacklisted_no_blacklisted():
    """If no names are blacklisted, all should be returned."""
    names = ["ethanol", "water", "methanol"]
    result = filter_blacklisted(names)
    assert result == ["ethanol", "water", "methanol"]


def test_filter_blacklisted_with_custom_blacklist():
    """A custom blacklist set can be provided."""
    custom = frozenset({"ethanol", "water"})
    names = ["ethanol", "methanol", "water"]
    result = filter_blacklisted(names, custom)
    assert result == ["methanol"]


def test_filter_blacklisted_empty_blacklist_returns_all():
    """An empty blacklist should return all names unchanged."""
    names = ["ethanol", "water"]
    result = filter_blacklisted(names, frozenset())
    assert result == ["ethanol", "water"]
