# pylint: disable=missing-docstring
from .misc import split_camel_case


def test_split_camel_case():
    assert split_camel_case("DamBreak") == ["Dam", "Break"]
    assert split_camel_case("CamelCamelCase") == ["Camel", "Camel", "Case"]
    assert split_camel_case("Camel2Camel2Case") == ["Camel2", "Camel2", "Case"]
    assert split_camel_case("getHTTPResponseCode") == [
        "get", "HTTP", "Response", "Code"
    ]
    assert split_camel_case("get2HTTPResponseCode") == [
        "get2", "HTTP", "Response", "Code"
    ]
    assert split_camel_case("HTTPResponseCode") == ["HTTP", "Response", "Code"]
    assert split_camel_case("HTTPResponseCodeXYZ") == [
        "HTTP", "Response", "Code", "XYZ"
    ]
    assert split_camel_case("") == []
