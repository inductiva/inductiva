"""Test the format_utils module."""

import datetime
from unittest.mock import patch
from pytest import mark
import pytest

from inductiva.utils import format_utils


@mark.parametrize("time_input, expected_time", [(60, "0:01:00"),
                                                (40, "0:00:40"),
                                                (3600, "1:00:00"),
                                                (3660, "1:01:00"),
                                                (3615, "1:00:15"),
                                                (645, "0:10:45")])
def test_seconds_formatter__valid_time__formatted_time(time_input,
                                                       expected_time):
    """Check that the seconds formatter returns the correct string."""
    assert format_utils.seconds_formatter(time_input) == expected_time


tabular_dict = {"A": ["aa", "aaa"], "B": [1, 11], "C": ["cc", "ccc"]}

tabular_rows = [["aa", 1, "cc"], ["aaa", 11, "ccc"]]

tabular_headers = ["A", "Bb", "Cc"]

tabular_formatters = {"A": [lambda x: x.upper()], "Z": [lambda x: x.lower()]}


def test_get_tabular_data__input_dict_no_headers_no_formatters__returns_dict():
    """Check if get_tabular_data function returns the correct dictionary
    when passed a dictionary and no headers or formatters."""
    assert format_utils.get_tabular_data(tabular_dict) == (
        tabular_dict, list(tabular_dict.keys()))


def test_get_tabular_data__input_dict_no_formatters__returns_dict():
    """Check if get_tabular_data function returns the correct dictionary
    when passed a dictionary, headers and no formatters."""
    assert format_utils.get_tabular_data(
        tabular_dict,
        tabular_headers) == (tabular_dict, list(tabular_dict.keys()))


def test_get_tabular_data__input_dict__returns_dict():
    """Check if get_tabular_data function returns the correct dictionary
    when passed a dictionary, headers and formatters."""
    tabular_dict_result = {
        "A": [
            tabular_formatters["A"][0](element) for element in tabular_dict["A"]
        ],
        "B": tabular_dict["B"],
        "C": tabular_dict["C"]
    }
    assert format_utils.get_tabular_data(
        tabular_dict, tabular_headers,
        formatters=tabular_formatters) == (tabular_dict_result,
                                           list(tabular_dict.keys()))


def test_get_tabular_data__input_dict_no_headers__returns_dict():
    """Check if get_tabular_data function returns the correct dictionary
    when passed a dictionary, formatters and no headers."""
    tabular_dict_result = {
        "A": [
            tabular_formatters["A"][0](element) for element in tabular_dict["A"]
        ],
        "B": tabular_dict["B"],
        "C": tabular_dict["C"]
    }
    assert format_utils.get_tabular_data(
        tabular_dict,
        formatters=tabular_formatters) == (tabular_dict_result,
                                           list(tabular_dict.keys()))


def test_get_tabular_data__input_list_no_headers_no_formatters__returns_dict():
    """Check if get_tabular_data function returns the correct dictionary
    when passed a list and no headers or formatters."""
    assert format_utils.get_tabular_data(tabular_rows) == (tabular_rows, [])


def test_get_tabular_data__input_list_no_formatters__returns_dict():
    """Check if get_tabular_data function returns the correct dictionary
    when passed a list, headers and no formatters."""
    tabular_dict_result = {
        tabular_headers[0]: [element[0] for element in tabular_rows],
        tabular_headers[1]: [element[1] for element in tabular_rows],
        tabular_headers[2]: [element[2] for element in tabular_rows],
    }
    assert format_utils.get_tabular_data(
        tabular_rows, tabular_headers) == (tabular_dict_result, tabular_headers)


def test_get_tabular_data__input_list__returns_dict():
    """Check if get_tabular_data function returns the correct dictionary
    when passed a list, headers and formatters."""
    tabular_dict_result = {
        tabular_headers[0]: [
            tabular_formatters[tabular_headers[0]][0](element[0])
            for element in tabular_rows
        ],
        tabular_headers[1]: [element[1] for element in tabular_rows],
        tabular_headers[2]: [element[2] for element in tabular_rows],
    }
    assert format_utils.get_tabular_data(
        tabular_rows, tabular_headers,
        formatters=tabular_formatters) == (tabular_dict_result, tabular_headers)


def test_get_tabular_data__input_list_no_headers__returns_dict():
    """Check if get_tabular_data function returns the correct dictionary
    when passed a list, formatters and no headers."""
    assert format_utils.get_tabular_data(
        tabular_rows, formatters=tabular_formatters) == (tabular_rows, [])


class TestEnumClass(format_utils.CaseInsensitiveEnum):
    ITEM1 = "item1"
    ITEM2 = "item2"


def test_case_insensitive_enum_lower():
    item = TestEnumClass("item2")

    assert item.value == "item2"


def test_case_insensitive_enum_upper():
    item = TestEnumClass("ITEM1")

    assert item.value == "item1"


def test_identity_formatter():
    assert format_utils.no_formatter("test") == "test"


@mark.parametrize("size,result", [(0, "0 B"), (1, "1 B"), (1000, "1.00 KB"),
                                  (1000**2, "1.00 MB"), (1000**3, "1.00 GB"),
                                  (1000**4, "1.00 TB"), (1000**5, "1.00 PB")])
def test_bytes_formatter(size, result):
    assert format_utils.bytes_formatter(size) == result


@mark.parametrize("text, emphasis, result",
                  [("test", format_utils.Emphasis.RED, "\033[31mtest\033[0m"),
                   ("test", format_utils.Emphasis.GREEN, "\033[32mtest\033[0m"),
                   ("test", format_utils.Emphasis.BOLD, "\033[1mtest\033[0m")])
def test_emphasis_formatter(text, emphasis, result):
    assert format_utils.emphasis_formatter(text, emphasis) == result


def test_emphasis_formatter_exception():
    with pytest.raises(ValueError) as e_info:
        format_utils.emphasis_formatter("test", "Gray")
        assert "Desired emphasis is not supported." in str(e_info.value)


def test_get_ansi_formatter_ansi_enabled():
    with patch("inductiva.ansi_enabled", new=True):

        assert format_utils.get_ansi_formatter(
        ) is format_utils.emphasis_formatter


def test_get_ansi_formatter_ansi_disabled():
    with patch("inductiva.ansi_enabled", new=False):

        assert format_utils.get_ansi_formatter() is format_utils.no_formatter


@mark.parametrize("date,result", [("2021-01-01T00:00:00", "01/01, 00:00:00"),
                                  ("2021-01-01T12:00:00", "01/01, 12:00:00"),
                                  ("2021-01-01T23:59:59", "01/01, 23:59:59"),
                                  (None, None)])
def test_datetime_formatter(date, result):
    assert format_utils.datetime_formatter(date) == result


@mark.parametrize("date,result", [("2021-01-02T00:00:00", "2021-01-02 00:00"),
                                  ("2021-01-03T12:00:00", "2021-01-03 12:00"),
                                  ("2021-01-04T23:59:59", "2021-01-04 23:59"),
                                  (None, None)])
def test_datetime_formatter_ymd_hm(date, result):
    assert format_utils.datetime_formatter_ymd_hm(date) == result


@mark.parametrize(
    "delta,result",
    [(datetime.timedelta(seconds=60), "1 minute"),
     (datetime.timedelta(seconds=40), "40 seconds"),
     (datetime.timedelta(seconds=3600), "1 hour"),
     (datetime.timedelta(days=2), "2 days"),
     (datetime.timedelta(days=1, seconds=3660), "1 day, 1 hour and 1 minute")])
def test_timedelta_formatter(delta, result):
    assert format_utils.timedelta_formatter(delta) == result


@mark.parametrize("text,ident,result",
                  [("test", 2, "  test"), ("test", 0, "test"),
                   ("test", 4, "    test"),
                   ("line1\nline2", 2, "  line1\n  line2"),
                   ("line1\nline2", 0, "line1\nline2"),
                   ("line1\nline2", 4, "    line1\n    line2")])
def test_table_indenter(text, ident, result):
    # pylint: disable=protected-access
    assert format_utils._table_indenter(text, ident) == result


def test_get_tabular_str_dict():
    d = {"A": ["aa", "aaa"], "B": [1, 11], "C": ["cc", "ccc"]}
    res = format_utils.get_tabular_str(d)
    print(res)
    assert ("A" in res and "B" in res and "C" in res and "aa" in res and
            "aaa" in res and "cc" in res and "ccc" in res and "1" in res and
            "11" in res)


def test_get_tabular_str_dict_missingval():
    d = {"A": [None, "aaa"], "B": [1, None], "C": [None, "ccc"]}
    res = format_utils.get_tabular_str(d)

    assert "n/a" in res and res.count("n/a") == 3


def test_get_tabular_str_list_of_lists():
    headers = ["A", "B", "C"]
    rows = [["aa", 1, "cc"], ["aaa", 11, "ccc"]]
    res = format_utils.get_tabular_str(rows, headers)

    assert ("A" in res and "B" in res and "C" in res and "aa" in res and
            "aaa" in res and "cc" in res and "ccc" in res and "1" in res and
            "11" in res)


def test_get_tabular_str_list_of_lists_missingval():
    headers = ["A", "B", "C"]
    rows = [[None, 1, None], ["aaa", None, "ccc"]]
    res = format_utils.get_tabular_str(rows, headers)

    assert "n/a" in res and res.count("n/a") == 3


def test_get_tabular_str_list_of_lists_formatters():
    header_formatters = [lambda x: x.upper()]
    headers = ["a", "b", "c"]
    rows = [[None, 1, None], ["aaa", None, "ccc"]]
    res = format_utils.get_tabular_str(rows,
                                       headers,
                                       header_formatters=header_formatters)

    assert "A" in res and "B" in res and "C" in res


@mark.parametrize("amount, result", [
    (100000.0, "100000.00 US$"),
    (1000.0, "1000.00 US$"),
    (1.0, "1.00 US$"),
    (0.3123456, "0.31 US$"),
    (0.1234567, "0.12 US$"),
    (0.01234567, "0.012 US$"),
    (0.001234567, "0.0012 US$"),
    (0.0001234567, "0.00012 US$"),
    (0.00001234567, "0.000012 US$"),
])
def test_currency_formatter(amount, result):
    assert format_utils.currency_formatter(amount) == result
