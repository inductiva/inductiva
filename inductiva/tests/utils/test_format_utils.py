"""Test the format_utils module."""

from pytest import mark

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
