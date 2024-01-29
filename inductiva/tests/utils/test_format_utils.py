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


@mark.parametrize("input_dict, expected_str", [
    ({
        "Name": ["Ana", "Bob"],
        "Age": [25, 32],
        "Country": ["Portugal", "Spain"]
    }, "Name      Age  Country\n"
     "------  -----  ---------\n"
     "Ana        25  Portugal\n"
     "Bob        32  Spain"),
    ({
        "Name": ["Ana"],
        "Age": [25, 32],
        "Country": ["Portugal", "Spain"]
    }, "Name      Age  Country\n"
     "------  -----  ---------\n"
     "Ana        25  Portugal\n"
     "Bob        32  Spain"),
    ({
        "Name": [],
        "Age": [],
        "Country": []
    }, "Name    Age    Country"
     "------  -----  ---------"),
])
def test_get_tabular_str_test_dict(input_dict, expected_str):
    """Check if get_tabular_str function returns
    the correct string when passed a dictionary."""
    assert format_utils.get_tabular_str(input_dict) == expected_str


@mark.parametrize("input_list_rows", "input_list_headers", "expected_str",
                  [([["Ana", 25, "Portugal"], ["Bob", 32, "Spain"]
                    ], ["Name", "Age", "Country"], "Name      Age  Country\n"
                    "------  -----  ---------\n"
                    "Ana        25  Portugal\n"
                    "Bob        32  Spain"),
                   ([["Ana", 25, "Portugal"], ["Bob", 32, "Spain"]
                    ], None, "---  --  --------\n"
                    "Ana  25  Portugal\n"
                    "Bob  32  Spain\n"
                    "---  --  --------")])
def test_get_tabular_str_test_lists(input_list_rows, input_list_headers,
                                    expected_str):
    """Check if get_tabular_str function returns
    the correct string when passed a list of lists."""
    assert format_utils.get_tabular_str(input_list_rows,
                                        input_list_headers) == expected_str
