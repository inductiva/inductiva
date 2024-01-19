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
