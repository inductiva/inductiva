"""Test the format_utils module."""

from pytest import mark

from inductiva.utils import format_utils


@mark.parametrize("time_input, expected_time", [(60, "01m 00s"),
                                                (40, "00m 40s"),
                                                (3600, "1h 00m 00s"),
                                                (3660, "1h 01m 00s"),
                                                (3615, "1h 00m 15s"),
                                                (645, "10m 45s")])
def test_seconds_formatter__valid_time__formatted_time(time_input,
                                                       expected_time):
    """Check that the seconds formatter returns the correct string."""
    assert format_utils.seconds_formatter(time_input) == expected_time
