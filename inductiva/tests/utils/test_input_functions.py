"""Test the input_functions method."""

from inductiva.utils.input_functions import user_confirmation_prompt


def test_input_functions_confirm(monkeypatch):

    monkeypatch.setattr("builtins.input", lambda _: "yes")

    l = [1, 2, 3]

    result = user_confirmation_prompt(l, "all", "unlisted", "listed", True)

    assert result is True


def test_input_functions_not_confirm(monkeypatch):

    monkeypatch.setattr("builtins.input", lambda _: "no")

    l = [1, 2, 3]

    result = user_confirmation_prompt(l, "all", "unlisted", "listed", True)

    assert result is False


def test_input_functions_big_list(monkeypatch):

    monkeypatch.setattr("builtins.input", lambda _: "no")

    l = [1, 2, 3, 4, 5, 6]

    result = user_confirmation_prompt(l, "all", "unlisted", "listed", True)

    assert result is False
