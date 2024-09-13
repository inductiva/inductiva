"""Functions to handle input from the user."""

from inductiva import constants
from ..localization import translator as __


def user_confirmation_prompt(items: list, all_msg: str, unlisted_msg: str,
                             listed_msg: str, is_all: bool) -> bool:
    """Prompt the user for confirmation to proceed with an action.

    Args:
        is_all (bool): Whether the action will affect all things.
        items (list): List of things that will be affected by
            the action.
        all_msg (str): Custom message when we select all things.
            Example: with all_things and all_msg="kill all tasks"
            we will print: "You are about to kill all tasks."
        unlisted_msg (str): Custom message when the list of
            things is too big.
            Example: with items with 100 elements and
            unlisted_msg="kill 100 tasks" we will print:
            "You are about to kill 100 tasks."
        listed_msg (str): Custom message when the list of things
            is small and we are about to list them all.
            Example: with items with 3 elements and
            listed_msg="kill the following tasks" we will print:
            "You are about to kill the following tasks:"
    Returns:
        bool: Whether the user has confirmed the action.
    """
    if is_all:
        print(all_msg)
    else:
        if len(items) > constants.MAX_CONFIRMATION_LINES:
            print(unlisted_msg)
        else:
            print(listed_msg)
            for thing in items:
                print(f"  - {thing}")
    prompt = input(__("user-prompt-confirmation"))
    confirm = prompt.lower() in ["y", "ye", "yes"]
    return confirm
