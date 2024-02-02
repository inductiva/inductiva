"""Functions to handle input from the user."""

from inductiva import constants
from ..localization import translator as __


def user_confirmation_prompt(all_thins: bool, list_of_things: list,
                             all_custom_message: str,
                             big_list_custom_message: str,
                             list_costum_message: str) -> bool:
    """Prompt the user for confirmation to proceed with an action.

    Args:
        all_thins (bool): Whether the action will affect all things.
        list_of_things (list): List of things that will be affected by
            the action.
        all_custom_message (str): Custom message when we select all things.
            Example: with all_things and all_custom_message="kill all tasks"
            we will print: "You are about to kill all tasks."
        big_list_custom_message (str): Custom message when the list of
            things is too big.
            Example: with list_of_things with 100 elements and
            big_list_custom_message="kill 100 tasks" we will print:
            "You are about to kill 100 tasks."
        list_costum_message (str): Custom message when the list of things
            is small and we are about to list them all.
            Example: with list_of_things with 3 elements and
            list_costum_message="kill the following tasks" we will print:
            "You are about to kill the following tasks:"
    Returns:
        bool: Whether the user has confirmed the action.
    """
    if all_thins:
        print(f"{__("user-prompt-prefix")} {all_custom_message}.")
    else:
        if len(list_of_things) > constants.MAX_CONFIRMATION_LINES:
            print(f"{__("user-prompt-prefix")} {big_list_custom_message}.")
        else:
            print(f"{__("user-prompt-prefix")} {list_costum_message}:")
            for thing in list_of_things:
                print(f"  - {thing}")
    prompt = input("Are you sure you want to proceed (y/[N])? ")
    confirm = prompt.lower() in ["y", "ye", "yes"]
    return confirm
