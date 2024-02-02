"""Functions to handle input from the user."""

from inductiva import constants


def user_confirmation_prompt(all_thins: bool, list_of_things: list,
                             all_custom_message: str,
                             big_list_custom_message: str,
                             list_costum_message: str) -> bool:
    """Prompt the user for confirmation to proceed with an action.

    Args:
        all_thins (bool): Whether the action will affect all things.
        list_of_things (list): List of things that will be affected by
            the action.
    Returns:
        bool: Whether the user has confirmed the action.
    """
    if all_thins:
        print(f"You are about to {all_custom_message}.")
    else:
        if len(list_of_things) > constants.MAX_CONFIRMATION_LINES:
            print(f"You are about to {big_list_custom_message}.")
        else:
            print(f"You are about to {list_costum_message}:")
            for thing in list_of_things:
                print(f"  - {thing}")
    prompt = input("Are you sure you want to proceed (y/[N])? ")
    confirm = prompt.lower() in ["y", "ye", "yes"]
    return confirm
