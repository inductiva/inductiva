"""Utilities to inspect and validate function arguments."""

import inspect
from typeguard import check_type # Beefed-up isinstance


def precondition(**conditions):
    """Decorator to validate type and value of function arguments.
    
    This decorator checks if the values of the arguments of a function match
    its respective annotations. Moreover, one can specify additional conditions
    that the arguments must satisfy. If the conditions are not met, the
    decorator raises a ValueError with the respective condition to be fulfilled.
    
    Each condition is based on a lambda function and shall be equalled to the
    name of the argument. E.g.:

        @precondition(y=lambda y: y == 10)
        def fun(x:bool, y: bool):
            return y and x
    
    Note that not all arguments need not have a condition.
    Further, if no annotation is given to an argument, the decorator will not
    check its type and/or raise any error.
    """
    def decorator(func):
        """Wrapper function to validate the function arguments."""
        signature = inspect.signature(func)
        params = signature.parameters

        def wrapper(*args, **kwargs):
            """Wrapper to validate the function arguments."""

            # Set mapping of arguments to the function's parameters.
            values_to_params_map = signature.bind(*args, **kwargs)

            # Iterate over each argument of the function
            for arg, value in values_to_params_map.arguments.items():
                # Check if each argument has an annotation
                if params[arg].annotation is not inspect.Parameter.empty:
                    # Validate the value passed conforms to annotation
                    if not check_type(value, params[arg].annotation):
                        raise TypeError(
                            f"{arg}={value} is not an instance of "
                            f"{params[arg].annotation}")
                    
                if arg in conditions and value is not None:
                    # Validate that the condition is fulfilled. 
                    if not conditions[arg](value):
                        # Convert lambda function into a string and extract
                        # the condition to be fulfilled.
                        condition = inspect.getsourcelines(
                            conditions[arg])[0][0]
                        condition = condition.strip().split(":")[1]
                        raise ValueError(
                            f"{arg}={value} does not satisfy the condition"
                            f"{condition}")
            return func(*args, **kwargs)
        return wrapper
    return decorator
