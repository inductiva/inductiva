"""Utils related to metaprogramming."""
import inspect


def get_type_annotations(function_ptr) -> dict:
    """Get type annotations of a given function.

    Args:
        function_ptr: Function to analyse.

    Return:
        Returns a dict with type annotations. "return" key specifies
        the type annotation of the return.
    """
    return inspect.getfullargspec(function_ptr).annotations


def get_method_name(function_ptr) -> str:
    """Constructs the name of a method supported by the Web API.

    Example: If `function_ptr` is the function `sum` placed in the
    `inductiva.math` module, it returns the string `math.sum`.

    Args:
        function_ptr: Function to analyse.

    Return:
        Returns name of the method.
    """
    module_name = function_ptr.__module__
    module_name_segments_execept_toplevel = module_name.split(".")[1:]
    module_name = ".".join(module_name_segments_execept_toplevel)

    return f"{module_name}.{function_ptr.__name__}"
