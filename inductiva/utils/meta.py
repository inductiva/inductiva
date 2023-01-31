"""
Utils related to metaprogramming.
"""
import inspect


def get_type_annotations(function_ptr) -> dict:
    """
    Get type annotations of a given function.
    """
    return inspect.get_annotations(function_ptr)


def get_method_name(function_ptr) -> str:
    """
    Constructs the name of a method name
    """
    module_name = function_ptr.__module__
    module_name_segments_execept_toplevel = module_name.split(".")[1:]
    module_name = ".".join(module_name_segments_execept_toplevel)

    return f"{module_name}.{function_ptr.__name__}"