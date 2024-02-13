"""
Module for generating table-like representations of a collections of items.
"""
from typing import Any, Optional, Callable, Iterable, Mapping
from dataclasses import dataclass

import tabulate

ValueCallable = Optional[Callable[[Any, Any], Any]]
HeaderCallable = Optional[Callable[[str], Any]]


@dataclass(frozen=True)
class Col:
    """
    Class for defining formatters for table columns. This class is to be used
    as a descriptor in classes that inherit from BaseTabulator, and provides
    the information needed to format a column in a table representation of a
    collection of items.

    Attributes:
        name (str): the name of the column.
        attr_name (str): the name of the attribute of the items to be used as
            the values of the column.
        attr_default (Any): the default value to be used if the attribute is not
            present in an item.
        value_formatter: a function that takes the value of the attribute and
            the item, and returns the formatted value.
        header_formatter: a function that takes the name of the column and
            returns the formatted header.

    Example:
        In the following example, the Col class is used to define how a column
        named "Person name and age" is to be formatted:

        class Person:
            name: str
            age: int

        f = Col("person name and age",
                "age",
                value_formatter=lambda x,p: f"{p.name} is {x} years old",
                header_formatter=str.capitalize)

        When applied to a Person("John", 25) instance in the context of the
        BaseTabulator class, the resulting column will have the header
        "Person name and age" (note the capitalization) and each cell will
        contain a string of the form "John is 25 years old".
    """
    name: str
    attr_name: str
    attr_default: Any = None
    value_formatter: ValueCallable = None
    header_formatter: HeaderCallable = None

    def __call__(self,
                 item: Any,
                 default_formatter: ValueCallable = None) -> Any:
        """Apply the value formatter to the item.

        The method applies the value formatter to the attribute of the
        given item with the name attr_name, and returns the formatted value.
        If the instance"s value formatter is not defined, the given formatter
        will be used. If the latter is not defined, the raw value will be
        returned.

        Args:
            item: the item to be formatted.
            default_formatter: the default formatter to be used if the
                value_formatter is not defined.
        """
        value = getattr(item, self.attr_name, self.attr_default)
        # fallback to the identity method when no valid formatter is available
        fmtr = self.value_formatter or default_formatter or self._id
        return fmtr(value, item)

    def get_formatted_header(self,
                             default_formatter: HeaderCallable = None) -> Any:
        """Get the formatted header of the column.

        The method applies the header formatter to the name of the column, and
        returns the formatted header. If the instance"s header formatter is not
        defined, the given formatter will be used. If the latter is not defined,
        the raw header will be returned.

        Args:
            default_formatter: the default formatter to be used if the
                header_formatter is not defined.
        """
        fmtr = self.header_formatter or default_formatter or self._id
        return fmtr(self.name)

    def _id(self, item: Any) -> Any:
        return item


class Tabulator(type):
    """Template metaclass for Table generating classes.

    Tabulators, i.e. classes that provide table-like representations of
    collections of items (lists and tuples), are classes that define a set of
    formatters, i.e. Col instances, that describe how each column in the
    generated table is to be generated. The purpose of this metaclass is to
    collect all Col instances into a format than can be easily used by the
    downstream classes. It also enables the definition of default
    formatters for values (attributes of the items) and headers.
    """

    def __new__(mcs, cls, bases, classdict):
        header_fmtr = classdict.get("default_header_formatter", None)
        if header_fmtr is not None:
            classdict["default_header_formatter"] = staticmethod(header_fmtr)
        elif not bases:
            # assign an unset default_header_formatter to the base class
            classdict["default_header_formatter"] = None

        value_fmtr = classdict.get("default_value_formatter", None)
        if value_fmtr is not None:
            classdict["default_value_formatter"] = staticmethod(value_fmtr)
        elif not bases:
            # assign an unset default_value_formatter to the base class
            classdict["default_value_formatter"] = None

        # collect all F instances into the "formatters" attribute of the class
        classdict["columns"] = {
            cls_attr: value
            for cls_attr, value in classdict.items()
            if isinstance(value, Col)
        }

        # construct and return the class
        return super().__new__(mcs, cls, bases, classdict)


class BaseTabulator(metaclass=Tabulator):
    """Base class for table generating classes.

    This class provides a default implementation of the to_dict method, which
    generates a dictionary from a collection of items using Col(umns),
    and a __call__ method that uses the to_dict method to generate a table
    representation of the items iterable using the tabulate library."""

    def __init__(self, tablefmt: str = "simple"):
        self.tablefmt = tablefmt

    def to_dict(self,
                items: Iterable[Iterable[Any]]) -> Mapping[str, Iterable[Any]]:
        """Generate a dictionary from a collection of items using formatters.

        The keys of the dictionary are the headers of the table, and the cells
        are lists of the values of the attributes of the items, formatted
        according to the formatters defined in the class.

        Args:
            items: an iterable of items to be formatted.
        """
        def_header_fmtr = self.default_header_formatter
        def_value_fmtr = self.default_value_formatter
        if not self.columns:
            print("Warning: no columns to display")
        return {
            fmtr.get_formatted_header(def_header_fmtr): [
                fmtr(item, def_value_fmtr) for item in items
            ] for fmtr in self.columns.values()
        }

    def __call__(self,
                 items: Iterable[Iterable[Any]],
                 tablefmt: Optional[str] = None) -> str:
        """Generate a table representation of the items iterable using tabulate.

        The method uses the to_dict method to generate a dictionary from the
        the items iterable, with keys and values formatted according to the
        columns defined in the class, and then uses the tabulate library to
        generate a stringified table representation of the dictionary.

        Args:
            items: an iterable of items to be formatted.
            tablefmt: the format of the table, as accepted by tabulate.
                Used to override the value of the `tablefmt` attribute
                used when instantiating the class
        """
        table = self.to_dict(items)
        return tabulate.tabulate(table,
                                 headers="keys",
                                 tablefmt=tablefmt or self.tablefmt)


class TabulatedList(list):
    """A list subclass that applies a tabulator to its representation."""

    def __init__(self, iterable=(), tabulator=None):
        """Initialize the TabulatedList instance.

        Args:
            iterable: the iterable of items to be wrapped in a TabulatedList
                instance.
            tabulator: a tabulator class that provides a table representation
                of the items in the list.
        """
        super().__init__(iterable)
        self.tabulator = tabulator

    def __repr__(self) -> str:
        """Return a string representation of the TabulatedList instance.

        If a tabulator is defined, it is used to generate a table. Otherwise,
        the default list representation is used.
        """
        if self.tabulator:
            return self.tabulator(self)
        return super().__repr__()

    def _repr_html_(self) -> str:
        """Return an HTML representation of the TabulatedList instance.

        The method is used by Jupyter to display the object in a notebook using
        HTML markdown.
        """
        return self.tabulator(self, format="unsafehtml")


def tabulated(tabulator: BaseTabulator = None):
    """Decorator to apply a tabulator to the return value of a function.

    The decorator takes a tabulator class as an argument, and returns a
    decorator that applies the tabulator to the return value of the decorated
    function. If the return value is a list, it is wrapped in a TabulatedList
    instance that uses the given tabulator class to generate a table
    representation of the list. If the return value is not a list, the tabulator
    is not applied, and the return value is passed through unchanged.

    Args:
        tabulator: a class that inherits from BaseTabulator, and provides a
            table representation of a collection of items.
    """

    def decorator(func):

        if tabulator is None:
            return func

        if not isinstance(tabulator, BaseTabulator):
            raise TypeError("tabulator must be a subclass of BaseTabulator")

        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)
            if isinstance(result, list):
                return TabulatedList(result, tabulator)
            return result

        return wrapper

    return decorator


if __name__ == "__main__":

    @dataclass
    class Person:
        name: str
        age: int
        birthplace: str
        gender: str

    class PersonTabulator(BaseTabulator):
        default_header_formatter = lambda x: f"<{x.upper()}>"
        default_value_formatter = lambda x, _: str(x).lower()
        column1 = Col("First Name",
                      "name",
                      value_formatter=lambda x, i: f"{x} ({i.gender})".upper())
        column2 = Col("Second Name",
                      "none",
                      attr_default="N/A",
                      header_formatter=str.lower)
        column3 = Col("Age", "age", value_formatter=lambda x, _: f"{x} years")

    @tabulated(PersonTabulator("grid"))
    def get():
        return [
            Person("John", 25, "New York", "M"),
            Person("Jane", 30, "Los Angeles", "F"),
            Person("Joe", 35, "Chicago", "M")
        ]

    people = get()
    print(f"{isinstance(people, TabulatedList)=}")
    print(f"{isinstance(people, list)=}")
    print(f"{len(people)=}")

    print("\nusing formatter:")
    print(people)
    print("\nafter unsetting formatter:")
    people.formatter = None
    print(people)
