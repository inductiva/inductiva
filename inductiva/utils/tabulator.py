"""
Module for generating table representations of collections of items.
"""
from typing import Any, Optional
from dataclasses import dataclass

import tabulate


@dataclass(frozen=True)
class F:
    """
    Class for defining formatters for table columns. This class is to be used
    as a descriptor in classes that inherit from BaseTabulator, and provides
    the information needed to format a column in a table representation of a
    collection of items.

    Attributes:
        name: the name of the column.
        attr_name: the name of the attribute of the items to be used as the
            values of the column.
        attr_default: the default value to be used if the attribute is not
            present in an item.
        value_formatter: a function that takes the value of the attribute and
            the item, and returns the formatted value.
        header_formatter: a function that takes the name of the column and
            returns the formatted header.

    Example:
        In the following example, the F class is used to define how a column
        named "Person name and age" is to be formatted:

        class Person:
            name: str
            age: int

        f = F("person name and age",
              "age",
              value_formatter=lambda x,p: f'{p.name} is {x} years old',
              header_formatter=str.capitalize)

        When applied to a Person("John", 25) instance in the context of the
        BaseTabulator class, the resulting column will have the header
        "Person name and age" and each cell will contain a string of the form
        "John is 25 years old".
    """
    name: str
    attr_name: str
    attr_default: Any = None
    value_formatter: Optional[callable] = None
    header_formatter: Optional[callable] = None


class Tabulator(type):
    """Template metaclass for Table generating classes.

    Tabulators, i.e. classes that provide table-like representations of
    collections of items (lists and tuples), are classes that define a set of
    formatters, i.e. F instances, that define how each column in the
    generated table is to be produced. The purpose of this metaclass is to
    collect all F instances into a format than can be easily used by the
    downstream classes. It also enables the definition of default
    formatters for values (attributes of the items) and headers.
    """

    def __new__(metacls, cls, bases, classdict):
        header_fmtr = classdict.get('default_header_formatter', None)
        if header_fmtr is not None:
            classdict["default_header_formatter"] = staticmethod(header_fmtr)
        elif not bases:
            # assign an unset default_header_formatter to the base class
            classdict["default_header_formatter"] = None

        value_fmtr = classdict.get('default_value_formatter', None)
        if value_fmtr is not None:
            classdict["default_value_formatter"] = staticmethod(value_fmtr)
        elif not bases:
            # assign an unset default_value_formatter to the base class
            classdict["default_value_formatter"] = None

        # collect all F instances into the "formatters" attribute of the class
        classdict["formatters"] = {cls_attr:value
                                   for cls_attr, value in classdict.items()
                                   if isinstance(value, F)}

        # construct and return the class
        return super().__new__(metacls, cls, bases, classdict)



class BaseTabulator(metaclass=Tabulator):
    """Base class for table generating classes.

    This class provides a default implementation of the to_dict method, which
    generates a dictionary from a collection of items using the formatters,
    and a __call__ method that uses the to_dict method to generate a table
    representation of the items iterable using the tabulate library."""

    def to_dict(self, items):
        """Generate a dictionary from a collection of items using the formatters.

        The keys of the dictionary are the headers of the table, and the cells
        are lists of the values of the attributes of the items, formatted
        according to the formatters defined in the class.
        If a formatter does not provide a formatting function for a given
        attribute, the default value formatter is used. If the default value
        formatter is not defined, the raw value is used. Similarly, if a
        formatter does not provide a formatting function for the header, the
        default header formatter is used. If the default header formatter is
        not defined, the raw header is used.

        Args:
            items: an iterable of items to be formatted.
        """
        table = {}
        for _, formatter in self.formatters.items():
            header_fmtr = formatter.header_formatter or \
                          self.default_header_formatter
            value_fmtr = formatter.value_formatter or \
                         self.default_value_formatter
            attr_default = formatter.attr_default
            attr_name = formatter.attr_name

            if header_fmtr is None:
                header = formatter.name
            else:
                header = header_fmtr(formatter.name)

            values = (getattr(item, attr_name, attr_default) for item in items)
            if value_fmtr:
                table[header] = [value_fmtr(value, item)
                                 for value, item in zip(values, items)]
            else:
                table[header] = list(values)

        return table

    def __call__(self, items, format="simple"):
        """Generate a table representation of the items iterable using tabulate.

        The method uses the to_dict method to generate a dictionary from the
        the items iterable, with keys and values formatted according to the
        formatters defined in the class, and then uses the tabulate library to
        generate a stringified table representation of the dictionary.

        Args:
            items: an iterable of items to be formatted.
            format: the format of the table, as accepted by tabulate.
        """
        table = self.to_dict(items)
        return tabulate.tabulate(table, headers=table.keys(), tablefmt=format)



class TabulatedList(list):
    """A list subclass that applies a tabulator to its representation."""
    def __init__(self, iterable=(), formatter=None):
        """Initialize the TabulatedList instance.

        Args:
            iterable: the iterable of items to be wrapped in a TabulatedList
                instance.
            formatter: a tabulator class that provides a table representation
                of the items in the list.
        """
        super().__init__(iterable)
        self.formatter = formatter

    def __repr__(self):
        """Return a string representation of the TabulatedList instance.

        If a formatter is defined, it is used to generate a table. Otherwise,
        the default list representation is used.
        """
        if self.formatter:
            return self.formatter(self)
        return super().__repr__()

    def _repr_html_(self):
        """Return an HTML representation of the TabulatedList instance.

        The method is used by Jupyter to display the object in a notebook using
        HTML markdown.
        """
        return self.formatter(self, format="html")

def tabulated(tabulator: BaseTabulator=None):
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
        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)
            if isinstance(result, list):
                return TabulatedList(result, tabulator())
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
        default_value_formatter = lambda x,_: str(x).lower()
        column1 = F('First Name', "name",
                    value_formatter=lambda x,i: f"{x} ({i.gender})".upper())        
        column2 = F('Second Name', "none", attr_default='N/A',
                    header_formatter=str.lower)
        column3 = F("Age", "age",
                    value_formatter=lambda x,_: f'{x} years')

    @tabulated(PersonTabulator)
    def get():
        return [Person('John', 25, 'New York', 'M'),
                Person('Jane', 30, 'Los Angeles', 'F'),
                Person('Joe', 35, 'Chicago', 'M')]

    items = get()
    print(f"{isinstance(items, TabulatedList)=}")
    print(f"{isinstance(items, list)=}")
    print(f"{len(items)=}")

    print("\nusing formatter:")
    print(items)
    print("\nafter unsetting formatter:")
    items.formatter = None
    print(items)

