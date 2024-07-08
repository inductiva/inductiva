from base_tabulator import BaseTabulator, Col, tabulated
from dataclasses import dataclass


@dataclass
class Person:
    name: str
    age: int
    birthplace: str
    gender: str


class PersonTabulator(BaseTabulator):
    default_header_formatter = lambda x: str.upper(x)
    default_value_formatter = lambda x,_: str.lower(x)

    fname = Col('First Name', "name", description="the first name of the person")
    age = Col('Age', "age", value_formatter=lambda x,p: f"{x} years ({p.gender})")
    birth = Col('Birth place', "birthplace")
    gender = Col('Gender', "gender")
    treta = Col('Treta', "treta", "n/a")


@tabulated(PersonTabulator, tablefmt="plain")
def get():
    return [Person('John', 25, 'New York', 'M'),
            Person('Jane', 30, 'Los Angeles', 'F'),
            Person('Joe', 35, 'Chicago', 'M')]

@tabulated(PersonTabulator, tablefmt="simple")
def listx():
    return [Person('John2', 25, 'New York', 'M'),
            Person('Jane2', 30, 'Los Angeles', 'F'),
            Person('Joe2', 35, 'Chicago', 'M')]


items = get()
# print(items.tabulator.columns)
print(items.tabulator.columns)
print({name:col.description for name, col in items.tabulator.columns.items()})
print(items)
print("****")

items.tabulator.fname.disable()
getattr(items.tabulator, "treta").disable()

print(items)

def on_pre_tabulate(x):
    x.age += 100
items.tabulator.on_pre_tabulate = on_pre_tabulate

print(items)

# items2 = items[1:]
# print(items2)
# print(type(items))
# print(type(items[0]))


# print(items.tabulator.columns.keys())
# print(PersonTabulator.columns.keys())


# items2 = listx()
# print(items2)
