"""
This module implements the TypeCheckedDataClass
"""
from dataclasses import dataclass, fields


@dataclass
class TypeCheckedDataClass:
    """
    This class is the base class of dataclasses that wants to check the type of their fields

    Inspired from : https://stackoverflow.com/questions/54863458/force-type-conversion-in-python-dataclass-init-method  # pylint:disable=line-too-long
    """
    def __post_init__(self):
        # todo implements the type checking recursively to hold for 
        # Optional[Union[str, float]] for example 
        for field in fields(self):
            value = getattr(self, field.name)
            try:
                # try to take into account typing generics
                if field.type._name == 'Union':
                    if not any (isinstance(value, _type) for _type in field.type.__args__):
                        raise ValueError(f'During filling of class {self.__class__.__name__}:\n'
                                        f'{field.name} should be of type {field.type}, '
                                        f'but was initialized with {repr(value)}')
                elif field.type._name == 'List':
                    if not all (isinstance(_val, field.type.__args__[0]) for _val in value):
                        raise ValueError(f'During filling of class {self.__class__.__name__}:\n'
                                        f'{field.name} should be of type {field.type}, '
                                        f'but was initialized with {repr(value)}')
            except AttributeError:
                # the field type is not a typing generic and thus has not _name attribute
                if not isinstance(value, field.type):
                    raise ValueError(f'During filling of class {self.__class__.__name__}:\n'
                                     f'{field.name} should be of type {field.type}, '
                                     f'but was initialized with {repr(value)}')
                
