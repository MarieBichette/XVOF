"""
This module implements the TypeCheckedDataClass
"""
from dataclasses import dataclass, fields
from functools import partial
from operator import le, lt
from typing import Union, List


@dataclass
class TypeCheckedDataClass:
    """
    This class is the base class of dataclasses that wants to check the type of their fields

    Inspired from : https://stackoverflow.com/questions/54863458/force-type-conversion-in-python-dataclass-init-method  # pylint:disable=line-too-long
    """
    def _raise_type_error(self, field_name, field_type, value):
        """
        A uniform way of throwing a type error
        """
        raise ValueError(f'During filling of class {self.__class__.__name__}:\n'
                         f'{field_name} should be of type {field_type}, '
                         f'but was initialized with {repr(value)}')

    def __post_init__(self):
        # TODO: implements the type checking recursively to hold for
        # Optional[Union[str, float]] for example 
        for _field in [f_ for f_ in fields(self) if f_.init]:
            value = getattr(self, _field.name)
            try:
                # try to take into account typing generics
                # import ipdb; ipdb.set_trace()
                if _field.type.__origin__ == Union:  # case of Union or Optional types
                    if not any(isinstance(value, _type) for _type in _field.type.__args__):
                        self._raise_type_error(_field.name, _field.type, value)
                elif _field.type.__origin__ in (list, List):  # case of List types
                    if not all(isinstance(_val, _field.type.__args__[0]) for _val in value):
                        self._raise_type_error(_field.name, _field.type, value)
                else:
                    raise RuntimeError(f"Unable to determine type of {value}")
            except AttributeError:
                # the _field type is not a typing generic and thus has not __origin__ attribute
                if not isinstance(value, _field.type):
                    self._raise_type_error(_field.name, _field.type, value)

    def _ensure_predicat(self, predicat, error_string, *args):
        for _field in [f_ for f_ in fields(self) if f_.init]:
            val = getattr(self, _field.name)
            if _field.name in args and val and not predicat(val):
                raise ValueError(error_string.format(val, _field.name))

    def _ensure_positivity(self, *args):
        predicat = partial(le, 0)
        error_string = ("{} < 0!\n"
                        "Please provide a value greater or equal to 0 for the coefficient {:s}")
        self._ensure_predicat(predicat, error_string, *args)

    def _ensure_strict_positivity(self, *args):
        predicat = partial(lt, 0)
        error_string = ("{} <= 0!\n"
                        "Please provide a value greater than 0 for the coefficient {:s}")
        self._ensure_predicat(predicat, error_string, *args)

    def _ensure_value_in(self, field_name, authorized_values):
        val = getattr(self, field_name)
        if val not in authorized_values:
            raise ValueError(f"{val} not in {authorized_values}!"
                             f"Please provide a value in {authorized_values} for "
                             f"the coefficient {field_name}")

    def _ensure_list_value_in(self, field_name, authorized_values):
        field_name_list = getattr(self, field_name)
        for elem in field_name_list:
            if elem not in authorized_values:
                raise ValueError(f"{elem} not in {authorized_values}!"
                                 f"Please provide a value in {authorized_values} for "
                                 f"the coefficient {field_name}")

    def _ensure_defined(self, field_name, class_name, json_path):
        val = getattr(self, field_name)
        if val is None:
            raise ValueError(f"Missing value {json_path} in data file. "
                             f"Cannot build the {field_name} field for class {class_name}")
