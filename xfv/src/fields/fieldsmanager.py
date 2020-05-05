# -*- coding: utf-8 -*-
"""
Implementing field manager class

:todo: Use Singleton metaclass
"""
from collections import OrderedDict
from xfv.src.fields.enrichedfield import EnrichedField


class FieldManager(OrderedDict):
    """
    Field manager class
    """

    def __init__(self):
        super(FieldManager, self).__init__()

    def __setitem__(self, key, value):
        """
        Set a field in the manager if the field doesn't yet exist or if it is an enriched field

        :param key: name of the field
        :param value: Field object
        """
        if (key not in list(self.keys()) or
                isinstance(value, EnrichedField) and not isinstance(self[key], EnrichedField)):
            super(FieldManager, self).__setitem__(key, value)
        else:
            raise KeyError("The field {:s} already exists in FieldManager!".format(key))

    def __str__(self):
        """
        :return: information about the contents of the manager
        """
        msg = " "
        # msg = "Initial fields have been created"
        # msg = "FieldManager contents :" + os.linesep
        # msg += os.linesep.join(("{:s} <-> {:s}".format(name, field)
        #                        for name, field in self.items()))
        return msg

    def move_classical_to_enriched_fields(self, size):
        """
        Turn all classical fields into enriched ones
        """
        for name, field in list(self.items()):
            self[name] = EnrichedField(size, field.current_value, field.new_value)

    def increment_fields(self):
        """
        Increment all the fields registered in the manager
        """
        for field in list(self.values()):
            field.increment_values()
