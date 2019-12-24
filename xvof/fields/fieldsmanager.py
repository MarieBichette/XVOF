# -*- coding: iso-8859-1 -*-
"""
Implementing field manager class

:todo: Use Singleton metaclass
"""
from collections import OrderedDict
from xvof.fields.enrichedfield import EnrichedField


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
        if key not in self.keys() or isinstance(value, EnrichedField) and not isinstance(self[key], EnrichedField):
            super(FieldManager, self).__setitem__(key, value)
        else:
            raise KeyError("Le champ {:s} existe déjà dans le gestionnaire!".format(key))

    def __str__(self):
        """
        :return: informations about the contents of the manager
        """
        msg = " "
        # msg = "Initial fields have been created"
        # msg = "FieldManager contents :" + os.linesep
        # msg += os.linesep.join(("{:s} <-> {:s}".format(name, field) for name, field in self.items()))
        return msg

    def moveClassicalToEnrichedFields(self, size):
        """
        Turn all classical fields into enriched ones
        """
        for name, field in self.items():
            self[name] = EnrichedField(size, field.current_value, field.new_value)

    def incrementFields(self):
        """
        Increment all the fields registered in the manager
        """
        for field in self.values():
            field.increment_values()
