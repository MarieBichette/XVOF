# -*- coding: iso-8859-1 -*-
"""
Implementing Singleton metaclass
From https://stackoverflow.com/questions/50065276/clearing-a-metaclass-singleton"

>>> class Spam(object):
...     __metaclass__ = Singleton
...     def __init__(self):
...         print "Creating spam"
>>> a = Spam()
Creating spam
>>> b = Spam()
>>> a is b
True
>>> c = Spam()
>>> a is c
True
"""


class Singleton(type):
    _instances = {}

    # Each of the following functions use cls instead of self
    # to emphasize that although they are instance methods of
    # Singleton, they are also *class* methods of a class defined
    # with Singleton
    def __call__(cls, *args, **kwargs):
        if cls not in Singleton._instances:
            Singleton._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return Singleton._instances[cls]

    def clear(cls):
        try:
            del Singleton._instances[cls]
        except KeyError:
            pass