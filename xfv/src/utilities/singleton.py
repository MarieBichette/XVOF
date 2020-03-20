# -*- coding: iso-8859-1 -*-
"""
Implementing Singleton metaclass
From Python Cookbook 3rd edition Chapter 9.13 "Using metaclass to control instance creation"

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

# class Singleton(type):
#     """
#     A metaclass implementing singleton pattern
#     """
#
#     def __init__(self, *args, **kwargs):
#         self.__instance = None
#         super(Singleton, self).__init__(*args, **kwargs)
#
#     def __call__(self, *args, **kwargs):
#         if self.__instance is None:
#             self.__instance = super(Singleton, self).__call__(*args, **kwargs)
#         return self.__instance
#
#     @classmethod
#     def destroy(cls):
#         print("Suppression de la classe singleton")
#         del cls._instances[cls]
