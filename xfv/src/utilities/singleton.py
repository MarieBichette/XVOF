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
    """
    A metaclass implementing singleton pattern
    """

    def __init__(self, *args, **kwargs):
        """
        Construction of singleton
        :param args:
        :param kwargs:
        """
        self.__instance = None
        super(Singleton, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        """
        Call of the singleton
        :param args:
        :param kwargs:
        :return:
        """
        if self.__instance is None:
            self.__instance = super(Singleton, self).__call__(*args, **kwargs)
        return self.__instance

    def clear(self):
        """
        Delete the singleton
        :return:
        """
        del self.__instance
        self.__instance = None


if __name__ == "__main__":
    import doctest
    doctest.testmod()
