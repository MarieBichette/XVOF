#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-

'Test des mask enriched'
class MyClass() :
    def __init__(self, valeur):
        self._valeur = valeur

    def multiplication(self, xval):
        self._val = self._valeur * xval
        return self._val

    @property
    def resultat(self):
        return self._val

    @property
    def modifres(self):
        res = self._val
        res += 3
        return res


classe = MyClass(4)
multipl = classe.multiplication(2)
print multipl # return 8

result = classe.resultat
print result #return 8

modif = classe.modifres
print modif #return 11
print classe.resultat #return 8
multipl2 = classe.multiplication(2)
print multipl2  #return 8