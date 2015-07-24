#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Divers choses telles que les propriétés
"""

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ########### IMPORTATIONS DIVERSES  ####################
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
from collections import namedtuple

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ###### DEFINITION DES CLASSES & FONCTIONS  ###############
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
numerical_props = namedtuple("numerical_props",
                             ["a_pseudo", "b_pseudo", "cfl"])

geometrical_props = namedtuple("geometrical_props", ["section"])

material_props = namedtuple("material_props",
                            ["pression_init", "energie_init", "rho_init", "eos"])

properties = namedtuple("properties", ["numeric", "material", "geometric"])
