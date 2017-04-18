#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-

from collections import namedtuple

Case = namedtuple("Case", ["case_name", "simulation", "directory_name", "label", "color"])

# Default
default = Case("default", "xfem", "", "XFEM", "blue")
default_ref = Case("default_ref", "reference", "", "REFERENCE", "red")

# Reference cases
eps_050_ref_wilkins = Case("ref_wilkins_05", "reference", "eps_050_ref_wilkins/", "REF : lump wilkins, $\epsilon= 0.5$", "red")
eps_050_ref_3x3 = Case("ref_3x3_05", "reference", "eps_050_ref_3x3/", "REF : matrice complete, $\epsilon =0.5$", "orange")
eps_025_ref_3x3 = Case("ref_3x3_025", "reference", "eps_025_ref_3x3/", "REF : matrice complete, $\epsilon = 0.25$", "orange")
eps_025_ref_wilkins = Case("ref_wilkins_025", "reference", "eps_025_ref_wilkins/", "REF : lump wilkins, $\epsilon= 0.25$", "red")

# XFEM cases
eps_050 = Case("eps_05","xfem","case_eps_050/", "XFEM : $\epsilon = 0.5$", "blue")
eps_025 = Case("eps_25","xfem","case_eps_025/", "XFEM : $\epsilon = 0.25$", "blue")



#--------------------------------------------------------------------------
# mass_1 = Case("mass_1", "xfem", "case_mass_1/", "complete mass matrix", "blue")
# mass_1 = Case("mass_1", "xfem", "case_mass_1/", "XFEM", "blue")
# mass_2 = Case("mass_2", "xfem", "case_mass_2/", "complete sans couplage", "gold")
# mass_3 = Case("mass_3", "xfem", "case_mass_3/", "lump(classic) couplage", "orange")
# mass_4 = Case("mass_4", "xfem", "case_mass_4/", "lump(classic) sans couplage", "darkorchid")
# mass_5 = Case("mass_5", "xfem", "case_mass_5/", "lump(enrich) couplage", "green")
# mass_6 = Case("mass_6", "xfem", "case_mass_6/", "lump(enrich) sans couplage", "deeppink")
# mass_7 = Case("mass_7", "xfem", "case_mass_7/", "lump(all) couplage", "steelblue")
# mass_8 = Case("mass_8", "xfem", "case_mass_8/", "lump(all) sans couplage", "olive")

# analytical_inverse = Case("analytical_inverse", "xfem", "case_analytical_inverse_complete_matrix/",
#                           "analytical inverse", "orange")


