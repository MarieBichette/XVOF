#!/usr/bin/env python2.7
#  -*- coding: iso-8859-1 -*-

from collections import namedtuple
import xlrd


Case = namedtuple("Case", ["case_name", "directory_name", "label", "color", "linestyle", "marker"])

equivalence_python = dict()
equivalence_python["Tiret"] = "--"
equivalence_python["Ligne"] = "-"
equivalence_python["Point"] = "."
equivalence_python["Diamond"] = "d"

def build_cases_from_case_names(case_name_list, directory_for_case_definition = "//home/marie/Documents/cases.xls"):
    """
    Retourne une liste de cas à post traiter à partir d'une liste de nom de cas
    :param case_name_list: liste des noms de cas à traiter
    :return: liste des cas associés
    """
    case_list = []

    # Ouverture du fichier XLS
    wb = xlrd.open_workbook(directory_for_case_definition)

    sh = wb.sheet_by_name(u'CaseList')

    for case_name in case_name_list:
        print "Searching for case {:}".format(case_name)
        # Pour chaque nom de cas : recherche de la ligne
        rownum = 0
        ligne = sh.row_values(rownum, start_colx=0, end_colx=None)

        # Recherche de la ligne du fichier excel qui contient le nom demandé
        while (rownum < sh.nrows) & (ligne[0] != case_name):
            ligne = sh.row_values(rownum, start_colx=0, end_colx=None)
            print ligne
            rownum += 1


        # Vérification qu'on a trouvé la bonne ligne dans le tableau xls
        if rownum == sh.nrows:
            raise ValueError("""
            Impossible de retrouver le case {:} dans le fichier de données {:}""".format(case_name,
                                                                                         directory_for_case_definition))
        # Création du cas à partir des données du fichier excel
        new_case = Case(ligne[0], ligne[1], ligne[2], ligne[3],
                        equivalence_python[ligne[4]], equivalence_python[ligne[5]])

        # Ajout du cas à la liste de retour
        case_list.append(new_case)
    return case_list

if __name__ == '__main__':
    case_name_list = ["czm", "default"]
    build_cases_from_case_names(case_name_list)