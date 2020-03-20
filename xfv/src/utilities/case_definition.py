#!/usr/bin/env python2.7
#  -*- coding: utf-8 -*-

import json
import sys
import os


class Case:

    def __init__(self, case_name, directory_name, label, color, linestyle, marker):
        self.case_name = case_name
        self.directory_name = directory_name
        self.label = label
        self.color = color
        self.linestyle = linestyle
        self.marker= marker


class CaseManager:

    def __init__(self):
        self.path_base = os.path.abspath(os.path.curdir)
        self.list_cases = os.path.join(self.path_base, "JSON_CASES/liste_des_cas.dat")

    def file_for_case(self, case_name):
        """
        file name for json file
        :param case_name:
        :return:
        """
        filename = case_name + '.json'
        path = os.path.join(self.path_base, "JSON_CASES")
        json_file = os.path.join(path, filename)
        return json_file

    def create_case(self, case_name, directory_name, label, color, linestyle, marker):
        """
        Create the case and register it in the case tracker file
        :param case_name:
        :param directory_name:
        :param label:
        :param color:
        :param linestyle:
        :param marker:
        :return:
        """
        # Création d'un dictionnaire à enregister dans le gestionnaire des cas json
        data = {'{:}'.format(case_name): {'directory_name': directory_name,
                                            'label': label,
                                            'color': color,
                                            'linestyle': linestyle,
                                            'marker': marker}}

        # Enregistrement du dictionnaire dans le fichier json
        fichier = self.file_for_case(case_name)
        if os.path.isfile(fichier):
            print "Overriting existing file : {:}...".format(fichier)

        with open(fichier, "w") as f:
            json.dump(data, f, indent=4, sort_keys=True,
                      separators=(',', ': '), ensure_ascii=False)

        # Enregistrement du cas dans le .dat contenant tous les noms de cas
        with open(self.list_cases, "a+") as f:
            all_cases = f.read()
            if (case_name not in all_cases):
                f.write(case_name)
                f.write("\n")
        print "Case {:} correctly registered".format(case_name)

    def find_case(self, case_name):
        """
        Build the case associated with the casename in arguments
        :param case_name: case name to build
        :return: Case
        """
        try:
            with open(self.file_for_case(case_name), "r") as f_json:
                data = json.load(f_json)
                case_dict = data[case_name]
                case = Case(case_name, case_dict["directory_name"],
                            case_dict["label"], case_dict["color"],
                            case_dict["linestyle"], case_dict["marker"])
        except:
            print "Le cas {:} n'a pas été trouvé.".format(case_name)
            choix = raw_input("Afficher tous les cas possibles ? (y/n) ")
            if choix == "y":
                print "Les cas disponibles sont :"
                with open(self.list_cases, "r") as f:
                    liste_cas = f.read()
                    print(liste_cas)
            exit(0)
        return case

    def modify_case(self, case_name, key, new_value):
        """
        Modifie la valeur de l'entrée key pour le case_name (update du fichier json)
        :param case_name:
        :param key:
        :param new_value:
        :return:
        """
        case = self.find_case(case_name)
        # import ipdb ; ipdb.set_trace()
        setattr(case, key, new_value)
        print "Le cas {:} a été modifié.".format(case_name)
        print "Sauvegarde des nouvelles valeurs"
        self.delete_case(case_name)
        self.create_case(case.case_name, case.directory_name, case.label, case.color, case.linestyle, case.marker)

    def print_info(self, case_name):
        """
        Affiche les informations sur le case case_name
        :param case_name:
        :return:
        """
        case = self.find_case(case_name)
        print u"Nom du cas : {:} ".format(case_name)
        print u"Répertoire : {:} ".format(case.directory_name)
        print u"Légende : {:} ".format(case.label)
        print u"Couleur : {:} ".format(case.color)
        print u"Linestyle : {:} ".format(case.linestyle)
        print u"Marqueur : {:} ".format(case.marker)

    def delete_case(self, case_name):
        """
        Supprime le cas
        :param case_name:
        :return:
        """
        # suppression su fichier cas.json
        os.remove(self.file_for_case(case_name))
        # suppression dans la liste des cas
        case_name_to_remove = case_name + "\n"
        with open(case_mng.list_cases, "r") as f:
            liste_cas = f.readlines()
            liste_cas.sort()
        liste_cas.remove(case_name_to_remove)
        with open(case_mng.list_cases, "w") as f_tri:
            for case in liste_cas:
                f_tri.write(case)
        print "Case {:} deleted".format(case_name)


if __name__ == '__main__':
    case_mng = CaseManager()
    msg = """Ce programme attend une option en argument :
            - avail Afficher les cas disponibles
            - create : Créer un cas
            - modify : Modifier le cas
            - info : Afficher les informations sur un cas
            - delete : Supprimer un cas
            - -h : Afficher l'aide"""

    # si le programme est lancé dans un terminal : lire les options sur quoi faire :
    if len(sys.argv) != 2:
        print msg
        raise IndexError()
    option = sys.argv[1]

    # L'option avail doit afficher les cas existants enregistrés dans le fichier liste_des_cas.dat
    if option == 'avail':
        # réécriture du fichier avec tri :
        with open(case_mng.list_cases, "r") as f:
            liste_cas = f.readlines()
            liste_cas.sort()
        with open(case_mng.list_cases, "w") as f_tri:
            for case in liste_cas:
                f_tri.write(case)
        # Affichage du fichier trié
        with open(case_mng.list_cases, "r") as f:
            liste_cas = f.read()
            print(liste_cas)

    # L'option create doit créer un cas avec les paramètres en arguments
    elif option == "create":
        print "Cette fonctionnalité va maintenant vous guider dans la création d'un nouveau cas : "
        name = raw_input("- Nom du cas : ")
        directory = raw_input("- Répertoire (à partir du répertoire xfv.src) : ")
        print "Options pour les tracés : "
        label = raw_input("- Légende : ")
        color = raw_input("- Couleur : ")
        linestyle = raw_input("- Type de trait (symbole Python) : ")
        marker = raw_input("- Marqueur (symbole Python) : ")
        case_mng.create_case(name, directory, label, color, linestyle, marker)

    elif option == "modify":
        case_name = raw_input("- Nom du cas à modifier : ")
        key = raw_input("Propriété à modifier : ")
        while key not in ["case_name", "directory_name", "label", "color", "linestyle", "marker"]:
            print ("""Propriété non reconnue. Choisir entre case_name, directory_name, label,
                                color, linestyle, marker""")
            key = raw_input("Propriété à modifier : ")
        value = raw_input("Nouvelle valeur : ")
        case_mng.modify_case(case_name, key, value)

    elif option == "delete":
        case_name = raw_input("- Nom du cas à supprimer : ")
        confirmation = raw_input("Confirmation [O/n] : ")
        if confirmation == "O":
            case_mng.delete_case(case_name)

    elif option == "info":
        case_name = raw_input("- Nom du cas à afficher : ")
        print "--- Informations sur le cas"
        case_mng.print_info(case_name)

    else:
        print msg