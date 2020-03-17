#!/usr/bin/env python2.7
"""
A collection of methods for building True fields from classical and enriched ones
"""
import numpy as np


def build_node_true_field(classical_field, enriched_field, node_status, enrichment_type=None):
    """
    Build the node true field based on the node status

    :param classical_field: field of classical values
    :param enriched_field: field of enriched values
    :param node_status: boolean mask where True indicates an enriched item
    :param enrichment_type: type of enrichment
    :return: the node true field
    """
    true_field = np.copy(classical_field)
    if not enrichment_type == "Hansbo":
        # Hansbo : vitesse nodale = vitesse classique => rien a faire. On verifie quand meme
        raise ValueError(
            """Don't know how to build true fields with enrichment type {}. Possibilities are Hansbo and Moes""".format(
                enrichment_type))
    return true_field.reshape((len(true_field), 1))


def build_cell_true_field(classical_field, enriched_field, enrichment_type):
    """
    Build the cell true field based on the cell status

    :param classical_field: field of classical values
    :param enriched_field: field of enriched values
    :param enrichment_type: type of enrichment. Moes ou Hansbo (utile pour reconstruction)
    :return the cell true field

    >>> import numpy as np
    >>> a = np.array([1., 2., 1.])
    >>> b = np.array([0., 0.5, 0.])
    >>> s = np.array([False, True, False])
    >>> build_cell_true_field(a, b, s).tolist()
    [1.0, 1.5, 2.5, 1.0]
    >>> b = np.array([0.25, 0., 0.])
    >>> s = np.array([True, False, False])
    >>> build_cell_true_field(a, b, s).tolist()
    [0.75, 1.25, 2.0, 1.0]
    >>> a = np.array([1., -2., 2., 3., -7, 10.])
    >>> b = np.array([0., -2., 0., 0., 2., 0.])
    >>> s = np.array([False, True, False, False, True, False])
    >>> build_cell_true_field(a, b, s).tolist()
    [1.0, 0.0, -4.0, 2.0, 3.0, -9.0, -5.0, 10.0]
    >>> a = np.array([-3., 2., 1., -3., -5, 9.])
    >>> b = np.array([0., -2., 3., 0., 0., 0.])
    >>> s = np.array([False, True, True, False, False, False])
    >>> build_cell_true_field(a, b, s).tolist()
    [-3.0, 4.0, 0.0, -2.0, 4.0, -3.0, -5.0, 9.0]
    """
    # ------------- extract sigma xx if field = tensor
    # if len(classical_field.value.flatten()) > classical_field.shape[0]:  # si tableau de dimension > 1
    #     # correspond a un field de type tenseur (3 valeurs enregistrees pour chaque cell)
    #     classical_field = classical_field[:, 0]  # on ne garde que la composante xx

    # idem pour le champ enrichi (operation transparente pour les champs scalaires)
    # enriched_field = enriched_field[:, 0:2]  # on ne garde que cell_id et la composante xx
    # -------------------------------------------------------------------

    if len(classical_field.flatten()) == classical_field.shape[0]:
        # si le champ est un champ scalaire, on rentre dans la boucle et on reshape le np.array pour avoir
        # deux dimensions (dont une egale a 1). Cela permet d'avoirla compatibilite avec le type des coordonnees
        # des items pour np.concatenate
        classical_field = classical_field.reshape((len(classical_field), 1))

    true_field = np.copy(classical_field)
    offset = 0

    # on recupere les cell_id dans la premiere colonne de la db
    # la methode de reconstruction des champs necessite que les indices des cells enrichies soient triees
    enriched_field = np.sort(enriched_field, 0)
    enriched_id = enriched_field[:, 0]

    if enrichment_type == 'Hansbo':
        for stable_index in enriched_id:
            # import ipdb ; ipdb.set_trace()
            moving_index = int(stable_index + offset)
            # true_field = np.insert(true_field, moving_index + 1, enriched_field[offset, 1])
            true_field = np.insert(true_field, moving_index + 1, enriched_field[offset, 1:])
            # reshape car le np.insert a casse la structure tenseur [n, 3].
            # Le np.reshape sert a retrouve une shape de [n+1, 3] apres insertion
            true_field = true_field.reshape(-1, classical_field.shape[1])
            offset += 1
    else:
        raise ValueError(
            """Don't know how to build true fields with enrichment type {}. Possibilities are Hansbo and Moes""".format(
                enrichment_type))

    return true_field


if __name__ == "__main__":
    import doctest
    doctest.testmod()
