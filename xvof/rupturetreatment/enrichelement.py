#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe définissant un traitement de la rupture en enrichissant l'élément concerné
"""
from xvof.element.element1dupgraded import Element1dUpgraded
from xvof.node.node1dupgraded import Node1dUpgraded
from xvof.rupturetreatment.rupturetreatment import RuptureTreatment


class EnrichElement(RuptureTreatment):
    """
    Un traitement de rupture qui enrichit l'élément rompu
    """
    __never_enriched = True
    def __init__(self, position_rupture):
        RuptureTreatment.__init__(self)
        self.__position_rupture = position_rupture

    def applyTreatment(self, cell, *args, **kwargs):
        if EnrichElement.__never_enriched:
            if (not isinstance(cell, Element1dUpgraded)):
                print "Enrichissement de la maille : {}".format(cell)
                enrich_element = Element1dUpgraded(cell, self.__position_rupture)
                kwargs["MAILLES"][cell.indice] = enrich_element
                topologie = kwargs['TOPOLOGIE']
                enrich_nodes_indices = topologie.getNodesBelongingToCell(enrich_element)
                enrich_nodes = []
                for i in enrich_nodes_indices:
                    enrich_nodes.append(kwargs["NOEUDS"][i])
                    kwargs["NOEUDS"][i] = Node1dUpgraded(kwargs["NOEUDS"][i])
                enrich_nodes = sorted(enrich_nodes, key=lambda m : m.coordt)
                enrich_nodes[0].position_relative = -1
                enrich_nodes[1].position_relative = +1
                print "Remplacement des noeuds concernés : {}".format(enrich_nodes)
                raw_input()
                [node_l, node_r] = enrich_nodes
                kwargs["MAILLES"][cell.indice - 1].noeuds[1] = node_l
                kwargs["MAILLES"][cell.indice + 1].noeuds[0] = node_r
                kwargs["NOEUDS"][node_l.index] = node_l
                kwargs["NOEUDS"][node_r.index] = node_r
#                 node_l.elements_voisins[1] = enrich_element
#                 node_r.elements_voisins[0] = enrich_element
                EnrichElement.__never_enriched = False
        kwargs["MAILLES_ROMPUES"].remove(cell)
