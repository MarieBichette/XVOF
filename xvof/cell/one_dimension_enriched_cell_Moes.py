# -*- coding: iso-8859-1 -*-
"""
Implementing the Element1dEnriched class for Htilde enrichment
"""
import numpy as np
from xvof.cell.one_dimension_enriched_cell import OneDimensionEnrichedCell
from xvof.cell.one_dimension_cell import OneDimensionCell
from xvof.data.data_container import DataContainer
from xvof.fields.enrichedfield import from_geometry_to_classic_field, from_geometry_to_enrich_field
from xvof.discontinuity.discontinuity import Discontinuity


# noinspection PyArgumentList,PyArgumentList,PyArgumentList,PyArgumentList,PyArgumentList,PyArgumentList
class OneDimensionMoesEnrichedCell(OneDimensionEnrichedCell):
    """
    A collection of 1d enriched elements. Treatment for Moes enrichment
    """
    def __init__(self, number_of_elements):
        super(OneDimensionMoesEnrichedCell, self).__init__(number_of_elements)

    @property
    def pressure_field(self):
        """
        :return: pressure field
        :rtype: np.array
        """
        compteur_decalage = 1
        p = np.copy(self.pressure.current_value)
        for ind in np.where(self.enriched):
            for disc in Discontinuity.discontinuity_list():
                ind_enriched = ind[0]
                # modif de pression classique -> pression gauche
                p[ind_enriched] = self.pressure.current_value[ind_enriched] - disc.additional_dof_pressure.current_value
                # on insère la pression droite dans le vecteur p
                p_droite = self.pressure.current_value[ind_enriched] + disc.additional_dof_pressure.current_value
                p = np.insert(p, ind_enriched + compteur_decalage, p_droite)
                compteur_decalage += 1
        return p

    @property
    def density_field(self):
        """
        :return: density field
        :rtype: np.array
        """
        compteur_decalage = 1
        rho = np.copy(self.density.current_value)
        for ind in np.where(self.enriched):
            for disc in Discontinuity.discontinuity_list():
                ind_enriched = ind[0]
                rho[ind_enriched] = self.density.current_value[ind_enriched] - disc.additional_dof_density.current_value
                rho_droite = self.density.current_value[ind_enriched] + disc.additional_dof_density.current_value
                rho = np.insert(rho, ind_enriched + compteur_decalage, rho_droite)
                compteur_decalage += 1
        return rho

    @property
    def energy_field(self):
        """
        :return: energy field
        :rtype: np.array
        """
        compteur_decalage = 1
        e = np.copy(self.energy.current_value)
        for ind in np.where(self.enriched):
            for disc in Discontinuity.discontinuity_list():
                ind_enriched = ind[0]
                # modif de pression classique -> pression gauche
                e[ind_enriched] = self.energy.current_value[ind_enriched] - disc.additional_dof_energy.current_value
                # on insère la pression droite dans le vecteur p
                e_droite = self.energy.current_value[ind_enriched] + disc.additional_dof_energy.current_value
                e = np.insert(e, ind_enriched + compteur_decalage, e_droite)
                compteur_decalage += 1
        return e

    @property
    def artificial_viscosity_field(self):
        """
        :return: artificial viscosity field
        :rtype: np.array
        """
        compteur_decalage = 1
        q = np.copy(self.pseudo.current_value)
        for ind in np.where(self.enriched):
            for disc in Discontinuity.discontinuity_list():
                ind_enriched = ind[0]
                # modif de pression classique -> pression gauche
                q[ind_enriched] = \
                    self.pseudo.current_value[ind_enriched] - disc.additional_dof_artificial_viscosity.current_value
                # on insère la pression droite dans le vecteur p
                q_droite = \
                    self.pseudo.current_value[ind_enriched] + disc.additional_dof_artificial_viscosity.current_value
                q = np.insert(q, ind_enriched + compteur_decalage, q_droite)
                compteur_decalage += 1
        return q

    @property
    def deviatoric_stress_field(self):
        return np.zeros([self.number_of_cells, 3])

    @property
    def stress_xx_field(self):
        return -(self.pressure_field + self.artificial_viscosity_field)

    def compute_enriched_elements_new_pressure(self):
        """
        Calcule les pressions + énergie interne + vitesse du son dans les parties gauche et droite des éléments enrichis
        puis décomposition en pression classique et pression enrichie
        Méthode générale
        :return:
        """

        if self.enriched.any():
            for disc in Discontinuity.discontinuity_list():
                mask = self.enriched  # mask ok 1 discontinuity
                # Préparation des valeurs :
                density_left = self.density.current_value[mask] - disc.additional_dof_density.current_value
                density_left_new = self.density.new_value[mask] - disc.additional_dof_density.new_value
                pressure_left = self.pressure.current_value[mask] - disc.additional_dof_pressure.current_value
                pressure_left_new = self.pressure.new_value[mask] - disc.additional_dof_pressure.new_value
                energy_left = self.energy.current_value[mask] - disc.additional_dof_energy.current_value
                energy_left_new = self.energy.new_value[mask] - disc.additional_dof_energy.new_value
                pseudo_left = self.pseudo.current_value[mask] - disc.additional_dof_artificial_viscosity.current_value
                cson_left_new = self.sound_velocity.new_value[mask] - disc.additional_dof_sound_velocity.new_value
                #
                density_right = self.density.current_value[mask] + disc.additional_dof_density.current_value
                density_right_new = self.density.new_value[mask] + disc.additional_dof_density.new_value
                pressure_right = self.pressure.current_value[mask] + disc.additional_dof_pressure.current_value
                pressure_right_new = self.pressure.new_value[mask] + disc.additional_dof_pressure.new_value
                energy_right = self.energy.current_value[mask] + disc.additional_dof_energy.current_value
                energy_right_new = self.energy.new_value[mask] + disc.additional_dof_energy.new_value
                pseudo_right = self.pseudo.current_value[mask] + disc.additional_dof_artificial_viscosity.current_value
                cson_right_new = self.sound_velocity.new_value[mask] + disc.additional_dof_sound_velocity.new_value
                # Appel de l'eos :
                energy_new_left_value, pressure_new_left_value, sound_velocity_new_left_value = \
                    OneDimensionCell.apply_equation_of_state(self, density_left, density_left_new, pressure_left,
                                                             pressure_left_new, energy_left, energy_left_new,
                                                             pseudo_left, cson_left_new)

                energy_new_right_value, pressure_new_right_value, sound_velocity_new_right_value = \
                    OneDimensionCell.apply_equation_of_state(self, density_right, density_right_new, pressure_right,
                                                             pressure_right_new, energy_right, energy_right_new,
                                                             pseudo_right, cson_right_new)
                # Assignation des résultats :
                self.pressure.new_value[mask] = from_geometry_to_classic_field(pressure_new_left_value,
                                                                               pressure_new_right_value)
                disc.additional_dof_pressure.new_value = from_geometry_to_enrich_field(pressure_new_left_value,
                                                                                  pressure_new_right_value)
                self.energy.new_value[mask] = from_geometry_to_classic_field(energy_new_left_value,
                                                                             energy_new_right_value)
                disc.additional_dof_energy.new_value = from_geometry_to_enrich_field(energy_new_left_value,
                                                                                energy_new_right_value)
                self.sound_velocity.new_value[mask] = from_geometry_to_classic_field(sound_velocity_new_left_value,
                                                                                     sound_velocity_new_right_value)
                disc.additional_dof_sound_velocity.new_value = \
                    from_geometry_to_enrich_field(sound_velocity_new_left_value, sound_velocity_new_right_value)

    def compute_enriched_elements_new_part_size(self, time_step, vecteur_vitesse_noeuds):
        """
        Calcule les nouvelles longueurs des parties gauche et droite des éléments enrichis
        puis transformation classique /enrichi
        :param time_step: time step
        :param vecteur_vitesse_noeuds: vitesse des noeuds enrichie (ddl classique)
        """
        for disc in Discontinuity.discontinuity_list():
            epsilon = disc.position_in_ruptured_element
            u1 = vecteur_vitesse_noeuds[disc.mask_in_nodes]
            u2 = vecteur_vitesse_noeuds[disc.mask_out_nodes]
            u1s = disc.additional_dof_velocity_new[0]
            u2s = disc.additional_dof_velocity_new[1]
            u1h = u1 - u1s
            u2h = u2 + u2s
            ug = (u2 - u2s) * epsilon + (u1 - u1s) * (1. - epsilon)
            ud = (u2 + u2s) * epsilon + (u1 + u1s) * (1. - epsilon)
            OneDimensionEnrichedCell.compute_new_left_right_size(time_step, disc, u1h, u2h, ug, ud)

    def compute_enriched_elements_new_density(self, topologie):
        """
        Calcule les nouvelles densités pour les éléments enrichis à partir de la conservation de la masse
        puis transformation classique / enrichi
        """
        for disc in Discontinuity.discontinuity_list():
            mask_in = topologie.cells_in_contact_with_node[disc.mask_in_nodes][0][1]
            density_left = self.density.current_value[mask_in] - disc.additional_dof_density.current_value
            density_right = self.density.current_value[mask_in] + disc.additional_dof_density.current_value
            density_left_new, density_right_new = \
                OneDimensionEnrichedCell.compute_new_left_right_density(density_left, density_right, disc)
            self.density.new_value[mask_in] = from_geometry_to_classic_field(density_left_new, density_right_new)
            disc.additional_dof_density.new_value = from_geometry_to_enrich_field(density_left_new, density_right_new)

    def compute_enriched_elements_new_pseudo(self, delta_t, topologie):
        """
        Calcule les nouvelles pseudo viscosités gauche et droite pour les éléments enrichis à partir de la methode
        compute_new_pseudo de OneDimensionCell avec les nouvelles valeurs enrichies
        :param delta_t: time_step
        :param topologie : topologie. donne la connectivité de chaque noeud / cell
        """
        for disc in Discontinuity.discontinuity_list():
            mask_in = topologie.cells_in_contact_with_node[disc.mask_in_nodes][0][1]
            # Partie gauche :
            density_left = self.density.current_value[mask_in] - disc.additional_dof_density.current_value
            density_left_new = self.density.new_value[mask_in] - disc.additional_dof_density.new_value
            sound_velocity_left = \
                self.sound_velocity.current_value[mask_in] - disc.additional_dof_sound_velocity.current_value
            pseudo_left_new = OneDimensionCell.compute_pseudo(delta_t, density_left, density_left_new,
                                                              disc.left_part_size.new_value, sound_velocity_left,
                                                              DataContainer().numeric.a_pseudo,
                                                              DataContainer().numeric.b_pseudo)
            # Partie droite :
            density_right = self.density.current_value[mask_in] + disc.additional_dof_density.current_value
            density_right_new = self.density.new_value[mask_in] + disc.additional_dof_density.new_value
            sound_velocity_right = \
                self.sound_velocity.current_value[mask_in] + disc.additional_dof_sound_velocity.current_value
            pseudo_right_new = OneDimensionCell.compute_pseudo(delta_t, density_right, density_right_new,
                                                              disc.right_part_size.new_value, sound_velocity_right,
                                                              DataContainer().numeric.a_pseudo,
                                                              DataContainer().numeric.b_pseudo)
            self.pseudo.new_value[mask_in] = from_geometry_to_classic_field(pseudo_left_new, pseudo_right_new)
            disc.additional_dof_artificial_viscosity.new_value = \
                from_geometry_to_enrich_field(pseudo_left_new, pseudo_right_new)

    def compute_enriched_elements_new_time_step(self, topologie):
        """
        Calcule les nouveaux pas de temps (qui dépendentde la taille des éléments pour les éléments enrichis à partir
        de la methode compute_new_time_step de OneDimensionCell avec les nouvelles valeurs enrichies pour les parties
        gauche et droite de l'élément enrichi
        """
        for disc in Discontinuity.discontinuity_list():
            mask_in = topologie.cells_in_contact_with_node[disc.mask_in_nodes][0][1]
            cfl = DataContainer().numeric.cfl
            cfl_pseudo = DataContainer().numeric.cfl_pseudo
            # Partie gauche
            density_left = self.density.current_value[mask_in] - disc.additional_dof_density.current_value
            density_left_new = self.density.new_value[mask_in] - disc.additional_dof_density.new_value
            sound_velocity_left_new = \
                self.sound_velocity.new_value[mask_in] - disc.additional_dof_sound_velocity.new_value
            pseudo_left = self.pseudo.current_value[mask_in] - disc.additional_dof_artificial_viscosity.current_value
            pseudo_left_new = self.pseudo.new_value[mask_in] - disc.additional_dof_artificial_viscosity.new_value

            dt_g = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, density_left, density_left_new,
                                                      disc.left_part_size.new_value, sound_velocity_left_new,
                                                      pseudo_left, pseudo_left_new)
            # Partie droite
            density_right = self.density.current_value[mask_in] + disc.additional_dof_density.current_value
            density_right_new = self.density.new_value[mask_in] + disc.additional_dof_density.new_value
            sound_velocity_right_new = \
                self.sound_velocity.new_value[mask_in] + disc.additional_dof_sound_velocity.new_value
            pseudo_right = self.pseudo.current_value[mask_in] + disc.additional_dof_artificial_viscosity.current_value
            pseudo_right_new = self.pseudo.new_value[mask_in] + disc.additional_dof_artificial_viscosity.new_value

            dt_d = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, density_right, density_right_new,
                                                      disc.right_part_size.new_value, sound_velocity_right_new,
                                                      pseudo_right, pseudo_right_new)

            self._dt[mask_in] = min(dt_g[mask_in], dt_d[mask_in])


