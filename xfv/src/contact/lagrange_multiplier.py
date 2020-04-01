# -*- coding: utf-8 -*-
"""
Definition of Contact management
"""


class LagrangianMultiplierContact(object):
    """
    A class for contact management using Lagrangian multipliers
    """
    def __init__(self, disc):
        """
        Constructeur :
        :param sigma_0 : cohesive strength
        :param separation_1 : ouverture au changement de pente
        :param stress_1 : contrainte au changement de pente
        :param delta_2 : ouverture critique (sigma = 0)
        """
        raise NotImplementedError("Méthode pas installée dans le code")
        self.disc = disc
        self._has_contact = False
        self._lambda_multiplier = 0
        self._node_mass_g = 0
        self._node_mass_d = 0

    @property
    def has_contact(self):
        return self._has_contact

    def compute_contact(self, node_velocity, time_step):
        """

        :param node_velocity:
        :param time_step:
        :return:
        """
        Mleft = self.disc.mass_matrix_enriched.get_mass_matrix_left()
        Mright = self.disc.mass_matrix_enriched.get_mass_matrix_right()
        M1 = Mleft[0, 0]
        M2enr = Mleft[1, 1]
        M2 = Mright[2, 2]
        M1enr = Mright[3, 3]
        epsilon = self.disc.discontinuity_position
        self._node_mass_g = M1 * M2enr / (M2enr * (1-epsilon) + M1 * epsilon)
        self._node_mass_d = M1enr * M2 / (M2 * (1-epsilon) + M1enr * epsilon)
        ug = (1-epsilon) * node_velocity[self.disc.mask_in_nodes] + epsilon * self.disc.additional_dof_velocity_new[0]
        ud = (1-epsilon) * self.disc.additional_dof_velocity_new[1] + epsilon * node_velocity[self.disc.mask_out_nodes]
        denom = time_step * (self._node_mass_g + self._node_mass_d) / (self._node_mass_g * self._node_mass_d)
        self._lambda_multiplier = (ug - ud).flatten() / denom

    def apply_contact(self, node_velocity, time_step):
        """

        :param node_velocity:
        :param time_step:
        :return:
        """
        node_velocity[self.disc.mask_in_nodes] -= time_step * self._lambda_multiplier / self._node_mass_g
        self.disc.additional_dof_velocity_new[1] -= time_step * (self._lambda_multiplier / self._node_mass_g).flatten()
        node_velocity[self.disc.mask_out_nodes] += time_step * self._lambda_multiplier / self._node_mass_d
        self.disc.additional_dof_velocity_new[0] += time_step * (self._lambda_multiplier / self._node_mass_d).flatten()

    def check_contact(self, node_position):
        """

        :param node_position:
        :return:
        """
        self.disc.compute_discontinuity_new_opening(node_position)
        opening_old = self.disc.discontinuity_opening.current_value
        opening_new = self.disc.discontinuity_opening.new_value
        tol = 1.e-16
        if opening_new < tol or opening_old < tol:
            self._has_contact = True
        else:
            self._has_contact = False
