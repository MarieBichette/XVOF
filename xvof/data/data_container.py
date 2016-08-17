"""
Implementing the DataContainer class
"""

from collections import namedtuple

import lxml.etree as et

from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.rupturetreatment.enrichelement import EnrichElement
from xvof.rupturetreatment.imposedpressure import ImposedPressure
from xvof.utilities.singleton import Singleton

numerical_props = namedtuple("numerical_props", ["a_pseudo", "b_pseudo", "cfl", "cells_number"])

geometrical_props = namedtuple("geometrical_props", ["section", "length"])

material_props = namedtuple("material_props", ["pression_init", "temp_init", "rho_init", "energie_init", "eos",
                                               "damage_treatment", "damage_treatment_value"])

time_props = namedtuple("time_props", ['initial_time_step', 'final_time'])

output_props = namedtuple("output_props", ['number_of_images', 'images_dump', 'images_show'])


class DataContainer(object):
    """
    Contains the data read from the datafile
    """
    __metaclass__ = Singleton

    def __init__(self, datafile_name="XDATA.xml"):
        '''
        Constructor
        '''
        self.__datadoc = et.parse(datafile_name)
        self.numeric = numerical_props(*self.__fillInNumericalProperties())
        self.geometric = geometrical_props(*self.__fillInGeometricalProperties())
        self.material = material_props(*self.__fillInMaterialProperties())
        self.time = time_props(*self.__fillInTimeProperties())
        self.output = output_props(*self.__fillInOutputProperties())

    def __fillInNumericalProperties(self):
        """
        :return: the numerical properties :
            - linear and quadratic artifical viscosity coefficients;
            - CFL;
            - number of cells
        :rtype: tuple(float, float, float, int)
        """
        b_pseudo = float(self.__datadoc.find('numeric-parameters/linear-pseudo').text)
        a_pseudo = float(self.__datadoc.find('numeric-parameters/quadratic-pseudo').text)
        cfl = float(self.__datadoc.find('numeric-parameters/cfl').text)
        cells_number = float(self.__datadoc.find('numeric-parameters/number-of-elements').text)
        return a_pseudo, b_pseudo, cfl, cells_number

    def __fillInGeometricalProperties(self):
        """
        :return: the geometric properties:
            - area of the cell;
            - length of the rod
        :rtype: tuple(float, float)
        """
        section = float(self.__datadoc.find('geometry/section').text)
        length = float(self.__datadoc.find('geometry/length').text)
        return section, length

    def __fillInMaterialProperties(self):
        """
        :return: the material properties:
            - initialization pressure, temperature, density and internal energy
            - eos
        :rtype: tuple(float, float, float, float, EquationOfStateBase)
        """
        init_pressure = float(self.__datadoc.find('matter/initialization/initial-pressure').text)
        init_temperature = float(self.__datadoc.find('matter/initialization/initial-temperature').text)
        init_density = float(self.__datadoc.find('matter/initialization/initial-density').text)
        init_internal_energy = float(self.__datadoc.find('matter/initialization/initial-internal-energy').text)
        if self.__datadoc.find('matter/equation-of-state').text == 'Mie-Gruneisen':
            eos = MieGruneisen()
        else:
            raise ValueError("Only MieGruneisen's equation of state is available")
        dmg_treatment_name = str(self.__datadoc.find('matter/damage-treatment/name').text)
        if dmg_treatment_name == "ImposedPressure":
            dmg_treatment = ImposedPressure
        elif dmg_treatment_name == "Enrichment":
            dmg_treatment = EnrichElement
        else:
            raise ValueError("Only 'ImposedPressure' or 'Enrichment' are possible values")
        dmg_treatment_value = float(self.__datadoc.find('matter/damage-treatment/value').text)
        return init_pressure, init_temperature, init_density, init_internal_energy, eos, \
               dmg_treatment, dmg_treatment_value

    def __fillInTimeProperties(self):
        """
        :return: time properties :
            - initial time step
            - final time
        :return: tuple(float, float)
        """
        initial_time_step = float(self.__datadoc.find('time-management/initial-time-step').text)
        final_time = float(self.__datadoc.find('time-management/final-time').text)
        return initial_time_step, final_time

    def __fillInOutputProperties(self):
        """
        :return: output properties:
            - number of images
            - is dump of images required?
            - is display of images required?
        :tuple(int, bool, bool)
        """
        number_of_images = int(self.__datadoc.find('output/number-of-images').text)
        images_dump = self.__datadoc.find('output/dump-images').text.lower() == 'true'
        images_show = self.__datadoc.find('output/show-images').text.lower() == 'true'
        return number_of_images, images_dump, images_show

    def hasExternalSolver(self):
        if self.__datadoc.find('numeric-parameters/external-solver-library') is not None:
            return True
    
    def getExternalSolverPath(self):
        return self.__datadoc.find('numeric-parameters/external-solver-library').text
