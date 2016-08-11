"""
Implementing the DataContainer class
"""

from collections import namedtuple

import lxml.etree as et

from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.utilities.singleton import Singleton

numerical_props = namedtuple("numerical_props", ["a_pseudo", "b_pseudo", "cfl"])

geometrical_props = namedtuple("geometrical_props", ["section", "length"])

material_props = namedtuple("material_props", ["pression_init", "temp_init", "rho_init", "energie_init", "eos"])


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

    def __fillInNumericalProperties(self):
        """
        :return: the numerical properties (linear and quadratic artifical viscosity parameters and CFL)
        :rtype: tuple
        """
        b_pseudo = float(self.__datadoc.find('numeric-parameters/linear-pseudo').text)
        a_pseudo = float(self.__datadoc.find('numeric-parameters/quadratic-pseudo').text)
        cfl = float(self.__datadoc.find('numeric-parameters/cfl').text)
        return a_pseudo, b_pseudo, cfl

    def __fillInGeometricalProperties(self):
        """
        :return: the geometric properties (area of the element, length of the rod)
        :rtype: tuple
        """
        section = float(self.__datadoc.find('geometry/section').text)
        length = float(self.__datadoc.find('geometry/length').text)
        return section, length

    def __fillInMaterialProperties(self):
        init_pressure = float(self.__datadoc.find('matter/initialization/initial-pressure').text)
        init_temperature = float(self.__datadoc.find('matter/initialization/initial-temperature').text)
        init_density = float(self.__datadoc.find('matter/initialization/initial-density').text)
        init_internal_energy = float(self.__datadoc.find('matter/initialization/initial-internal-energy').text)
        if self.__datadoc.find('matter/equation-of-state').text == 'Mie-Gruneisen':
            eos = MieGruneisen()
        else:
            raise ValueError("Only MieGruneisen's equation of state is available")
        return init_pressure, init_temperature, init_density, init_internal_energy, eos

    def getFinalTime(self):
        return float(self.__datadoc.find('time-management/final-time').text)
    
    def getInitialTimeStep(self):
        return float(self.__datadoc.find('time-management/initial-time-step').text)

    def getNumberOfElements(self):
        return float(self.__datadoc.find('numeric-parameters/number-of-elements').text)

    def hasExternalSolver(self):
        if self.__datadoc.find('numeric-parameters/external-solver-library') is not None:
            return True
    
    def getExternalSolverPath(self):
        return self.__datadoc.find('numeric-parameters/external-solver-library').text
    
    def getNumerOfImages(self):
        return int(self.__datadoc.find('output/number-of-images').text)
    
    def hasImagesDump(self):
        if (self.__datadoc.find('output/dump-images').text).lower() == 'true':
            return True
        else:
            return False

    def hasImagesShow(self):
        if (self.__datadoc.find('output/show-images').text).lower() == 'true':
            return True
        else:
            return False