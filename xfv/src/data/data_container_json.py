# -*- coding: utf-8 -*-
"""
Implementing the DataContainer class
"""

from collections import namedtuple
import json
from pathlib import Path
import sys

from xfv.src.equationsofstate.miegruneisen import MieGruneisen
from xfv.src.plasticitycriterion.vonmises import VonMisesCriterion
from xfv.src.rupturecriterion.halfrodcomparison import HalfRodComparisonCriterion
from xfv.src.rupturecriterion.minimumpressure import MinimumPressureCriterion
from xfv.src.rupturecriterion.damage_criterion import DamageCriterion
from xfv.src.rupturecriterion.maximalstress import MaximalStressCriterion
from xfv.src.rheology.constantshearmodulus import ConstantShearModulus
from xfv.src.rheology.constantyieldstress import ConstantYieldStress
from xfv.src.custom_functions.constant_value import ConstantValue
from xfv.src.custom_functions.march_table import MarchTable
from xfv.src.custom_functions.ramp import Ramp
from xfv.src.custom_functions.two_steps import TwoSteps
from xfv.src.custom_functions.successive_ramp import SuccessiveRamp
from xfv.src.cohesive_model.bilinear_cohesive_law import BilinearCohesiveZoneModel
from xfv.src.cohesive_model.linear_cohesive_law import LinearCohesiveZoneModel
from xfv.src.cohesive_model.trilinear_cohesive_law import TrilinearCohesiveZoneModel
from xfv.src.cohesive_model.progressive_unloading_model import ProgressiveUnloadingModel
from xfv.src.cohesive_model.zero_force_unloading_model import ZeroForceUnloadingModel
from xfv.src.cohesive_model.loss_of_stiffness_unloading_model import LossOfStiffnessUnloadingModel
from xfv.src.utilities.singleton import Singleton

NumericalProps = namedtuple("NumericalProps", ["a_pseudo", "b_pseudo", "cfl", "cfl_pseudo"])

BoundaryConditionsProps = namedtuple("BoundaryConditionsProps", ["left_BC", "right_BC"])

GeometricalProps = namedtuple("GeometricalProps", ["section", "initial_interface_position"])

MaterialProps = namedtuple("MaterialProps", ["initial_values", "constitutive_model",
                                             "failure_model", "damage_model"])

InitialValues = namedtuple("InitialValues", ["velocity_init", "pression_init", "temp_init",
                                             "rho_init", "energie_init",
                                             "yield_stress_init", "shear_modulus_init"])

ConstitutiveModel = namedtuple("ConstitutiveModel", ["eos", "elasticity_model",
                                                     "plasticity_model", "plasticity_criterion"])

FailureModel = namedtuple("FailureModel", ["failure_treatment", "failure_treatment_value",
                                           "type_of_enrichment", "lump_mass_matrix",
                                           "failure_criterion", "failure_criterion_value"])

DamageModel = namedtuple("DamageModel", ["cohesive_model", "name"])


TimeProps = namedtuple("TimeProps", ['initial_time_step', 'final_time', 'is_time_step_constant',
                                     "time_step_reduction_factor_for_failure"])

OutputProps = namedtuple("OutputProps", ['number_of_images', 'dump', 'databases'])

DatabaseProps = namedtuple("DatabaseProps",
                           ["identifier", "path", "time_period", "iteration_period"])


class DataContainerJson(metaclass=Singleton):
    """
    This class provides access to all the data found in the json datafile
    """
    def __init__(self, datafile_path=None):
        """
        Constructor
        """
        datafile_path = Path(datafile_path)
        print("Opening data file : {}".format(datafile_path.resolve()))
        self._datafile_dir = datafile_path.resolve().parent
        with datafile_path.open('r') as fi:
            self.__datadoc = json.load(fi)
        self.numeric = NumericalProps(*self.__fill_in_numerical_props())
        self.geometric = GeometricalProps(*self.__fill_in_geometrical_props())
        self.time = TimeProps(*self.__fill_in_time_props())
        self.output = OutputProps(*self.__fill_in_output_props())
        # self.boundary_condition = BoundaryConditionsProps(*self.__fill_in_bc_props())

    def __fill_in_numerical_props(self):
        """
        Returns the quantities needed to fill numerical properties:
            - coefficient of linear artificial viscosity
            - coefficient of quadratic artificial viscosity
            - CFL coefficient
            - CFL coefficient of artificial viscosity
        """
        params = self.__datadoc['numeric-parameters']
        return (params['linear-pseudo'], params['quadratic-pseudo'],
                params['cfl'], params['cfl-pseudo'])

    def __fill_in_geometrical_props(self):
        """
        Returns the quantities needed to fill geometrical properties
            - area of section of the cell
            - position of the interface between target and projectile
              (2 materials case)
        """
        params = self.__datadoc['geometry']
        section = params['section']
        initial_interface_position = params.get('initial-interface-position', 0.)
        return section, initial_interface_position

    def __fill_in_time_props(self):
        """
        Returns the quantities needed to fill time properties
            - initial time step
            - final time
            - is time step constant
            - time step reducation factor for failure
        """
        params = self.__datadoc['time-management']
        initial_time_step = params['initial-time-step']
        final_time = params['final-time']
        cst_dt = params.get('constant-time-step', False)
        time_step_reduction = params.get('time-step-reduction-factor-for-failure')
        return initial_time_step, final_time, cst_dt, time_step_reduction

    def __fill_in_output_props(self):
        """
        Returns the quantities needed to fill output properties
            - number of images
            - cell / node selected for extraction of time history
            - is display of times figures required?
            - list of output database properties
        """
        params = self.__datadoc['output']
        number_of_images = params['number-of-images']
        dump = params['dump-images']

        # Databases
        db_prop_l = []
        for elem in params['database']:
            identi = elem['identifier']
            database_path = elem['path']
            iteration_period = elem.get('iteration-period')
            time_period = elem.get('time-period')
            db_props = DatabaseProps(identi, database_path, time_period,
                                     iteration_period)
            db_prop_l.append(db_props)
        return number_of_images, dump, db_prop_l


if __name__ == "__main__":
    data = DataContainerJson(datafile_path='XDATA.json')
    print(data.numeric)
    print(data.geometric)
    print(data.time)
    print(data.output)