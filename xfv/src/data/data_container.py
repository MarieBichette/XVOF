#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementing the DataContainer class
"""

from collections import namedtuple
import os.path
import json
import lxml.etree as et

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
from xfv.src.cohesive_model_unloading.progressive_unloading_model import ProgressiveUnloadingModel
from xfv.src.cohesive_model_unloading.zero_force_unloading_model import ZeroForceUnloadingModel
from xfv.src.cohesive_model_unloading.loss_of_stiffness_unloading_model import LossOfStiffnessUnloadingModel
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

OutputProps = namedtuple("OutputProps", ['number_of_images', 'dump', 'cells_numbers',
                                         'nodes_numbers', 'databases'])

DatabaseProps = namedtuple("DatabaseProps",
                           ["identifier", "path", "time_period", "iteration_period",
                            "cell_indexes", "node_indexes"])


class DataContainer(object, metaclass=Singleton):  # pylint: disable=too-many-instance-attributes
    """
    Contains the data read from the datafile
    """

    def __init__(self, datafile_path=None):
        """
        Constructor
        """
        print("Opening data file : {:s}".format(os.path.abspath(datafile_path)))
        self.project_dir = os.path.commonprefix([os.path.abspath(datafile_path),
                                                 os.path.abspath(__file__)])
        self.__datadoc = et.parse(datafile_path)
        self.numeric = NumericalProps(*self.__fill_in_numerical_props())
        self.geometric = GeometricalProps(*self.__fill_in_geometrical_props())
        self.time = TimeProps(*self.__fill_in_time_props())
        self.output = OutputProps(*self.__fill_in_output_props())
        self.boundary_condition = BoundaryConditionsProps(*self.__fill_in_bc_props())

        # test si on a un projectile et une cible ou si on a qu'un bloc matter
        try:
            _ = self.__datadoc.find("matter/projectile/").text
            self.data_contains_a_projectile = True
        except AttributeError:
            self.data_contains_a_projectile = False
        try:
            _ = self.__datadoc.find("matter/target/").text
            self.data_contains_a_target = True
        except AttributeError:
            self.data_contains_a_target = False

        self.material_projectile = None
        if self.data_contains_a_projectile:
            self.material_projectile = MaterialProps(
                *self.__fill_in_material_props("matter/projectile"))

        self.material_target = None
        if self.data_contains_a_target:
            self.material_target = MaterialProps(*self.__fill_in_material_props("matter/target"))

        if not (self.data_contains_a_projectile or self.data_contains_a_target):
            self.material_target = MaterialProps(*self.__fill_in_material_props("matter"))
            # astuce pour initilialiser un projectile fictif et ne pas changer la structure du code
            self.material_projectile = self.material_target

    def __fill_in_bc_props(self):
        """
        :return: the pressure to be applied on the boundary of geometry
        """
        left_boundary = self.__get_boundary_condition_def('boundary-conditions/left-boundary/')
        right_boundary = self.__get_boundary_condition_def('boundary-conditions/right-boundary/')
        return left_boundary, right_boundary

    def __fill_in_numerical_props(self):
        """
        :return: the numerical properties :
            - linear and quadratic artifical viscosity coefficients;
            - CFL and CFL for pseudo;
        :rtype: tuple(float, float, float, int)
        """
        b_pseudo = float(self.__datadoc.find('numeric-parameters/linear-pseudo').text)
        a_pseudo = float(self.__datadoc.find('numeric-parameters/quadratic-pseudo').text)
        cfl = float(self.__datadoc.find('numeric-parameters/cfl').text)
        cfl_pseudo = float(self.__datadoc.find('numeric-parameters/cfl-pseudo').text)
        return a_pseudo, b_pseudo, cfl, cfl_pseudo

    def __fill_in_geometrical_props(self):
        """
        :return: the geometric properties:
            - area of the cell;
            - position of the interface between target and projectile for 2 materials simulations
        :rtype: tuple(float, float)
        """
        section = float(self.__datadoc.find('geometry/section').text)
        initial_interface_position = 0.
        try:
            initial_interface_position = float(
                self.__datadoc.find('geometry/initial-interface-position').text)
        except AttributeError:
            pass

        return section, initial_interface_position

    def __check_material_consistency(self, material="matter/target"):
        """
        Check the coherence of the data
        """
        # Vérification de la cohérence du jeu de données :
        material_eos = self.__datadoc.find(material + '/equation-of-state/' + 'coefficients').text
        material_eos = material_eos.split("_")[0]
        material_init = self.__datadoc.find(material + '/initialization/' + 'init-thermo').text
        material_init = material_init.split("_")[0]
        material_rheology = None
        try:
            material_rheology = self.__datadoc.find(material + '/rheology/' + 'coefficients').text
            material_rheology = material_rheology.split("_")[0]
        except AttributeError:
            pass

        if not material_eos == material_init:
            raise ValueError("""Incohérence dans le jeu de données. """
                             """Les matériaux renseignés pour l'init ({:}) et
                                l'eos ({:}) sont différents""".format(
                                    material_init, material_eos))
        if material_rheology is not None and not material_eos == material_rheology:
            raise ValueError("""Incohérence dans le jeu de données. Les matériaux renseignés
                                 pour l'eos ({:}) et la rheologie ({:}) sont différents""".format(
                                     material_eos, material_rheology))

    def __fill_in_material_props(self, material="matter/target"):
        """
        :return: the material properties:
            - initialization pressure, temperature, density and internal energy,
              yield stress and shear modulus
            - constitutive model : eos, elasticity model (None si pas d'élasticité),
                                   plasticity model(None par défault),
                                   plasticity criterion
            - failure treatment : bool to activate failure, failure treatment method,
                                  failure treatment value (position of discontinuity
                                  in failured element or imposed pressure value) ,
                                  typeof enrichment chosen, mass lumping technique
            - damage treatment :

        """
        init = InitialValues(*self.__get_initial_values(material))

        behavior = ConstitutiveModel(
            self.__get_equation_of_state_props(material),
            *self.__get_rheology_props(material))

        failure_treatment, failure_treatment_value, type_of_enrichment, lump_mass_matrix = \
            self.__get_failure_props(material)
        failure_criterion, failure_criterion_value = self.__get_failure_criterion_props(material)

        failure = FailureModel(failure_treatment, failure_treatment_value, type_of_enrichment,
                               lump_mass_matrix, failure_criterion, failure_criterion_value)

        if failure.failure_treatment is not None and failure_criterion is None:
            raise ValueError("Failure criterion expected. "
                             "No failure criterion is specified or "
                             "specified criterion not understood")

        damage = DamageModel(*self.__get_damage_props(material))

        if failure.failure_treatment is not None and damage.cohesive_model is not None:
            if (failure_criterion_value != damage.cohesive_model.cohesive_strength and
                    isinstance(failure_criterion, MaximalStressCriterion)):
                print("Failure criterion value and cohesive strength have different value. " \
                      "This may result in errors in the future")

        return init, behavior, failure, damage

    def __fill_in_time_props(self):
        """
        :return: time properties :
            - initial time step
            - final time
            - is time step constant
        :rtype: tuple(float, float, bool)
        """
        initial_time_step = float(self.__datadoc.find('time-management/initial-time-step').text)
        final_time = float(self.__datadoc.find('time-management/final-time').text)
        cst_dt = False
        if self.__datadoc.find('time-management/constant-time-step') is not None:
            cst_dt = self.__datadoc.find('time-management/constant-time-step').text == "True"
        try:
            time_step_reduction = float(
                self.__datadoc.find('time-management/time-step-reduction-factor-for-failure').text)
        except AttributeError:
            time_step_reduction = None
        return initial_time_step, final_time, cst_dt, time_step_reduction

    def __fill_in_output_props(self):
        """
        :return:
            - number of images
            -cell_number / node_number : cell / node selected for extraction of time history
            - is display of times figures required?
            - list of output database properties
        :tuple(int, [int], [int], bool, [DatabaseProps])
        """
        number_of_images = int(self.__datadoc.find('output/number-of-images').text)
        dump = self.__datadoc.find('output/dump-images').text.lower() == "true"

        # todo : virer ce type de sortie
        try:
            cell_numbers = self.__datadoc.find('output/cell-for-time-figure').text
            cell_numbers = cell_numbers.split(',')
        except ValueError:  # la ligne correspondante ne contient pas de cell id (de type int)
            cell_numbers = None
        except AttributeError:  # la ligne correspondante est absente ou comment�e
            cell_numbers = None
        try:
            str_node_numbers = self.__datadoc.find('output/node-for-time-figure').text
            node_numbers = str_node_numbers.split(',')
        except ValueError:  # la ligne correspondante ne contient pas de cell id (de type int)
            node_numbers = None
        except AttributeError:  # la ligne correspondante est absente ou comment�e
            node_numbers = None
        #end todo

        # Databases
        db_prop_l = []
        for elem in self.__datadoc.iterfind('output/database'):
            identi = elem.find('identifier').text
            database_path = elem.find('path').text
            iteration_period, time_period = None, None
            cell_indexes, node_indexes = None, None
            if elem.find('iteration-period') is not None:
                iteration_period = int(elem.find('iteration-period').text)
            else:
                time_period = float(elem.find('time-period').text)
            # indices spécifiques : pas implémenté et pose pb pour reconstruction des champs
            # if elem.find('cell-indexes') is not None:
            #     cell_indexes = [int(ind) for ind in elem.find('cell-indexes').text.split(',')]
            # if elem.find('node-indexes') is not None:
            #     node_indexes = [int(ind) for ind in elem.find('node-indexes').text.split(',')]
            db_props = DatabaseProps(identi, database_path, time_period,
                                     iteration_period, cell_indexes, node_indexes)
            db_prop_l.append(db_props)
        return number_of_images, dump, cell_numbers, node_numbers, db_prop_l

    def __get_initial_values(self, matter):
        """
        Reads the XDATA file and find the initialization quantities for the material matter
        :param matter: material to be considered
        :type : string
        :return: initial thermodynamical quantities
        :rtype: float x5
        """
        repertoire_base = matter + '/initialization/'
        velocity = float(self.__datadoc.find(repertoire_base + 'initial-velocity').text)

        json_name = self.__datadoc.find(repertoire_base + 'init-thermo').text
        json_path = os.path.join(self.project_dir, "data/CONSTITUTIVE_MODEL/" + json_name)
        with open(json_path, 'r') as json_fid:
            coef = json.load(json_fid)
            coef = coef["InitThermo"]
            density = float(coef["initial_density"])
            pressure = float(coef["initial_pressure"])
            temperature = float(coef["initial_temperature"])
            internal_energy = float(coef["initial_internal_energy"])

        yield_stress, shear_modulus = self.__get_yield_stress_and_shear_modulus(matter)
        return (velocity, pressure, temperature, density, internal_energy,
                yield_stress, shear_modulus)

    def __get_equation_of_state_props(self, matter):
        """
        Creates the equation of state with XDATA input parameters for the material matter
        :param matter: material to be considered
        :type string
        :return: the equation of state
        :rtype EquationOfState
        """
        repertoire_base = matter + '/equation-of-state/'
        if self.__datadoc.find(repertoire_base + 'name').text == 'Mie-Gruneisen':
            json_name = self.__datadoc.find(repertoire_base + 'coefficients').text
            json_path = os.path.join(self.project_dir, "data/CONSTITUTIVE_MODEL/" + json_name)
            with open(json_path, 'r') as json_fid:
                coef = json.load(json_fid)
                coef = coef["MieGruneisen"]
                # Lecture des paramètres
                params_key = ("ref_sound_velocity", "s1", "s2", "s3", "ref_density",
                              "coefficient_gruneisen", "param_b", "ref_internal_energy")
                params = [float(coef[p]) for p in params_key]
            # Création de l'équation d'état
            return MieGruneisen(*params)
        raise NotImplementedError("Only MieGruneisen equation of state is implemented for now")

    def __get_failure_props(self, matter):
        """
        Creates the failure model with XDATA input parameters for the material matter
        :param matter: material to be considered
        :type string
        :return: has_failure : booleen pour activer la rupture
                failure_model : mod�le de traitement de la rupture
                failure_treatment_value : position of disconitnuity in cracked element
                                          or imposed pressure
                type_of_enrichment : Hansbo
                lump_mass_matrix : lumping strategy
        :rtype bool, str, float, str, str
        """

        repertoire_base = matter + "/failure/failure-treatment/"

        failure_treatment = None
        try:
            failure_treatment = str(self.__datadoc.find(repertoire_base + 'name').text)
        except AttributeError:
            pass

        if failure_treatment is not None:
            if failure_treatment not in ["ImposedPressure", "Enrichment"]:
                raise ValueError("Only 'ImposedPressure' or 'Enrichment' are possible values")
            failure_treatment_value = float(self.__datadoc.find(repertoire_base + 'value').text)

        if failure_treatment == "Enrichment":
            type_of_enrichment = str(
                self.__datadoc.find(repertoire_base + 'type-of-enrichment').text)
            if type_of_enrichment not in ['Hansbo'] and failure_treatment == "Enrichment":
                raise ValueError("Unknown enrichment type. = Hansbo")
            else:
                lump_mass_matrix = str(
                    self.__datadoc.find(repertoire_base + 'lump-mass-matrix').text)
                if lump_mass_matrix.lower() not in ["menouillard", "somme", "none"]:
                    print ("""Don't recognize lumping technique. """
                           """Only those lumps are possible : menouillard, somme \n"""
                           """No lumping technique is applied. Mass matrix is consistent""")
        else:
            print("No failure treatment will be applied for {:}".format(matter))
            failure_treatment_value = 0.
            type_of_enrichment = None
            lump_mass_matrix = None

        return failure_treatment, failure_treatment_value, type_of_enrichment, lump_mass_matrix

    def __get_failure_criterion_props(self, matter):
        """
        Creates the failure criterion model with XDATA input parameters for the material matter
        :param matter: material to be considered
        :type string
        :return: failure_criterion : critère de rupture
                failure_criterion_value : valeur seuil
        :rtype Criterion, float
        """
        repertoire_base = matter + "/failure/failure-criterion/"
        failure_criterion = None
        failure_criterion_value = 0.
        try:
            failure_criterion_name = str(self.__datadoc.find(repertoire_base + 'name').text)
            failure_criterion_value = float(self.__datadoc.find(repertoire_base + 'value').text)
            if failure_criterion_name == "MinimumPressure":
                failure_criterion = MinimumPressureCriterion(failure_criterion_value)
            elif failure_criterion_name == "Damage":
                failure_criterion = DamageCriterion(failure_criterion_value)
            elif failure_criterion_name == "HalfRodComparison":
                failure_criterion = HalfRodComparisonCriterion(failure_criterion_value)
            elif failure_criterion_name == "MaximalStress":
                failure_criterion = MaximalStressCriterion(failure_criterion_value)
            else:
                print("Failure criterion is not in MinimumPressure, "
                      "Damage, HalfRodComparison, MaximalStress")
                exit(0)
        except AttributeError:
            pass

        return failure_criterion, failure_criterion_value

    def __get_damage_props(self, matter):
        """
        Creates the damage model with XDATA input parameters for the material matter
        :param matter (string): material to be considered
        :return: cohesive_model : mod�le de loi cohesive (classe cohesive_model)
        :rtype bool, str, float, float
        """
        repertoire_base = matter + '/failure/damage-treatment/cohesive-model/'

        cohesive_model = None
        cohesive_model_name = ""
        try:
            cohesive_model_name = self.__datadoc.find(repertoire_base + 'name').text
            if cohesive_model_name.lower() not in ["linear", "bilinear", "trilinear"]:
                raise ValueError("""Le modéle cohésif doit être parmi linear, """
                                 """bilinear, trilinear""")
        except AttributeError:
            return None, ""
        if cohesive_model_name != "":
            # Données de base :
            cohesive_strength = float(
                self.__datadoc.find(repertoire_base + 'coefficients/cohesive-strength').text)
            critical_separation = float(
                self.__datadoc.find(repertoire_base + 'coefficients/critical-separation').text)

            # Unloading model
            unloading_model = None
            unloading_model_name = self.__datadoc.find(
                repertoire_base + 'unloading-model/name').text
            if unloading_model_name.lower() not in ["zeroforceunloading", "progressiveunloading",
                                                    "lossofstiffnessunloading"]:
                raise ValueError(""" Le modéle de décharge doit être dans """
                                 """zeroforceunloading, progressiveunloading, """
                                 """lossofstiffnessunloading """)

            if unloading_model_name.lower() == "zeroforceunloading":
                slope = float(self.__datadoc.find(repertoire_base + 'unloading-model/slope').text)
                unloading_model = ZeroForceUnloadingModel(slope, cohesive_strength)
            elif unloading_model_name.lower() == "progressiveunloading":
                slope = float(self.__datadoc.find(repertoire_base + 'unloading-model/slope').text)
                unloading_model = ProgressiveUnloadingModel(slope, cohesive_strength)
            elif unloading_model_name.lower() == "lossofstiffnessunloading":
                slope = float(self.__datadoc.find(repertoire_base + 'unloading-model/slope').text)
                unloading_model = LossOfStiffnessUnloadingModel(slope, cohesive_strength)

            if cohesive_model_name.lower() == "bilinear":
                separation_1 = float(
                    self.__datadoc.find(
                        repertoire_base + 'coefficients/separation-at-point-1').text)
                contrainte_1 = float(
                    self.__datadoc.find(repertoire_base + 'coefficients/stress-at-point-1').text)
                cohesive_model = BilinearCohesiveZoneModel(
                    cohesive_strength, separation_1, contrainte_1,
                    critical_separation, unloading_model)

            elif cohesive_model_name.lower() == "trilinear":
                separation_1 = float(self.__datadoc.find(
                    repertoire_base + 'coefficients/separation-at-point-1').text)
                contrainte_1 = float(self.__datadoc.find(
                    repertoire_base + 'coefficients/stress-at-point-1').text)
                separation_2 = float(self.__datadoc.find(
                    repertoire_base + 'coefficients/separation-at-point-2').text)
                contrainte_2 = float(
                    self.__datadoc.find(repertoire_base + 'coefficients/stress-at-point-2').text)
                cohesive_model = TrilinearCohesiveZoneModel(
                    cohesive_strength, separation_1, contrainte_1, separation_2, contrainte_2,
                    critical_separation, unloading_model)
            else:  # loi lineaire
                cohesive_model = LinearCohesiveZoneModel(
                    cohesive_strength, critical_separation, unloading_model)
        print("Applying a damage model for " + matter + " which is : " + cohesive_model_name)
        print("Unloading model is : " + unloading_model_name)
        return cohesive_model, cohesive_model_name

    def __get_yield_stress_and_shear_modulus(self, matter):
        repertoire_base = matter + '/rheology/'
        try:
            json_name = self.__datadoc.find(repertoire_base + 'coefficients').text
        except AttributeError:
            return 0, 0
        json_path = os.path.join(self.project_dir, "data/CONSTITUTIVE_MODEL/" + json_name)
        with open(json_path, 'r') as json_fid:
            coef = json.load(json_fid)
            coef = coef["EPP"]
        yield_stress = float(coef["yield_stress"])
        shear_modulus = float(coef["shear_modulus"])
        return yield_stress, shear_modulus

    def __get_rheology_props(self, matter):
        """
        Reads the elasticity parameters for material matter
        :param matter: material to be considered
        :return: elasticity_model : mod�le d'�lasticit�
                 plasticity_model : mod�le de plasticit�
                 plasticity_criterion : m�thode de calcul pour le crit�re de plasticit�
        :rtype ShearModulus, YieldStress models + VonMisesCriterion
        """
        elasticity_model = None
        plasticity_model = None
        plasticity_criterion = None
        repertoire_base = matter + '/rheology/'
        yield_stress, shear_modulus = self.__get_yield_stress_and_shear_modulus(matter)
        try:
            # Elasticité
            elasticity_model_name = str(self.__datadoc.find(
                repertoire_base + 'elasticity-model').text)
            if elasticity_model_name != "Linear":
                raise ValueError("Model {:} not implemented. Choose Linear".
                                 format(elasticity_model_name))
            elasticity_model = ConstantShearModulus(shear_modulus)
        except AttributeError:
            print("Impossible de construire le modèle d'élasticité")
        try:
            # Plasticité
            plasticity_model_name = str(self.__datadoc.find(
                repertoire_base + 'plasticity-model').text)
            if plasticity_model_name != "EPP":
                raise ValueError("Model {:} not implemented. Choose EPP".
                                 format(plasticity_model_name))
            plasticity_model = ConstantYieldStress(yield_stress)
            plasticity_criterion_name = str(self.__datadoc.find(
                repertoire_base + 'plasticity-criterion').text)
            if plasticity_criterion_name == "VonMises":
                plasticity_criterion = VonMisesCriterion()
        except AttributeError:
            print("Impossible de construire le mod�le de plasticit�")
        return elasticity_model, plasticity_model, plasticity_criterion

    def hasExternalSolver(self):  # pylint: disable=invalid-name
        """
        Returns True if an external solver is required in the data
        """
        return self.__datadoc.find('numeric-parameters/external-solver-library') is not None

    def getExternalSolverPath(self):  # pylint: disable=invalid-name
        """
        Returs the external solver required in the data
        """
        return self.__datadoc.find('numeric-parameters/external-solver-library').text

    def __get_boundary_condition_def(self, info):
        """
        Creates the boundary condition class with pressure law
        :return: pressure_law
        """
        bc_law = None
        type_bc = str(self.__datadoc.find(info + 'type').text)

        if type_bc.lower() == "velocity":
            class_name = str(self.__datadoc.find(info + 'bc-law').text)

            if class_name.lower() == "constant":
                value = float(self.__datadoc.find(info + 'value').text)
                bc_law = ConstantValue(value)
                bc_law.register_velocity()
            else:
                raise ValueError("""Mauvais type de loi de pression en entrée """
                                 """pour {:} avec un type {:}."""
                                 """Les possibilit�s sont : (Constant)""".format(info, type_bc))

        elif type_bc.lower() == "pressure":
            class_name = str(self.__datadoc.find(info + 'bc-law').text)

            if class_name.lower() == "constant":
                value = float(self.__datadoc.find(info + 'value').text)
                bc_law = ConstantValue(value)

            elif class_name.lower() == "twostep":
                params = [float(self.__datadoc.find(info + tag).text)
                          for tag in ('value1', 'value2', 'time-activation')]
                bc_law = TwoSteps(*params)

            elif class_name.lower() == "ramp":
                params = [float(self.__datadoc.find(info + tag).text)
                          for tag in ('value1', 'value2',
                                      'time-activation-value-1',
                                      'time-activation-value-2')
                         ]
                bc_law = Ramp(*params)

            elif class_name == "marchtable":
                _file = str(self.__datadoc.find(info + 'value').text)
                bc_law = MarchTable(_file)

            elif class_name == "creneauramp":
                f_ramp_tags = ('initial-value', 'plateau-value',
                               'start-first-ramp-time', 'reach-value2-time')
                s_ramp_tags = ('plateau-value', 'end-value',
                               'start-second-ramp-time', 'reach-value3-time')
                f_ramp_params = [float(self.__datadoc.find(info + tag).text)
                                 for tag in f_ramp_tags]
                s_ramp_params = [float(self.__datadoc.find(info + tag).text)
                                 for tag in s_ramp_tags]
                first_ramp = Ramp(*f_ramp_params)
                second_ramp = Ramp(*s_ramp_params)
                bc_law = SuccessiveRamp(first_ramp, second_ramp)
            else:
                raise ValueError("""Mauvais type de loi de pression en entrée """
                                 """pour {:} avec un type {:}."""
                                 """Les possibilités sont : """
                                 """(Constant|TwoStep|Ramp|MarchTable|CreaneauRamp)"""
                                 .format(info, type_bc))

            bc_law.register_pressure()
        else:
            raise ValueError
        return bc_law
