#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Implementing the DataContainer class
"""

from collections import namedtuple
import os.path
import json
import lxml.etree as et

from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.plasticitycriterion.vonmises import VonMisesCriterion
from xvof.rupturecriterion.halfrodcomparison import HalfRodComparisonCriterion
from xvof.rupturecriterion.minimumpressure import MinimumPressureCriterion
from xvof.rupturecriterion.damage_criterion import DamageCriterion
from xvof.rupturecriterion.maximalstress import MaximalStressCriterion
from xvof.rheology.constantshearmodulus import ConstantShearModulus
from xvof.rheology.constantyieldstress import ConstantYieldStress
from xvof.boundary_condition.constantvelocity import ConstantVelocity
from xvof.boundary_condition.constantpressure import ConstantPressure
from xvof.boundary_condition.march_table import MarchTablePressure
from xvof.boundary_condition.ramppressure import RampPressure
from xvof.boundary_condition.twostepspressure import TwoStepsPressure
from xvof.boundary_condition.creneau_ramp_pressure import CreneauRampPressure
from xvof.cohesive_model.bilinear_cohesive_law import BilinearCohesiveZoneModel
from xvof.cohesive_model.linear_cohesive_law import LinearCohesiveZoneModel
from xvof.cohesive_model.trilinear_cohesive_law import TrilinearCohesiveZoneModel
from xvof.cohesive_model.progressive_unloading_model import ProgressiveUnloadingModel
from xvof.cohesive_model.zero_force_unloading_model import ZeroForceUnloadingModel
from xvof.cohesive_model.loss_of_stiffness_unloading_model import LossOfStiffnessUnloadingModel
from xvof.utilities.singleton import Singleton

numerical_props = namedtuple("numerical_props", ["a_pseudo", "b_pseudo", "cfl", "cfl_pseudo"])

boundary_conditions_props = namedtuple("boundary_conditions_props", ["left_BC", "right_BC"])

geometrical_props = namedtuple("geometrical_props", ["section", "initial_interface_position"])

material_props = namedtuple("material_props", ["initial_values", "constitutive_model", "failure_model", "damage_model"])

initial_values = namedtuple("initial_values", ["velocity_init", "pression_init", "temp_init", "rho_init", "energie_init",
                                               "yield_stress_init", "shear_modulus_init"])

constitutive_model = namedtuple("constitutive_model", ["eos", "elasticity_model",
                                                       "plasticity_model", "plasticity_criterion"])

failure_model = namedtuple("failure_model", ["failure_treatment", "failure_treatment_value",
                                               "type_of_enrichment", "lump_mass_matrix",
                                               "failure_criterion", "failure_criterion_value"])

damage_model = namedtuple("damage_model", ["cohesive_model", "name"])


time_props = namedtuple("time_props", ['initial_time_step', 'final_time', 'is_time_step_constant',
                                       "time_step_reduction_factor_for_failure"])

output_props = namedtuple("output_props", ['number_of_images', 'dump', 'cells_numbers', 'nodes_numbers', 'databases'])

database_props = namedtuple("database_props",
                            ["identifier", "path", "time_period", "iteration_period", "cell_indexes", "node_indexes"])


class DataContainer(object):
    """
    Contains the data read from the datafile
    """
    __metaclass__ = Singleton

    def __init__(self, datafile_path):
        """
        Constructor
        """
        print "Opening data file : {:s}".format(os.path.abspath(datafile_path))
        self.datafile_path = datafile_path
        self.project_dir = os.path.commonprefix([os.path.abspath(self.datafile_path),
                                                 os.path.abspath(__file__)])
        # self.project_dir = os.path.split(os.path.dirname(os.path.abspath(self.datafile_path)))[0]
        self.__datadoc = et.parse(datafile_path)
        self.numeric = numerical_props(*self.__fillInNumericalProperties())
        self.geometric = geometrical_props(*self.__fillInGeometricalProperties())
        self.time = time_props(*self.__fillInTimeProperties())
        self.output = output_props(*self.__fillInOutputProperties())
        self.boundary_condition = boundary_conditions_props(*self.__fillInBCProperties())

        # test si on a un projectile et une cible ou si on a qu'un bloc matter
        try:
            blabla = self.__datadoc.find("matter/projectile/").text
            self.data_contains_a_projectile = True
        except AttributeError:
            self.data_contains_a_projectile = False
        try:
            bla = self.__datadoc.find("matter/target/").text
            self.data_contains_a_target = True
        except AttributeError:
            self.data_contains_a_target = False

        if self.data_contains_a_projectile:
            self.material_projectile = material_props(*self.__fillInMaterialProperties("matter/projectile"))
        if self.data_contains_a_target:
            self.material_target = material_props(*self.__fillInMaterialProperties("matter/target"))

        if not (self.data_contains_a_projectile or self.data_contains_a_target):
            self.material_target = material_props(*self.__fillInMaterialProperties("matter"))
            self.material_projectile = self.material_target  # astuce pour initilialiser un projectile fictif et ne pas changer la structure du code

    def __fillInBCProperties(self):
        """
        :return: the pressure to be applied on the boundary of geometry
        """
        left_boundary = self.__getBoundaryConditionDefiniton('boundary-conditions/left-boundary/')
        right_boundary = self.__getBoundaryConditionDefiniton('boundary-conditions/right-boundary/')
        return left_boundary, right_boundary

    def __fillInNumericalProperties(self):
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

    def __fillInGeometricalProperties(self):
        """
        :return: the geometric properties:
            - area of the cell;
            - position of the interface between target and projectile for 2 materials simulations
        :rtype: tuple(float, float)
        """
        section = float(self.__datadoc.find('geometry/section').text)
        initial_interface_position = 0.
        try:
            initial_interface_position = float(self.__datadoc.find('geometry/initial-interface-position').text)
        except AttributeError:
            pass

        return section, initial_interface_position

    def __fillInMaterialProperties(self, material="matter/target"):
        """
        :return: the material properties:
            - initialization pressure, temperature, density and internal energy, yield stress and shear modulus
            - constitutive model : eos, elasticity model (None si pas d'élasticité), plasticity model(None par défault),
                                    plasticity criterion
            - failure treatment : bool to activate failure, failure treatment method,
                                  failure treatment value (position of discontinuity
                                  in failured element or imposed pressure value) , typeof enrichment chosen, mass lumping technique
            - damage treatment :

        """
        # Vérification de la cohérence du jeu de données :
        mat_eos = self.__datadoc.find(material + '/equation-of-state/' + 'coefficients').text
        material_eos = mat_eos.split("_")[0]
        mat_init = self.__datadoc.find(material + '/initialization/' + 'init-thermo').text
        material_init = mat_init.split("_")[0]
        material_rheology = None
        try:
            mat_rheo = self.__datadoc.find(material + '/rheology/' + 'coefficients').text
            material_rheology = mat_rheo.split("_")[0]
        except AttributeError:
            pass

        if not(material_eos == material_init):
            raise ValueError("""Incohérence dans le jeu de donnée. Les matériaux renseignés pour l'init ({:}) et
            l'eos ({:}) sont différents""".format(material_init, material_eos))
        if material_rheology is not None:
            if not(material_eos == material_rheology):
                raise ValueError("""Incohérence dans le jeu de donnée. Les matériaux renseignés
                l'eos ({:}) et la rheologie ({:}) sont différents""".format(material_eos, material_rheology))

        # Grandeurs initiales -----------------------------
        init = initial_values(*self.__getInitialValues(material))

        # Constitutive model ----------------------------------
        eos = self.__getEquationOfStateProperties(material)
        elasticity_model, plasticity_model, plasticity_criterion = self.__getRheologyProperties(material)
        behavior = constitutive_model(eos, elasticity_model, plasticity_model, plasticity_criterion)

        # failure treatment ------------------------------------
        failure_treatment, failure_treatment_value, type_of_enrichment, lump_mass_matrix = \
            self.__getFailureProperties(material)
        failure_criterion, failure_criterion_value = self.__getFailureCriterionProperties(material)

        failure = failure_model(failure_treatment, failure_treatment_value, type_of_enrichment,
                                lump_mass_matrix, failure_criterion, failure_criterion_value)

        if failure.failure_treatment is not None and failure_criterion is None:
            raise ValueError("Failure criterion expected. "
                             "No failure criterion is specified or specified criterion not understood")

        # Damage Law-----------------------------------
        damage = damage_model(*self.__getDamageProperties(material))

        if failure.failure_treatment is not None and damage.cohesive_model is not None:
            if failure_criterion_value != damage.cohesive_model.cohesive_strength and \
                            type(failure_criterion) == MaximalStressCriterion:
                print "Failure criterion value and cohesive strength have different value. " \
                      "This may result in errors in the future"

        # Return ---------------------------------------------------
        return (init, behavior, failure, damage)

    def __fillInTimeProperties(self):
        """
        :return: time properties :
            - initial time step
            - final time
            - is time step constant
        :rtype: tuple(float, float, bool)
        """
        initial_time_step = float(self.__datadoc.find('time-management/initial-time-step').text)
        final_time = float(self.__datadoc.find('time-management/final-time').text)
        if self.__datadoc.find('time-management/constant-time-step') is not None:
            cst_dt = True if self.__datadoc.find('time-management/constant-time-step').text == "True" else False
        else:
            cst_dt = False
        try:
            time_step_reduction = float(self.__datadoc.find('time-management/time-step-reduction-factor-for-failure').text)
        except AttributeError:
            time_step_reduction = None
        return initial_time_step, final_time, cst_dt, time_step_reduction

    def __fillInOutputProperties(self):
        """
        :return:
            - number of images
            -cell_number / node_number : cell / node selected for extraction of time history
            - is display of times figures required?
            - list of output database properties
        :tuple(int, [int], [int], bool, [database_props])
        """
        number_of_images = int(self.__datadoc.find('output/number-of-images').text)
        dump = True if self.__datadoc.find('output/dump-images').text.lower() == "true" else False

        # todo : virer ce type de sortie
        try:
            str_cell_numbers = self.__datadoc.find('output/cell-for-time-figure').text
            cell_numbers = str_cell_numbers.split(',')
        except ValueError:  # la ligne correspondante ne contient pas de cell id (de type int)
            cell_numbers = None
        except AttributeError:  # la ligne correspondante est absente ou commentée
            cell_numbers = None
        try:
            str_node_numbers = self.__datadoc.find('output/node-for-time-figure').text
            node_numbers = str_node_numbers.split(',')
        except ValueError:  # la ligne correspondante ne contient pas de cell id (de type int)
            node_numbers = None
        except AttributeError:  # la ligne correspondante est absente ou commentée
            node_numbers = None
        #end todo

        # Databases
        db_prop_l = []
        for el in self.__datadoc.iterfind('output/database'):
            identi = el.find('identifier').text
            database_path = el.find('path').text
            iteration_period, time_period = None, None
            cell_indexes, node_indexes = None, None
            if el.find('iteration-period') is not None:
                iteration_period = int(el.find('iteration-period').text)
            else:
                time_period = float(el.find('time-period').text)
            # indices spécifiques : pas implémenté et pose pb pour reconstruction des champs
            # if el.find('cell-indexes') is not None:
            #     cell_indexes = [int(ind) for ind in el.find('cell-indexes').text.split(',')]
            # if el.find('node-indexes') is not None:
            #     node_indexes = [int(ind) for ind in el.find('node-indexes').text.split(',')]
            db_props = database_props(identi, database_path, time_period, iteration_period, cell_indexes, node_indexes)
            db_prop_l.append(db_props)
        return number_of_images, dump, cell_numbers, node_numbers, db_prop_l

    def __getInitialValues(self, matter):
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
        json_path = os.path.join(self.project_dir, "CONSTITUTIVE_MODEL/" + json_name)
        with open(json_path, 'r') as json_fid:
            coef = json.load(json_fid)
            coef = coef["InitThermo"]
            density = float(coef["initial_density"])
            pressure = float(coef["initial_pressure"])
            temperature = float(coef["initial_temperature"])
            internal_energy = float(coef["initial_internal_energy"])

        yield_stress, shear_modulus = self.__getYieldStressAndShearModulusFromData(matter)
        return velocity, pressure, temperature, density, internal_energy, yield_stress, shear_modulus

    def __getEquationOfStateProperties(self, matter):
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
            json_path = os.path.join(self.project_dir, "CONSTITUTIVE_MODEL/" + json_name)
            with open(json_path, 'r') as json_fid:
                coef = json.load(json_fid)
                coef = coef["MieGruneisen"]
                # Lecture des paramètres
                czero = float(coef["ref_sound_velocity"])
                S1 = float(coef["s1"])
                S2 = float(coef["s2"])
                S3 = float(coef["s3"])
                rhozero = float(coef["ref_density"])
                grunzero = float(coef["coefficient_gruneisen"])
                b = float(coef["param_b"])
                ezero = float(coef["ref_internal_energy"])
            # Création de l'équation d'état
            eos = MieGruneisen(czero, S1, S2, S3, rhozero, grunzero, b, ezero)
        else:
            raise NotImplementedError("Only MieGruneisen equation of state is implemented for now")
        return eos

    def __getFailureProperties(self, matter):
        """
        Creates the failure model with XDATA input parameters for the material matter
        :param matter: material to be considered
        :type string
        :return: has_failure : booleen pour activer la rupture
                failure_model : modèle de traitement de la rupture
                failure_treatment_value : position of disconitnuity in cracked element or imposed pressure
                type_of_enrichment : Moes or Hansbo
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
            type_of_enrichment = str(self.__datadoc.find(repertoire_base + 'type-of-enrichment').text)
            if type_of_enrichment not in ['Moes', 'Hansbo'] and failure_treatment == "Enrichment":
                raise ValueError("Unknown enrichment type. Possibilities are Moes and Hansbo")
            else:
                lump_mass_matrix = str(self.__datadoc.find(repertoire_base + 'lump-mass-matrix').text)
                if lump_mass_matrix.lower() not in ["menouillard", "somme", "none"]:
                    print "Don't recognize lumping technique. Only those lumps are possible : menouillard, somme \n" \
                          "No lumping technique is applied. Mass matrix is consistent"
        else:
            print ("No failure treatment will be applied for {:}".format(matter))
            failure_treatment_value = 0.
            type_of_enrichment = None
            lump_mass_matrix = None

        return failure_treatment, failure_treatment_value, type_of_enrichment, lump_mass_matrix
    
    def __getFailureCriterionProperties(self, matter):
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
                print("Failure criterion is not in MinimumPressure, Damage, HalfRodComparison, MaximalStress")
                exit(0)
        except AttributeError:
            pass

        return failure_criterion, failure_criterion_value
        
    def __getDamageProperties(self, matter):
        """
        Creates the damage model with XDATA input parameters for the material matter
        :param matter (string): material to be considered
        :return: cohesive_model : modèle de loi cohesive (classe cohesive_model)
        :rtype bool, str, float, float
        """
        repertoire_base = matter + '/failure/damage-treatment/cohesive-model/'

        cohesive_model = None
        cohesive_model_name = ""
        try:
            cohesive_model_name = self.__datadoc.find(repertoire_base + 'name').text
            if cohesive_model_name.lower() not in ["linear", "bilinear", "trilinear"]:
                raise ValueError("""Le modèle cohésif doit être parmi linear, bilinear, trilinear""")
        except AttributeError:
            return None, ""
        if cohesive_model_name != "":
            # Données de base :
            cohesive_strength = float(self.__datadoc.find(repertoire_base + 'coefficients/cohesive-strength').text)
            critical_separation = float(self.__datadoc.find(repertoire_base + 'coefficients/critical-separation').text)

            # Unloading model
            unloading_model = None
            unloading_model_name = self.__datadoc.find(repertoire_base + 'unloading-model/name').text
            if unloading_model_name.lower() not in ["zeroforceunloading", "progressiveunloading",
                                                    "lossofstiffnessunloading"]:
                raise ValueError(""" Le modèle de décharge doit être dans "zeroforceunloading", "progressiveunloading",
                "lossofstiffnessunloading" """)

            if unloading_model_name.lower() == "zeroforceunloading":
                slope = float(self.__datadoc.find(repertoire_base + 'unloading-model/slope').text)
                unloading_model = ZeroForceUnloadingModel(slope, cohesive_strength)
            elif unloading_model_name.lower() == "progressiveunloading":
                slope = float(self.__datadoc.find(repertoire_base + 'unloading-model/slope').text)
                unloading_model = ProgressiveUnloadingModel(slope, cohesive_strength)
            elif unloading_model_name.lower() == "lossofstiffnessunloading":
                slope = float(self.__datadoc.find(repertoire_base + 'unloading-model/slope').text)
                unloading_model = LossOfStiffnessUnloadingModel(slope,cohesive_strength)

            if cohesive_model_name.lower() == "bilinear":
                separation_1 = float(self.__datadoc.find(repertoire_base + 'coefficients/separation-at-point-1').text)
                contrainte_1 = float(self.__datadoc.find(repertoire_base + 'coefficients/stress-at-point-1').text)
                cohesive_model = BilinearCohesiveZoneModel(cohesive_strength, separation_1, contrainte_1,
                                                           critical_separation, unloading_model)

            elif cohesive_model_name.lower() == "trilinear":
                separation_1 = float(self.__datadoc.find(repertoire_base + 'coefficients/separation-at-point-1').text)
                contrainte_1 = float(self.__datadoc.find(repertoire_base + 'coefficients/stress-at-point-1').text)
                separation_2 = float(self.__datadoc.find(repertoire_base + 'coefficients/separation-at-point-2').text)
                contrainte_2 = float(self.__datadoc.find(repertoire_base + 'coefficients/stress-at-point-2').text)
                cohesive_model = TrilinearCohesiveZoneModel(cohesive_strength, separation_1, contrainte_1,
                                                           separation_2, contrainte_2,
                                                           critical_separation, unloading_model)
            else:  # loi lineaire
                cohesive_model = LinearCohesiveZoneModel(cohesive_strength, critical_separation, unloading_model)
        print "Applying a damage model for " + matter + " which is : " + cohesive_model_name
        print "Unloading model is : " + unloading_model_name
        return cohesive_model, cohesive_model_name

    def __getYieldStressAndShearModulusFromData(self, matter):
        repertoire_base = matter + '/rheology/'
        try:
            json_name = self.__datadoc.find(repertoire_base + 'coefficients').text
        except AttributeError:
            return 0, 0
        json_path = os.path.join(self.project_dir, "CONSTITUTIVE_MODEL/" + json_name)
        with open(json_path, 'r') as json_fid:
            coef = json.load(json_fid)
            coef = coef["EPP"]
        yield_stress = float(coef["yield_stress"])
        shear_modulus = float(coef["shear_modulus"])
        return yield_stress, shear_modulus

    def __getRheologyProperties(self, matter):
        """
        Reads the elasticity parameters for material matter
        :param matter: material to be considered
        :return: elasticity_model : modèle d'élasticité
                 plasticity_model : modèle de plasticité
                 plasticity_criterion : méthode de calcul pour le critère de plasticité
        :rtype ShearModulus, YieldStress models + VonMisesCriterion
        """
        elasticity_model = None
        plasticity_model = None
        plasticity_criterion = None
        repertoire_base = matter + '/rheology/'
        yield_stress, shear_modulus = self.__getYieldStressAndShearModulusFromData(matter)
        try:
            # Elasticité
            elasticity_model_name = str(self.__datadoc.find(repertoire_base + 'elasticity-model').text)
            if elasticity_model_name != "Linear":
                raise ValueError("Model {:} not implemented. Choose Linear".format(elasticity_model_name))
            elasticity_model = ConstantShearModulus(shear_modulus)
        except AttributeError:
            print "Impossible de construire le modèle d'élasticité"
            pass
        try:
            # Plasticité
            plasticity_model_name = str(self.__datadoc.find(repertoire_base + 'plasticity-model').text)
            if plasticity_model_name != "EPP":
                raise ValueError("Model {:} not implemented. Choose EPP".format(plasticity_model_name))
            plasticity_model = ConstantYieldStress(yield_stress)
            plasticity_criterion_name = str(self.__datadoc.find(repertoire_base + 'plasticity-criterion').text)
            if plasticity_criterion_name == "VonMises":
                plasticity_criterion = VonMisesCriterion()
        except AttributeError:
            print "Impossible de construire le modèle de plasticité"
            pass
        return elasticity_model, plasticity_model, plasticity_criterion

    def __getMaterialDelimiters(self, matter):
        repertoire_base = matter + '/init-geometrical-limit/'
        left = float(self.__datadoc.find(repertoire_base + 'left').text)
        right = float(self.__datadoc.find(repertoire_base + 'right').text)
        return left, right

    def hasExternalSolver(self):
        if self.__datadoc.find('numeric-parameters/external-solver-library') is not None:
            return True
    
    def getExternalSolverPath(self):
        return self.__datadoc.find('numeric-parameters/external-solver-library').text

    def __getBoundaryConditionDefiniton(self, info):
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
                bc_law = ConstantVelocity(value)
            else:
                raise ValueError("""Mauvais type de loi de pression en entrée pour {:} avec un type {:}.
                                 Les possibilités sont : (Constant)""".format(info, type_bc))

        elif type_bc.lower() == "pressure":
            class_name = str(self.__datadoc.find(info + 'bc-law').text)

            if class_name.lower() == "constant":
                value = float(self.__datadoc.find(info + 'value').text)
                bc_law = ConstantPressure(value)

            elif class_name.lower() == "twostep":
                value1 = float(self.__datadoc.find(info + 'value1').text)  # 1st value
                value2 = float(self.__datadoc.find(info + 'value2').text)  # 2nd value
                time = float(self.__datadoc.find(info + 'time-activation').text)  # time
                bc_law = TwoStepsPressure(value1, value2, time)

            elif class_name.lower() == "ramp":
                value1 = float(self.__datadoc.find(info + 'value1').text)
                value2 = float(self.__datadoc.find(info + 'value2').text)
                time1 = float(self.__datadoc.find(info + 'time-activation-value1').text)
                time2 = float(self.__datadoc.find(info + 'time-activation-value2').text)
                # et une rampe entre time-activation-value1 et time-activation-value2
                bc_law = RampPressure(value1, value2, time1, time2)

            elif class_name == "marchtable":
                file = str(self.__datadoc.find(info + 'value').text)
                bc_law = MarchTablePressure(file)

            elif class_name == "creneauramp":
                value1 = float(self.__datadoc.find(info + 'initial-value').text)
                value2 = float(self.__datadoc.find(info + 'plateau-value').text)
                value3 = float(self.__datadoc.find(info + 'end-value').text)
                time1 = float(self.__datadoc.find(info + 'start-first-ramp-time').text)
                time2 = float(self.__datadoc.find(info + 'reach-value2-time').text)
                time3 = float(self.__datadoc.find(info + 'start-second-ramp-time').text)
                time4 = float(self.__datadoc.find(info + 'reach-value3-time').text)
                bc_law = CreneauRampPressure(value1, value2, value3, time1, time2, time3, time4)
            else:
                raise ValueError("""Mauvais type de loi de pression en entrée pour {:} avec un type {:}.
                                Les possibilités sont : (Constant|TwoStep|Ramp|MarchTable|CreaneauRamp)"""
                                 .format(info, type_bc))
        else:
            raise(ValueError, """Mauvais type de CL : les possibilités sont (pressure | velocity)""")
        return bc_law