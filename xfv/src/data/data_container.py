# -*- coding: utf-8 -*-
"""
Implementing the DataContainer class
"""
import json
from pathlib import Path
from typing import Dict, List, NamedTuple, Tuple, Optional, Union

from xfv.src.utilities.singleton import Singleton
from xfv.src.data.user_defined_functions_props import (
    UserDefinedFunctionPropsType, ConstantValueFunctionProps, TwoStepsFunctionProps,
    RampFunctionProps, MarchTableFunctionProps, SuccessiveRampFunctionProps)
from xfv.src.data.cohesive_model_props import (CohesiveZoneModelProps,
                                               LinearCohesiveZoneModelProps,
                                               BilinearCohesiveZoneModelProps,
                                               TrilinearCohesiveZoneModelProps)
from xfv.src.data.unloading_model_props import (UnloadingModelProps,
                                                ConstantStiffnessUnloadingProps,
                                                LossOfStiffnessUnloadingProps)
from xfv.src.data.contact_props import (ContactProps, PenaltyContactProps)
from xfv.src.data.equation_of_state_props import (EquationOfStateProps, MieGruneisenProps)
from xfv.src.data.yield_stress_props import (YieldStressProps, ConstantYieldStressProps)
from xfv.src.data.shear_modulus_props import (ShearModulusProps, ConstantShearModulusProps)
from xfv.src.data.plasticity_criterion_props import (PlasticityCriterionProps,
                                                     VonMisesCriterionProps)
from xfv.src.data.rupture_criterion_props import (RuptureCriterionProps,
                                                  DamageCriterionProps,
                                                  HalfRodComparisonCriterionProps,
                                                  MaximalStressCriterionProps,
                                                  MinimumPressureCriterionProps)


class NumericalProps(NamedTuple):  # pylint: disable=missing-class-docstring
    a_pseudo: float
    b_pseudo: float
    cfl: float
    cfl_pseudo: float


class GeometricalProps(NamedTuple):  # pylint: disable=missing-class-docstring
    section: float
    initial_interface_position: float


class TimeProps(NamedTuple):  # pylint: disable=missing-class-docstring
    initial_time_step: float
    final_time: float
    is_time_step_constant: bool
    time_step_reduction_factor_for_failure: Optional[float]


class DatabaseProps(NamedTuple):  # pylint: disable=missing-class-docstring
    identifier: str
    path: str
    time_period: Optional[float]
    iteration_period: Optional[int]


class OutputProps(NamedTuple):  # pylint: disable=missing-class-docstring
    number_of_images: int
    dump: bool
    databases: List[DatabaseProps]


class BoundaryType(NamedTuple):  # pylint: disable=missing-class-docstring
    type_bc: str
    law: UserDefinedFunctionPropsType


class BoundaryConditionsProps(NamedTuple):  # pylint: disable=missing-class-docstring
    left_BC: BoundaryType
    right_BC: BoundaryType


class DamageModelProps(NamedTuple):  # pylint: disable=missing-class-docstring
    cohesive_model: CohesiveZoneModelProps
    name: str


class ContactModelProps(NamedTuple):  # pylint: disable=missing-class-docstring
    contact_model: ContactProps
    name: str


class InitialValues(NamedTuple):  # pylint: disable=missing-class-docstring
    velocity_init: float
    pression_init: float
    temp_init: float
    rho_init: float
    energie_init: float
    yield_stress_init: float
    shear_modulus_init: float


class ConstitutiveModelProps(NamedTuple):  # pylint: disable=missing-class-docstring
    eos: EquationOfStateProps
    elasticity_model: ShearModulusProps
    plasticity_model: YieldStressProps
    plasticity_criterion: PlasticityCriterionProps


class FailureModelProps(NamedTuple):  # pylint: disable=missing-class-docstring
    failure_treatment: str
    failure_treatment_value: float
    type_of_enrichment: str
    lump_mass_matrix: str
    failure_criterion: RuptureCriterionProps
    failure_criterion_value: float


class MaterialProps(NamedTuple):  # pylint: disable=missing-class-docstring
    initial_values: InitialValues
    constitutive_model: ConstitutiveModelProps
    failure_model: FailureModelProps
    damage_model: DamageModelProps
    contact_model: ContactModelProps


class DataContainer(metaclass=Singleton):  # pylint: disable=too-few-public-methods, too-many-instance-attributes
    """
    This class provides access to all the data found in the json datafile
    """
    def __init__(self, datafile_path: Union[str, Path]):
        """
        Constructor
        """
        datafile_path = Path(datafile_path)
        print("Opening data file : {}".format(datafile_path.resolve()))
        self._datafile_dir = datafile_path.resolve().parent
        with datafile_path.open('r') as file_in:
            self.__datadoc = json.load(file_in)
        self.numeric = NumericalProps(*self.__fill_in_numerical_props())
        self.geometric = GeometricalProps(*self.__fill_in_geometrical_props())
        self.time = TimeProps(*self.__fill_in_time_props())
        self.output = OutputProps(*self.__fill_in_output_props())
        self.boundary_condition = BoundaryConditionsProps(*self.__fill_in_bc_props())

        # test si on a un projectile et une cible ou si on a qu'un bloc matter
        matter_data = self.__datadoc['matter']
        self.data_contains_a_projectile = 'projectile' in matter_data.keys()
        self.data_contains_a_target = 'target' in matter_data.keys()

        self.material_projectile = None
        if self.data_contains_a_projectile:
            self.material_projectile = MaterialProps(
                *self.__fill_in_material_props(matter_data['projectile']))

        if self.data_contains_a_target:
            self.material_target = MaterialProps(
                *self.__fill_in_material_props(matter_data['target']))
        else:
            self.material_target = MaterialProps(
                *self.__fill_in_material_props(matter_data))
            self.data_contains_a_target = True

        if not (self.data_contains_a_projectile or self.data_contains_a_target):
            self.material_projectile = self.material_target

    def __fill_in_numerical_props(self) -> Tuple[float, float, float, float]:
        """
        Returns the quantities needed to fill numerical properties:
            - coefficient of linear artificial viscosity
            - coefficient of quadratic artificial viscosity
            - CFL coefficient
            - CFL coefficient of artificial viscosity
        """
        params: Dict[str, float] = self.__datadoc['numeric-parameters']
        return (params['quadratic-pseudo'], params['linear-pseudo'],
                params['cfl'], params['cfl-pseudo'])

    def __fill_in_geometrical_props(self) -> Tuple[float, float]:
        """
        Returns the quantities needed to fill geometrical properties
            - area of section of the cell
            - position of the interface between target and projectile
              (2 materials case)
        """
        params: Dict[str, float] = self.__datadoc['geometry']
        section = params['section']
        initial_interface_position = params.get('initial-interface-position', 0.)
        return section, initial_interface_position

    def __fill_in_time_props(self) -> Tuple[float, float, bool, Optional[float]]:
        """
        Returns the quantities needed to fill time properties
            - initial time step
            - final time
            - is time step constant
            - time step reducation factor for failure
        """
        params = self.__datadoc['time-management']
        initial_time_step: float = params['initial-time-step']
        final_time: float = params['final-time']
        cst_dt: bool = params.get('constant-time-step', False)
        time_step_reduction: Optional[float] = params.get('time-step-reduction-factor-for-failure')
        return initial_time_step, final_time, cst_dt, time_step_reduction

    def __fill_in_output_props(self) -> Tuple[int, bool, List[DatabaseProps]]:
        """
        Returns the quantities needed to fill output properties
            - number of images
            - cell / node selected for extraction of time history
            - is display of times figures required?
            - list of output database properties
        """
        params = self.__datadoc['output']
        number_of_images: int = params['number-of-images']
        dump: bool = params['dump-images']

        # Databases
        db_prop_l = []
        for elem in params['database']:
            identi: str = elem['identifier']
            database_path: str = elem['path']
            iteration_period: Optional[int] = elem.get('iteration-period')
            time_period: Optional[float] = elem.get('time-period')
            db_props = DatabaseProps(identi, database_path, time_period,
                                     iteration_period)
            db_prop_l.append(db_props)
        return number_of_images, dump, db_prop_l

    def __fill_in_bc_props(self) -> Tuple[BoundaryType, BoundaryType]:
        """
        Returns the quantities needed to fill boundary conditions properties
            - left
            - right
        """
        params = self.__datadoc['boundary-conditions']
        left_boundary = self.__get_boundary_condition_def(params['left-boundary'])
        right_boundary = self.__get_boundary_condition_def(params['right-boundary'])
        return left_boundary, right_boundary

    @staticmethod
    def __get_boundary_condition_def(info) -> BoundaryType:
        """
        Returns the BoundaryType object corresponding to the info datas.
        """
        type_bc: str = info['type'].lower()
        func_name: str = info['bc-law'].lower()
        if func_name == 'constant':
            cvf = ConstantValueFunctionProps(info['value'])
            return BoundaryType(type_bc, cvf)

        if func_name == "twostep":
            tsf = TwoStepsFunctionProps(info['value1'], info['value2'], info['time-activation'])
            return BoundaryType(type_bc, tsf)

        if func_name == "ramp":
            raf = RampFunctionProps(info['value1'], info['value2'],
                                    info['time-activation-value-1'],
                                    info['time-activation-value-2'])
            return BoundaryType(type_bc, raf)

        if func_name == "marchtable":
            mtf = MarchTableFunctionProps(info['value'])
            return BoundaryType(type_bc, mtf)

        if func_name == "creneauramp":
            f_ramp = RampFunctionProps(
                info['initial-value'], info['plateau-value'],
                info['start-first-ramp-time'], info['reach-value2-time'])
            s_ramp = RampFunctionProps(
                info['plateau-value'], info['end-value'],
                info['start-second-ramp-time'], info['reach-value3-time'])
            srf = SuccessiveRampFunctionProps(f_ramp, s_ramp)
            return BoundaryType(type_bc, srf)

        raise ValueError(f"Unkown function type {func_name}."
                         "Please use one of [constant, twostep, ramp, marchtable, creneauramp]")

    @staticmethod
    def __get_damage_props(matter) -> Optional[Tuple[CohesiveZoneModelProps, str]]:
        """
        Returns the values needed to fill the damage model properties:
            - the cohesive model
            - the cohesive model name
        """
        try:
            params = matter['failure']['damage-treatment']['cohesive-model']
        except KeyError:
            return None

        cohesive_model_name = params['name'].lower()

        cohesive_strength = params['coefficients']['cohesive-strength']
        critical_separation = params['coefficients']['critical-separation']

        unloading_model_name = params['unloading-model']['name'].lower()
        unloading_model_slope = params['unloading-model']['slope']

        if unloading_model_name == "progressiveunloading":
            unloading_model_props: UnloadingModelProps = (
                ConstantStiffnessUnloadingProps(unloading_model_slope))
        elif unloading_model_name == "lossofstiffnessunloading":
            unloading_model_props = LossOfStiffnessUnloadingProps()
        else:
            raise ValueError(f"Unknwown unloading model name: {unloading_model_name} "
                             "Please choose among (progressiveunloading, "
                             " lossofstiffnessunloading)")

        if cohesive_model_name == "linear":
            cohesive_model_props: CohesiveZoneModelProps = LinearCohesiveZoneModelProps(
                cohesive_strength, critical_separation, unloading_model_props)
        elif cohesive_model_name == "bilinear":
            cohesive_model_props = BilinearCohesiveZoneModelProps(
                cohesive_strength, critical_separation, unloading_model_props,
                params['coefficients']['separation-at-point-1'],
                params['coefficients']['stress-at-point-1'])
        elif cohesive_model_name == "trilinear":
            cohesive_model_props = TrilinearCohesiveZoneModelProps(
                cohesive_strength, critical_separation, unloading_model_props,
                params['coefficients']['separation-at-point-1'],
                params['coefficients']['stress-at-point-1'],
                params['coefficients']['separation-at-point-2'],
                params['coefficients']['stress-at-point-2'])
        else:
            raise ValueError(f"Unknwon cohesive model: {cohesive_model_name} ."
                             "Please choose among (linear, bilinear, trilinear)")

        return cohesive_model_props, cohesive_model_name

    @staticmethod
    def __get_contact_props(matter) -> Optional[Tuple[ContactProps, str]]:
        """
        Returns the values needed to fill the contact model properties:
            - the contact model
            - the cohesive model name
        """
        try:
            params = matter['failure']['contact-treatment']['contact-model']
        except KeyError:
            return None
        contact_model_name = params['name'].lower()
        if contact_model_name == "penalty":
            penalty_stiffness: float = params['penalty-stiffness']
            contact_model_props: ContactProps = PenaltyContactProps(penalty_stiffness)
        else:
            raise ValueError(f"Unknwon contact model: {contact_model_name} ."
                             "Please choose among (penalty)")
        return contact_model_props, contact_model_name

    def __fill_in_material_props(self, material):
        """
        Returns the values needed to fill the material properties:
            - the initial values
            - the damage properties

        """
        init = InitialValues(*self.__get_initial_values(material))

        behavior = ConstitutiveModelProps(
            self.__get_equation_of_state_props(material),
            *self.__get_rheology_props(material))

        # TODO : rassembler tout ça !
        failure_treatment, failure_treatment_value, type_of_enrichment, lump_mass_matrix = (
            self.__get_failure_props(material))
        failure_criterion, failure_criterion_value = self.__get_failure_criterion_props(material)

        failure = FailureModelProps(failure_treatment, failure_treatment_value, type_of_enrichment,
                                    lump_mass_matrix, failure_criterion, failure_criterion_value)

        if failure.failure_treatment is not None and failure_criterion is None:
            raise ValueError("Failure criterion expected. "
                             "No failure criterion is specified or "
                             "specified criterion not understood")

        dmg_props = self.__get_damage_props(material)
        if dmg_props:
            damage = DamageModelProps(*dmg_props)
        else:
            damage = None

        if failure.failure_treatment is not None and damage and damage.cohesive_model is not None:
            if (failure_criterion_value != damage.cohesive_model.cohesive_strength and
                    isinstance(failure_criterion, MaximalStressCriterionProps)):
                print("Failure criterion value and cohesive strength have different value. "
                      "This may result in errors in the future")

        contact_props = self.__get_contact_props(material)
        if contact_props:
            contact = ContactModelProps(*contact_props)
        else:
            contact = None

        return init, behavior, failure, damage, contact

    def __get_initial_values(self, matter) -> Tuple[float, float, float, float,
                                                    float, float, float]:
        """
        Reads the XDATA file and find the initialization quantities for the material matter

        :param matter: material to be considered
        :return: initial thermodynamical quantities
        """
        params = matter['initialization']
        velocity = params['initial-velocity']

        json_path = params['init-thermo']
        with open(self._datafile_dir / json_path, 'r') as json_fid:
            coef = json.load(json_fid)
            coef = coef["InitThermo"]
            density = float(coef["initial_density"])
            pressure = float(coef["initial_pressure"])
            temperature = float(coef["initial_temperature"])
            internal_energy = float(coef["initial_internal_energy"])

        yield_stress, shear_modulus = self.__get_yield_stress_and_shear_modulus(matter)
        return (velocity, pressure, temperature, density, internal_energy,
                yield_stress, shear_modulus)

    def __get_yield_stress_and_shear_modulus(self, matter) -> Tuple[float, float]:
        """
        Returns the yield stress and shear modulus read from the json file specified
        under the coefficients key
        """
        params = matter.get('rheology')
        if not params:
            return 0., 0.

        coefficients_file = params.get('coefficients')
        if not coefficients_file:
            return 0., 0.

        with open(self._datafile_dir / coefficients_file, 'r') as json_fid:
            coef = json.load(json_fid)
            yield_stress = float(coef['EPP']["yield_stress"])
            shear_modulus = float(coef['EPP']["shear_modulus"])
        return yield_stress, shear_modulus

    def __get_equation_of_state_props(self, matter) -> EquationOfStateProps:
        """
        Returns the properties of the equation of state read from the datafile

        :param matter: material to be considered
        :return: the equation of state parameters
        """
        params = matter['equation-of-state']
        if params['name'] == 'Mie-Gruneisen':
            json_path = params['coefficients']
            with open(self._datafile_dir / json_path, 'r') as json_fid:
                coef = json.load(json_fid)
                coef = coef["MieGruneisen"]
                # Lecture des paramètres
                params_key = ("ref_sound_velocity", "s1", "s2", "s3", "ref_density",
                              "coefficient_gruneisen", "param_b", "ref_internal_energy")
                params = [float(coef[p]) for p in params_key]
            # Création de l'équation d'état
            return MieGruneisenProps(*params)

        raise NotImplementedError("Only MieGruneisen equation of state is implemented for now")

    def __get_rheology_props(self, matter) -> Tuple[Optional[ShearModulusProps],
                                                    Optional[YieldStressProps],
                                                    Optional[PlasticityCriterionProps]]:
        """
        Reads the elasticity parameters for material matter

        :param matter: material to be considered
        :return: elasticity_model : elasticity model
                 plasticity_model : plasticity model
                 plasticity_criterion : plasticity criterion
        """
        params = matter.get('rheology')
        if not params:
            return None, None, None

        yield_stress, shear_modulus = self.__get_yield_stress_and_shear_modulus(matter)
        elasticity_model_name = params['elasticity-model']
        if elasticity_model_name != "Linear":
            raise ValueError(f"Model {elasticity_model_name} not implemented. Choose Linear")
        elasticity_model = ConstantShearModulusProps(shear_modulus)

        plasticity_model_name = params.get('plasticity-model')
        if not plasticity_model_name:
            return elasticity_model, None, None

        if plasticity_model_name != "EPP":
            raise ValueError(f"Model {plasticity_model_name} not implemented. Choose EPP")
        plasticity_model = ConstantYieldStressProps(yield_stress)

        plasticity_criterion_name = params['plasticity-criterion']
        if plasticity_criterion_name == "VonMises":
            plasticity_criterion = VonMisesCriterionProps()
        else:
            raise ValueError(f"Plasticity criterion {plasticity_criterion_name}. Choose VonMises")
        return elasticity_model, plasticity_model, plasticity_criterion

    @staticmethod
    def __get_failure_props(matter) -> Tuple[Optional[str], Optional[float],
                                             Optional[str], Optional[str]]:
        """
        Returns the data needed to fill the FailureModel props

            - failure_model : rupture treatment model name
            - failure_treatment_value : position of disconitnuity in cracked element
            -                           or imposed pressure
            - type_of_enrichment : Hansbo
            - lump_mass_matrix : lumping strategy
        """
        failure_data = matter.get('failure')
        if not failure_data:
            return None, None, None, None

        failure_treatment_data = matter['failure']['failure-treatment']

        failure_treatment = failure_treatment_data.get('name')
        if failure_treatment is not None:
            if failure_treatment not in ["ImposedPressure", "Enrichment"]:
                raise ValueError(f"Unknown failure treatment {failure_treatment}."
                                 "Please choose among (ImposedPressure, Enrichment)")
            failure_treatment_value = failure_treatment_data['value']

        if failure_treatment == "Enrichment":
            type_of_enrichment = failure_treatment_data['type-of-enrichment']
            if type_of_enrichment not in ['Hansbo'] and failure_treatment == "Enrichment":
                raise ValueError(f"Unknown enrichment type {type_of_enrichment}. "
                                 "Only Hansbo is available for now")
            lump_mass_matrix = failure_treatment_data['lump-mass-matrix']
            if lump_mass_matrix.lower() not in ["menouillard", "somme", "none"]:
                print(f"Unknown lumping technique {lump_mass_matrix}. "
                      "Only those lumps are possible : menouillard, somme \n"
                      "No lumping technique is applied. Mass matrix is consistent")
        else:
            failure_treatment_value = 0.
            type_of_enrichment = None
            lump_mass_matrix = None

        return failure_treatment, failure_treatment_value, type_of_enrichment, lump_mass_matrix

    @staticmethod
    def __get_failure_criterion_props(matter) -> Tuple[
            Optional[RuptureCriterionProps], float]:
        """
        Returns the failure criterion properties needed to fill the failure model properties

        :return: failure_criterion : the rupture criterion
                failure_criterion_value : the threshold value
        """
        failure_data = matter.get('failure')
        if not failure_data:
            return None, 0.

        failure_criterion_data = failure_data['failure-criterion']

        failure_criterion_name = failure_criterion_data['name']
        failure_criterion_value = failure_criterion_data['value']

        if failure_criterion_name == "MinimumPressure":
            failure_criterion: RuptureCriterionProps = (
                MinimumPressureCriterionProps(failure_criterion_value))
        elif failure_criterion_name == "Damage":
            failure_criterion = DamageCriterionProps(failure_criterion_value)
        elif failure_criterion_name == "HalfRodComparison":
            failure_criterion = HalfRodComparisonCriterionProps(failure_criterion_value)
        elif failure_criterion_name == "MaximalStress":
            failure_criterion = MaximalStressCriterionProps(failure_criterion_value)
        else:
            raise ValueError(f"Unknown failure criterion {failure_criterion_name}. "
                             "Please choose among (MinimumPressure, Damage, "
                             "HalfRodComparison, MaximalStress")

        return failure_criterion, failure_criterion_value

    def hasExternalSolver(self):  # pylint: disable=invalid-name
        """
        Returns True if an external solver is required in the data
        """
        return 'external-solver-library' in self.__datadoc['numeric-parameters'].keys()

    def getExternalSolverPath(self):  # pylint: disable=invalid-name
        """
        Returns the external solver required in the data
        """
        return self.__datadoc['numeric-parameters']['external-solver-library']


if __name__ == "__main__":
    # pylint: disable=invalid-name
    data = DataContainer(datafile_path='XDATA.json')
    print(data.numeric)
    print(data.geometric)
    print(data.time)
    print(data.output)
    print(data.boundary_condition)
    print(data.material_projectile)
    print(data.material_target)

    left_bc = data.boundary_condition.left_BC.law.build_custom_func()
    print(left_bc.evaluate(0))
    print(left_bc.evaluate(1.5e-6))
    print(left_bc.evaluate(10))
    right_bc = data.boundary_condition.right_BC.law.build_custom_func()
    print(right_bc.evaluate(0))
    print(right_bc.evaluate(10))

    print(data.material_target)
    print(data.material_projectile)
