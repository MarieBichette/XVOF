# -*- coding: utf-8 -*-
"""
Implementing the DataContainer class
"""
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union, Any

from xfv.src.utilities.singleton import Singleton
from xfv.src.data.type_checked_dataclass import TypeCheckedDataClass
from xfv.src.data.user_defined_functions_props import (UserDefinedFunctionPropsType,
                                                       ConstantValueFunctionProps,
                                                       TwoStepsFunctionProps,
                                                       RampFunctionProps,
                                                       MarchTableFunctionProps,
                                                       SuccessiveRampFunctionProps)
from xfv.src.data.cohesive_model_props import (CohesiveZoneModelProps,
                                               LinearCohesiveZoneModelProps,
                                               BilinearCohesiveZoneModelProps,
                                               TrilinearCohesiveZoneModelProps)
from xfv.src.data.unloading_model_props import (UnloadingModelProps,
                                                ConstantStiffnessUnloadingProps,
                                                LossOfStiffnessUnloadingProps)
from xfv.src.data.contact_props import (ContactProps, PenaltyContactProps,
                                        LagrangianMultiplierProps)
from xfv.src.data.equation_of_state_props import (EquationOfStateProps, MieGruneisenProps)
from xfv.src.data.yield_stress_props import (YieldStressProps, ConstantYieldStressProps)
from xfv.src.data.shear_modulus_props import (ShearModulusProps, ConstantShearModulusProps)
from xfv.src.data.plasticity_criterion_props import (PlasticityCriterionProps,
                                                     VonMisesCriterionProps)
from xfv.src.data.rupture_criterion_props import (RuptureCriterionProps,
                                                  DamageCriterionProps,
                                                  PorosityCriterionProps,
                                                  HalfRodComparisonCriterionProps,
                                                  MaximalStressCriterionProps,
                                                  MinimumPressureCriterionProps)
from xfv.src.data.porosity_model_props import (PorosityModelProps,
                                               JohnsonModelProps)
from xfv.src.data.enriched_mass_matrix_props import (EnrichedMassMatrixProps,
                                                     ConsistentMassMatrixProps,
                                                     LumpMenouillardMassMatrixProps,
                                                     LumpSumMassMatrixProps)


@dataclass  # pylint: disable=missing-class-docstring
class NumericalProps(TypeCheckedDataClass):
    a_pseudo: float
    b_pseudo: float
    cfl: float
    cfl_pseudo: float
    consistent_mass_matrix_on_last_cells: bool

    def __post_init__(self):
        super().__post_init__()  # type checking first
        self._ensure_positivity('a_pseudo', 'b_pseudo', 'cfl', 'cfl_pseudo')


@dataclass  # pylint: disable=missing-class-docstring
class GeometricalProps(TypeCheckedDataClass):
    section: float
    initial_interface_position: float

    def __post_init__(self):
        super().__post_init__()  # type checking first
        self._ensure_strict_positivity('section')
        self._ensure_positivity('initial_interface_position')


@dataclass  # pylint: disable=missing-class-docstring
class TimeProps(TypeCheckedDataClass):
    initial_time_step: float
    final_time: float
    is_time_step_constant: bool
    time_step_reduction_factor_for_failure: Optional[float]

    def __post_init__(self):
        super().__post_init__()
        self._ensure_positivity('initial_time_step', 'time_step_reduction_factor_for_failure')
        self._ensure_strict_positivity('final_time')


@dataclass  # pylint: disable=missing-class-docstring
class DatabaseProps(TypeCheckedDataClass):
    identifier: str
    path: str
    time_period: Optional[float]
    iteration_period: Optional[int]

    def __post_init__(self):
        super().__post_init__()
        self._ensure_strict_positivity('time_period', 'iteration_period')
        if self.time_period is not None and self.iteration_period is not None:
            raise ValueError("Please provide one of (time-period, iteration-period) "
                             "but not both!")


ALL_VARIABLES = ["NodeVelocity", "NodeCoordinates", "CellSize", "Pressure", "Density",
                 "InternalEnergy", "SoundVelocity", "ArtificialViscosity", "Stress",
                 "DeviatoricStress", "EquivalentPlasticStrainRate", "PlasticStrainRate",
                 "Porosity", "CohesiveForce", "DiscontinuityOpening", "ShearModulus", "YieldStress"]


@dataclass  # pylint: disable=missing-class-docstring
class OutputProps(TypeCheckedDataClass):
    number_of_images: int
    dump: bool
    databases: List[DatabaseProps]
    variables: List[str]

    def __post_init__(self):
        super().__post_init__()
        self._ensure_positivity('number_of_images')
        self._ensure_list_value_in("variables", ALL_VARIABLES)


@dataclass  # pylint: disable=missing-class-docstring
class BoundaryType(TypeCheckedDataClass):
    type_bc: str
    law: UserDefinedFunctionPropsType

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        self._ensure_value_in('type_bc', ('velocity', 'pressure'))


@dataclass  # pylint: disable=missing-class-docstring
class BoundaryConditionsProps(TypeCheckedDataClass):
    left_BC: BoundaryType  # pylint: disable=invalid-name
    right_BC: BoundaryType  # pylint: disable=invalid-name


@dataclass  # pylint: disable=missing-class-docstring
class PorosityModProps(TypeCheckedDataClass):
    porosity_model: PorosityModelProps
    name: str


@dataclass  # pylint: disable=missing-class-docstring
class ContactModelProps(TypeCheckedDataClass):
    contact_model: ContactProps
    name: str
    # Do not check if name is among authorized values because
    # it is done in one of the DataContainer's method and
    # moving the test here, implies to allow the contact_model to be None


@dataclass  # pylint: disable=missing-class-docstring
class InitialValues(TypeCheckedDataClass):
    velocity_init: float
    pression_init: float
    temp_init: float
    rho_init: float
    energie_init: float
    yield_stress_init: float
    shear_modulus_init: float
    porosity_init: float

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        self._ensure_strict_positivity('rho_init', 'temp_init', 'porosity_init')
        self._ensure_positivity('yield_stress_init', 'shear_modulus_init')

        if self.porosity_init < 1.0:
            raise ValueError(f"{self.porosity_init} < 1.0\n"
                             "is not a possible value for initial porosity.")
        # TODO: one the initialization via eos will be done, check that only two fields
        # are not null (v and e or v and T)


@dataclass  # pylint: disable=missing-class-docstring
class ConstitutiveModelProps(TypeCheckedDataClass):
    eos: EquationOfStateProps
    elasticity_model: Optional[ShearModulusProps]
    plasticity_model: Optional[YieldStressProps]
    plasticity_criterion: Optional[PlasticityCriterionProps]


@dataclass  # pylint: disable=missing-class-docstring
class FailureModelProps(TypeCheckedDataClass):
    failure_treatment: Optional[str]
    failure_treatment_value: Optional[float]
    lump_mass_matrix: Optional[EnrichedMassMatrixProps]
    failure_criterion: Optional[RuptureCriterionProps]
    failure_criterion_value: Optional[float]
    failure_criterion_index: Optional[int]

    def __post_init__(self):
        super().__post_init__()  # typecheck first
        if self.failure_criterion is None and self.failure_treatment is not None:
            raise ValueError("A failure criterion is required if failure treatment is set")


@dataclass  # pylint: disable=missing-class-docstring
class MaterialProps(TypeCheckedDataClass):
    initial_values: InitialValues
    constitutive_model: ConstitutiveModelProps
    failure_model: FailureModelProps
    cohesive_model: Optional[CohesiveZoneModelProps]
    porosity_model: Optional[PorosityModProps]
    contact_model: Optional[ContactModelProps]


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

    def __fill_in_numerical_props(self) -> Tuple[float, float, float, float, bool]:
        """
        Returns the quantities needed to fill numerical properties:
            - coefficient of linear artificial viscosity
            - coefficient of quadratic artificial viscosity
            - CFL coefficient
            - CFL coefficient of artificial viscosity
        """
        params: Dict[str, float] = self.__datadoc['numeric-parameters']
        last_cells_consistent_mass_matrix: bool = params.get(
            'consistent-mass-matrix-on-last-cells', False)
        return (params['quadratic-pseudo'], params['linear-pseudo'],
                params['cfl'], params['cfl-pseudo'], last_cells_consistent_mass_matrix)

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

    def __fill_in_output_props(self) -> Tuple[int, bool, List[DatabaseProps], List[str]]:
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

        variables_l = []
        if params['variables'][0] == "All":
            variables_l = ALL_VARIABLES
        else:
            for var in params['variables']:
                variables_l.append(var)

        return number_of_images, dump, db_prop_l, variables_l

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
    def __get_cohesive_model_props(matter) -> Optional[CohesiveZoneModelProps]:
        """
        Returns the values needed to fill the cohesive zone model properties:
            - the cohesive model
            - the cohesive model name
        """
        try:
            params = matter['failure']['cohesive-model']
        except KeyError:
            return None

        cohesive_model_name = params['name'].lower()

        cohesive_strength = params['coefficients']['cohesive-strength']
        critical_separation = params['coefficients']['critical-separation']

        unloading_model_name = params['unloading-model']['name'].lower()
        unloading_model_slope: Optional[float] = params['unloading-model'].get('slope')

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

        return cohesive_model_props

    @staticmethod
    def __get_porosity_model_props(matter) -> Optional[Tuple[PorosityModelProps, str]]:
        """
        Returns the values needed to fill the damage model properties:
            - the porosity model
            - the porosity model name
        """
        try:
            params = matter['failure']['porosity-model']
        except KeyError:
            return None

        try:
            porosity_model_name = params['name'].lower()
            if porosity_model_name == "johnsonmodel":
                initial_porosity_for_johnson = params['coefficients']['initial-porosity']
                effective_strength_for_johnson = params['coefficients']['effective-strength']
                viscosity_for_johnson = params['coefficients']['viscosity']
                porosity_model_props: PorosityModelProps = JohnsonModelProps(
                    initial_porosity_for_johnson,
                    effective_strength_for_johnson,
                    viscosity_for_johnson)
        except:
            raise ValueError(f"No keyword 'name' for porosity model name: {porosity_model_name}."
                             "Please choose among (JohnsonModel)."
                             "For JohnsonModel, Keyword 'coefficients' must be defined. ")

        return porosity_model_props, porosity_model_name

    @staticmethod
    def __get_contact_props(matter) -> Optional[Tuple[ContactProps, str]]:
        """
        Returns the values needed to fill the contact model properties:
            - the contact model
            - the cohesive model name
        """
        try:
            params = matter['failure']['contact-treatment']
        except KeyError:
            return None
        contact_model_name = params['name'].lower()
        if contact_model_name == "penalty":
            penalty_stiffness: float = params['penalty-stiffness']
            contact_model_props: ContactProps = PenaltyContactProps(penalty_stiffness)
        elif contact_model_name == "lagrangianmultiplier":
            contact_model_props: ContactProps = LagrangianMultiplierProps()
        else:
            raise ValueError(f"Unknwon contact model: {contact_model_name} ."
                             "Please choose among (Penalty|LagrangianMultiplier)")
        return contact_model_props, contact_model_name

    def __fill_in_material_props(self, material) -> Tuple[InitialValues,
                                                          ConstitutiveModelProps,
                                                          FailureModelProps,
                                                          Optional[CohesiveZoneModelProps],
                                                          Optional[PorosityModProps],
                                                          Optional[ContactModelProps]]:
        """
        Returns the values needed to fill the material properties:
            - the initial values
            - the damage properties

        """
        # Initialisation
        init = InitialValues(*self.__get_initial_values(material))

        # Bulk material behavior
        behavior = ConstitutiveModelProps(
            self.__get_equation_of_state_props(material),
            *self.__get_rheology_props(material))

        # Failure treatment
        failure_treatment, failure_treatment_value, lump_mass_matrix = (
            self.__get_failure_props(material))
        failure_criterion, failure_criterion_value, failure_index = (
            self.__get_failure_criterion_props(material))

        failure = FailureModelProps(failure_treatment, failure_treatment_value,
                                    lump_mass_matrix, failure_criterion, failure_criterion_value,
                                    failure_index)

        # Surface degradation behavior
        cohesive_model: Optional[CohesiveZoneModelProps] = self.__get_cohesive_model_props(material)
        if failure.failure_treatment is not None and cohesive_model is not None:
            if (failure_criterion_value != cohesive_model.cohesive_strength and
                    isinstance(failure_criterion, MaximalStressCriterionProps)):
                print("Failure criterion value and cohesive strength have different value. "
                      "This may result in errors in the future")

        # Porosity model
        porosity_model_props = self.__get_porosity_model_props(material)
        porosity_model: Optional[PorosityModProps] = None
        if porosity_model_props:
            porosity_model = PorosityModProps(*porosity_model_props)

        # Contact between discontinuity boundaries treatment
        contact_props = self.__get_contact_props(material)
        if contact_props:
            contact: Optional[ContactModelProps] = ContactModelProps(*contact_props)
        else:
            contact = None

        return init, behavior, failure, cohesive_model, porosity_model, contact

    def __get_initial_values(self, matter) -> Tuple[float, float, float, float,
                                                    float, float, float, float]:
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

        porosity = self.__get_initial_porosity(matter)

        return (velocity, pressure, temperature, density, internal_energy,
                yield_stress, shear_modulus, porosity)

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

    def __get_initial_porosity(self, matter)->float:
        """
        Returns the initial porosity. We choose to initiate porosity to 1.0 (no porosity). 
        The value defined in the porosity model is assigned to the porosity variable in the porosity_model.
        """
        return 1.0

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
                # Read parameters
                params_key = ("ref_sound_velocity", "s1", "s2", "s3", "ref_density",
                              "coefficient_gruneisen", "param_b", "ref_internal_energy")
                params = [float(coef[p]) for p in params_key]
            # Returns the eos properties
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
    def __get_failure_props(matter: Any) -> Tuple[Optional[str], Optional[float],
                                                  Optional[EnrichedMassMatrixProps]]:
        """
        Returns the data needed to fill the FailureModel props

            - failure_model : rupture treatment model name
            - failure_treatment_value : position of discontinuity in cracked element
                                        or imposed pressure
            - lump_mass_matrix : lumping strategy
        """
        failure_data = matter.get('failure')
        if not failure_data:
            return None, None, None

        failure_treatment_data = matter['failure']['failure-treatment']
        failure_treatment = failure_treatment_data.get('name')

        if failure_treatment is not None:
            # Failure treatment is either Enrichment or ImposedPressure
            failure_treatment_value = failure_treatment_data['value']

            # Case enrichment :
            if failure_treatment == "Enrichment":
                # Choice of the enriched mass matrix lumping
                lump_name: str = failure_treatment_data['lump-mass-matrix']
                if lump_name.lower() == "menouillard":
                    lump_mass_matrix = LumpMenouillardMassMatrixProps()
                elif lump_name.lower() == "somme":
                    lump_mass_matrix = LumpSumMassMatrixProps()
                else:
                    print("No lump (menouillard|somme). Mass matrix is consistent")
                    lump_mass_matrix = ConsistentMassMatrixProps()

            # Case ImposedPressure :
            elif failure_treatment == "ImposedPressure":
                lump_mass_matrix = None

            else:
                raise ValueError(f"Unknown failure treatment {failure_treatment}."
                                 "Please choose among (ImposedPressure, Enrichment)")
        else:
            failure_treatment_value = 0.
            lump_mass_matrix = None
        return failure_treatment, failure_treatment_value, lump_mass_matrix

    @staticmethod
    def __get_failure_criterion_props(matter) -> Tuple[
            Optional[RuptureCriterionProps], Optional[float], Optional[int]]:
        """
        Returns the failure criterion properties needed to fill the failure model properties

        :return: failure_criterion : the rupture criterion
                failure_criterion_value : the threshold value
                failure_cell_index: the index of the cracked cell at the beginning of the simulation
        """
        failure_data = matter.get('failure')
        if not failure_data:
            return None, 0., 0

        failure_criterion_data = failure_data['failure-criterion']
        fail_crit_name: str = failure_criterion_data['name']
        fail_crit_value: Optional[float] = failure_criterion_data.get('value')
        failure_cell_index: Optional[int] = failure_criterion_data.get('index')

        if fail_crit_name == "MinimumPressure":
            failure_criterion = MinimumPressureCriterionProps(fail_crit_value)
        elif fail_crit_name == "Damage":
            failure_criterion = DamageCriterionProps(fail_crit_value)
        elif fail_crit_name == "Porosity":
            failure_criterion = PorosityCriterionProps(fail_crit_value)
        elif fail_crit_name == "HalfRodComparison":
            failure_criterion = HalfRodComparisonCriterionProps(failure_cell_index)
        elif fail_crit_name == "MaximalStress":
            failure_criterion = MaximalStressCriterionProps(fail_crit_value)
        else:
            raise ValueError(f"Unknown failure criterion {fail_crit_name}. "
                             "Please choose among (MinimumPressure, Damage, "
                             "HalfRodComparison, MaximalStress")

        return failure_criterion, fail_crit_value, failure_cell_index


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
