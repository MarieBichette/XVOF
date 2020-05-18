.. XVOF documentation master file, created by
   sphinx-quickstart on Thu May 14 12:55:03 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

XVOF Documentation
================================

Cell
--------------------------------
.. module:: cell
.. autoclass:: Cell
   :members:
.. autoclass:: OneDimensionCell
   :members:
.. autoclass:: OneDimensionHansboEnrichedCell
   :members:

Cohesive Zone Model
--------------------------------
.. module:: cohesive_model
.. autoclass:: CohesiveLaw
   :members:
.. autoclass:: CohesiveZoneModel
   :members:

Cohesive Zone Model : unloading options
-------------------------------------------
.. module:: cohesive_model_unloading
.. autoclass:: UnloadingModelBase
   :members:
.. autoclass:: ConstantStiffnessUnloading
   :members:
.. autoclass:: LossOfStiffnessUnloading
   :members:

Contact
--------------------------------
.. module:: contact
.. autoclass:: ContactBase
   :members:
.. autoclass:: LagrangianMultiplierContact
   :members:
.. autoclass:: PenaltyContact
   :members:

Customized Boundary Condition Functions
-------------------------------------------
.. module:: custom_functions
.. autoclass:: CustomFunction
   :members:
.. autoclass:: ConstantValue
   :members:
.. autoclass:: MarchTable
   :members:
.. autoclass:: Ramp
   :members:
.. autoclass:: SuccessiveRamp
   :members:
.. autoclass:: TwoSteps
   :members:

Data
--------------------------------
.. module:: data
.. autoclass:: DataContainer
   :members:


Discontinuity
--------------------------------
.. module:: discontinuity
.. autoclass:: Discontinuity
   :members:

Equation Of State
--------------------------------
.. module:: equationsofstate
.. autoclass:: EquationOfStateBase
   :members:
.. autoclass:: MieGruneisen
   :members:

Fields
--------------------------------
.. module:: fields
.. autoclass:: Field
   :members:
.. autoclass:: FieldManager
   :members:

FigureManager (animation)
--------------------------------
.. module:: figure_manager
.. autoclass:: PhysicFigure
   :members:
.. autoclass:: FigureManager
   :members:

Mass Matrix
--------------------------------
.. module:: mass_matrix
.. autoclass:: EnrichedMassMatrix
   :members:
.. autoclass:: EnrichedMassMatrixConsistent
   :members:
.. autoclass:: EnrichedMassMatrixLump
   :members:
.. autoclass:: EnrichedMassMatrixLumpMenouillard
   :members:
.. autoclass:: EnrichedMassMatrixLumpSum
   :members:
.. autoclass:: OneDimensionMassMatrix
   :members:
.. autoclass:: SymNDArray
   :members:
.. autofunction:: compute_wilkins_mass_matrix
.. autofunction:: multiplication_masse
.. autofunction:: inverse_masse
.. autofunction:: lump_matrix

Mesh
--------------------------------
.. module:: mesh
.. autoclass:: Mesh1dEnriched
   :members:
.. autoclass:: Topology
   :members:
.. autoclass:: Topology1D
   :members:

Nodes
--------------------------------
.. module:: node
.. autoclass:: Node
   :members:
.. autoclass:: OneDimensionNode
   :members:
.. autoclass:: OneDimensionHansboEnrichedNode
   :members:

Outputs
--------------------------------
.. module:: output_manager
.. autoclass:: OutputDatabase
   :members:
.. autoclass:: OutputManager
   :members:
.. autoclass:: OutputTimeControler
   :members:
.. autoclass:: OutputDatabaseExploit
   :members:
.. autofunction:: build_node_true_field
.. autofunction:: build_cell_true_field

Plasticity criterion
--------------------------------
.. module:: plasticitycriterion
.. autoclass:: PlasticityCriterion
   :members:
.. autoclass:: VonMisesCriterion
   :members:

Rheology
--------------------------------
.. module:: rheology
.. autoclass:: ShearModulus
   :members:
.. autoclass:: ConstantShearModulus
   :members:
.. autoclass:: YieldStress
   :members:
.. autoclass:: ConstantYieldStress
   :members:

Rupture criterion
--------------------------------
.. module:: rupturecriterion
.. autoclass:: RuptureCriterion
   :members:
.. autoclass:: HalfRodComparisonCriterion
   :members:
.. autoclass:: DamageCriterion
   :members:
.. autoclass:: MinimumPressureCriterion
   :members:
.. autoclass:: MaximalStressCriterion
   :members:

Rupture treatments
--------------------------------
.. module:: rupturetreatment
.. autoclass:: RuptureTreatment
   :members:
.. autoclass:: ImposedPressure
   :members:
.. autoclass:: EnrichElement
   :members:

Solver
--------------------------------
.. module:: solver
.. autoclass:: NewtonRaphsonBase
   :members:
.. autoclass:: NewtonRaphson
   :members:
.. autoclass:: ClassicalNewtonRaphsonIncrement
   :members:

Utilities
--------------------------------
.. module:: utilities
.. autofunction:: timeit_file
.. autofunction:: compute_second_invariant
.. autofunction:: compute_trace
.. autofunction:: captured_output
.. autoclass:: Singleton
   :members:

Index and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
