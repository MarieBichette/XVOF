.. XVOF documentation master file, created by
   sphinx-quickstart on Thu May 14 12:55:03 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

XVOF
================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Cell
================================
.. module:: cell
.. autoclass:: Cell
   :members:
.. autoclass:: OneDimensionCell
   :members:
.. autoclass:: OneDimensionHansboEnrichedCell
   :members:

Cohesive Zone Model
================================
.. module:: cohesive_model
.. autoclass:: CohesiveLaw
   :members:
.. autoclass:: CohesiveZoneModel
   :members:

Cohesive Zone Model : unloading options
=============================================
.. module:: cohesive_model_unloading
.. autoclass:: UnloadingModelBase
   :members:
.. autoclass:: ConstantStiffnessUnloading
   :members:
.. autoclass:: LossOfStiffnessUnloading
   :members:

Contact
================================
.. module:: contact
.. autoclass:: ContactBase
   :members:
.. autoclass:: LagrangianMultiplierContact
   :members:
.. autoclass:: PenaltyContact
   :members:

Customized Boundary Condition Functions
===========================================
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
==================
.. module:: data
.. autoclass:: DataContainer
   :members:


Discontinuity
===================
.. module:: discontinuity
.. autoclass:: Discontinuity
   :members:

Equation Of State
====================
.. module:: equationsofstate
.. autoclass:: EquationOfStateBase
   :members:
.. autoclass:: MieGruneisen
   :members:

.. module:: node
.. autoclass:: Node
   :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
