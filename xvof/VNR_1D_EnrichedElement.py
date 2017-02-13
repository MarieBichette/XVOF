#!/usr/bin/env python2.7
# -*- coding: iso-8859-15 -*-
"""
:todo: merge with VNR_1D_ClassicalElement
:todo: increase perf of external solver library
:todo: correct pb of variable time step
"""
import matplotlib

matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
from xvof.figure_manager.figure_manager import FigureManager
from xvof.pressurelaw.constantpressure import ConstantPressure
from xvof.pressurelaw.twostepspressure import TwoStepsPressure
from xvof.rupturecriterion.minimumpressure import MinimumPressureCriterion
from xvof.data.data_container import DataContainer
from xvof.mesh.mesh1denriched import Mesh1dEnriched
from xvof.output_manager.hdf5database import Hdf5Database
from xvof.output_manager.output_manager import OutputManager
import time

# SIMULATION PARAMETERS
data = DataContainer()
# TIME MANAGEMENT
FinalTime = data.time.final_time
InitialTimeStep = data.time.initial_time_step
# GEOMETRY
Length = data.geometric.length
# MATTER
InitialPressure = data.material.pression_init
#  NUMERIC
NumberOfElements = data.numeric.cells_number
# LOADING
LeftBoundaryPressure = TwoStepsPressure(15e+09, InitialPressure, 2.0e-06)
RightBoundaryPressure = ConstantPressure(InitialPressure)
RuptureCriterion = MinimumPressureCriterion(-7.0e+09)
RuptureTreatment = data.material.damage_treatment(data.material.damage_treatment_value)
# OUTPUT
ImagesNumber = data.output.number_of_images
TheOutputManager = OutputManager()
#  =================================================

if __name__ == '__main__':
    np.set_printoptions(formatter={'float': '{: 25.23g}'.format})
    #
    simulation_time = 0.
    step = 0
    dt = InitialTimeStep
    dt_staggered = dt
    dt_crit = 2 * dt
    # ---------------------------------------------#
    #         MESH CREATION                        #
    # ---------------------------------------------#
    coord_init = np.zeros([NumberOfElements + 1, 1], dtype=np.float64, order='C')
    coord_init[:, 0] = np.linspace(0, Length, NumberOfElements + 1)
    vit_init = np.zeros([NumberOfElements + 1, 1], dtype=np.float64, order='C')
    my_mesh = Mesh1dEnriched(initial_coordinates=coord_init, initial_velocities=vit_init)
    # ---------------------------------------------#
    #  FIGURES MANAGER SETUP                       #
    # ---------------------------------------------#
    TheFigureManager = FigureManager(my_mesh, dump=data.output.images_dump, show=data.output.images_show)
    if ImagesNumber != 0:
        TheFigureManager.set_time_controler(FinalTime / ImagesNumber)
        TheFigureManager.populate_figs()
    # ---------------------------------------------#
    #  OUTPUT MANAGER SETUP                        #
    # ---------------------------------------------#
    db1 = Hdf5Database("./simulation_output.hdf5")
    db2 = Hdf5Database("./simulation_output_2.hdf5")
    TheOutputManager.register_database_iteration_ctrl("GenericDatabase", db1, 10)
    TheOutputManager.register_database_time_ctrl("GenericDatabase2", db2, 1e-06)
    TheOutputManager.register_field("NodeCoordinates", my_mesh.nodes, "xt", "GenericDatabase", "GenericDatabase2")
    TheOutputManager.register_field("NodeStatus", my_mesh.nodes, "enriched", "GenericDatabase", "GenericDatabase2")
    TheOutputManager.register_field("CellStatus", my_mesh.cells, "enriched", "GenericDatabase")
    TheOutputManager.register_field("ClassicalNodeVelocity", my_mesh.nodes, "umundemi", "GenericDatabase", "GenericDatabase2")
    TheOutputManager.register_field("EnrichedNodeVelocity", my_mesh.nodes, "umundemi_enriched", "GenericDatabase", "GenericDatabase2")
    TheOutputManager.register_field("ClassicalPressure", my_mesh.cells.pressure, "new_value", "GenericDatabase")
    TheOutputManager.register_field("EnrichedPressure", my_mesh.cells.pressure, "new_enr_value", "GenericDatabase")
    TheOutputManager.register_field("ClassicalDensity", my_mesh.cells.density, "new_value", "GenericDatabase")
    TheOutputManager.register_field("EnrichedDensity", my_mesh.cells.density, "new_enr_value", "GenericDatabase")
    TheOutputManager.register_field("ClassicalInternalEnergy", my_mesh.cells.energy, "new_value", "GenericDatabase")
    TheOutputManager.register_field("EnrichedInternalEnergy", my_mesh.cells.energy, "new_enr_value", "GenericDatabase")
    TheOutputManager.register_field("ClassicalSoundVelocity", my_mesh.cells.sound_velocity, "new_value", "GenericDatabase")
    TheOutputManager.register_field("EnrichedSoundVelocity", my_mesh.cells.sound_velocity, "new_enr_value", "GenericDatabase")
    TheOutputManager.register_field("ClassicalArtificalViscosity", my_mesh.cells.pseudo, "new_value", "GenericDatabase")
    TheOutputManager.register_field("EnrichedArtificalViscosity", my_mesh.cells.pseudo, "new_enr_value", "GenericDatabase")
    # ---------------------------------------------#
    #         NODAL MASS COMPUTATION               #
    # ---------------------------------------------#
    my_mesh.compute_cells_sizes()
    my_mesh.compute_cells_masses()
    my_mesh.compute_nodes_masses()
    print "CALCULUS LAUNCHED!"
    compute_time = 0.
    while simulation_time < FinalTime:
        loop_begin_time = time.time()
        if step % 1000 == 0:
            msg = ("""Iteration {:<4d} -- Time : {:15.9g} seconds with"""
                   """ a time step of {:15.9g} seconds and a staggered time step of {:15.9g}\n"""). \
                format(step, simulation_time, dt, dt_staggered)
            print msg
        # ---------------------------------------------#
        #         NODES VELOCITIES COMPUTATION         #
        # ---------------------------------------------#
        my_mesh.compute_new_nodes_velocities(dt_staggered)
        # ---------------------------------------------#
        #         NODES COORDINATES COMPUTATION        #
        # ---------------------------------------------#
        my_mesh.compute_new_nodes_coordinates(dt)
        # ---------------------------------------------#
        #         CELLS VOLUMES COMPUTATION            #
        # ---------------------------------------------#
        my_mesh.compute_new_cells_sizes(dt)
        # ---------------------------------------------#
        #         CELLS DENSITIES COMPUTATION          #
        # ---------------------------------------------#
        my_mesh.compute_new_cells_densities()
        # ---------------------------------------------#
        #         CELLS PRESSURES COMPUTATION          #
        # ---------------------------------------------#
        my_mesh.compute_new_cells_pressures()
        # ---------------------------------------------#
        #              RUPTURE                         #
        # ---------------------------------------------#
        my_mesh.get_ruptured_cells(RuptureCriterion)
        my_mesh.apply_rupture_treatment(RuptureTreatment)
        # ---------------------------------------------#
        #         NODES FORCES COMPUTATION             #
        # ---------------------------------------------#
        my_mesh.compute_new_nodes_forces()
        # ---------------------------------------------#
        #         LOADING                              #
        # ---------------------------------------------#
        my_mesh.apply_pressure('left', LeftBoundaryPressure.evaluate(simulation_time))
        my_mesh.apply_pressure('right', RightBoundaryPressure.evaluate(simulation_time))
        # ---------------------------------------------#
        #         TIME STEP COMPUTATION                #
        # ---------------------------------------------#
        dt_crit = my_mesh.compute_new_time_step()
        # ---------------------------------------------#
        #         PSEUDOVISCOSITY COMPUTATION          #
        # ---------------------------------------------#
        my_mesh.compute_new_cells_pseudo_viscosity(dt)
        # ---------------------------------------------#
        #                INCREMENTATION                #
        # ---------------------------------------------#
        my_mesh.increment()
        simulation_time += dt
        if not data.time.is_time_step_constant:
            dt_staggered = 0.5 * (dt_crit + dt)
            dt = dt_crit
        step += 1
        loop_end_time = time.time()
        compute_time += loop_end_time - loop_begin_time
        # ---------------------------------------------#
        #                OUTPUT MANAGEMENT             #
        # ---------------------------------------------#
        TheOutputManager.update(simulation_time, step)
        TheFigureManager.update(simulation_time, step)
    print "Total time spent in compute operation is : {:15.9g} seconds".format(compute_time)
    #plt.show()
    TheOutputManager.finalize()
