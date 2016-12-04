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
#  =================================================

if __name__ == '__main__':
    np.set_printoptions(formatter={'float': '{: 25.23g}'.format})
    #
    simulation_time = 0.
    step = 0
    dt = InitialTimeStep
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
    my_fig_manager = FigureManager(my_mesh, dump=data.output.images_dump, show=data.output.images_show)
    if ImagesNumber != 0:
        delta_t_images = FinalTime / ImagesNumber
        my_fig_manager.populate_figs()
    else:
        delta_t_images = FinalTime * 2.0
    t_next_image = delta_t_images
    # ---------------------------------------------#
    #         NODAL MASS COMPUTATION               #
    # ---------------------------------------------#
    my_mesh.compute_cells_sizes()
    my_mesh.compute_nodes_masses()
    print "CALCULUS LAUNCHED!"
    compute_time = 0.
    while simulation_time < FinalTime:
        loop_begin_time = time.time()
        if step % 1000 == 0:
            msg = ("""Iteration {:<4d} -- Time : {:15.9g} seconds with"""
                   """ a time step of {:15.9g} seconds\n""").format(step, simulation_time, dt)
            print msg
        # ---------------------------------------------#
        #         NODES VELOCITIES COMPUTATION         #
        # ---------------------------------------------#
        my_mesh.compute_new_nodes_velocities(dt)
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
        dt = dt_crit
        step += 1
        loop_end_time = time.time()
        compute_time += loop_end_time - loop_begin_time
        # ---------------------------------------------#
        #                OUTPUT MANAGEMENT             #
        # ---------------------------------------------#
        if simulation_time > t_next_image:
            my_fig_manager.update_figs("t={:5.4g} us".format(simulation_time / 1.e-06))
            t_next_image += delta_t_images
    print "Total time spent in compute operation is : {:15.9g} seconds".format(compute_time)
    plt.show()
