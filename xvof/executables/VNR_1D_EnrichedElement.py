#!/usr/bin/env python2.7
"""
todo: to complete
"""
import os.path
import sys
import time
import matplotlib  # pylint: disable=unused-import
#matplotlib.use('Qt4Agg')


import matplotlib.pyplot as plt
import numpy as np

from xvof.src.figure_manager.figure_manager      import FigureManager
from xvof.src.data.data_container                import DataContainer
from xvof.src.mesh.mesh1denriched                import Mesh1dEnriched
from xvof.src.data.save_time_data                import CellTimeData, NodeTimeData
from xvof.src.output_manager.outputmanager       import OutputManager
from xvof.src.output_manager.outputdatabase      import OutputDatabase
from xvof.src.rupturetreatment.enrichelement     import EnrichElement
from xvof.src.rupturetreatment.imposedpressure   import ImposedPressure


def main():
    """
    Launch the program
    """
    # ------------------------------------------------------------------
    #             PARAMETERS INITIALIZATION
    # ------------------------------------------------------------------
    # ---- # DATA FILES
    #path = os.path.curdir
    path = sys.argv[1]  # prend en argument le chemin vers le cas
    data = DataContainer(os.path.join(path, "XDATA.xml"))
    meshfile = os.path.join(path, "mesh.txt")
    print "Running simulation for {:s}".format(os.path.normpath(os.path.abspath(path)))

    # ---- # TIME MANAGEMENT
    final_time = data.time.final_time
    initial_time_step = data.time.initial_time_step

    # ---- # LOADING
    left_boundary_condition = data.boundary_condition.left_BC
    right_boundary_condition = data.boundary_condition.right_BC

    # ---- # RUPTURE
    rupture_criterion = data.material_target.failure_model.failure_criterion
    rupture_treatment = None
    if data.material_target.failure_model.failure_treatment == "ImposedPressure":
        rupture_treatment = ImposedPressure(
            data.material_target.failure_model.failure_treatment_value)
    elif data.material_target.failure_model.failure_treatment == "Enrichment":
        rupture_treatment = EnrichElement(
            data.material_target.failure_model.failure_treatment_value)
        type_of_enrichment = data.material_target.failure_model.type_of_enrichment
        print "Enrichment method : {}".format(type_of_enrichment)

    # ---- # PLASTICITY
    target_plasticity_criterion = data.material_target.constitutive_model.plasticity_criterion
    projectile_plasticity_criterion = (
        data.material_projectile.constitutive_model.plasticity_criterion)

    # ---- # OUTPUT
    images_number = data.output.number_of_images
    the_output_mng = OutputManager()

    # Initialization for time figure plot
    cells_for_time_figure = []
    nodes_for_time_figure = []
    str_nodes_for_time_figure = data.output.nodes_numbers
    str_cells_for_time_figure = data.output.cells_numbers
    try:
        for cell_number in str_cells_for_time_figure:
            cells_for_time_figure.append(int(cell_number))
    except TypeError:
        pass

    try:
        for node_number in str_nodes_for_time_figure:
            nodes_for_time_figure.append(int(node_number))
    except TypeError:
        pass
    history_list = [CellTimeData(cell_id, path) for cell_id in cells_for_time_figure]
    history_list += [NodeTimeData(node_id, path) for node_id in nodes_for_time_figure]

    np.set_printoptions(formatter={'float': '{: 25.23g}'.format})
    #
    simulation_time = 0.
    step = 0
    dt = initial_time_step  # pylint: disable=invalid-name
    dt_staggered = dt / 2
    # Le premier increment pour la vitesse a un pas de temps de dt/2 pour
    # tenir compte du fait que les vitesses sont
    # init a t=0 et pas t = -1/2 dt (pour garderl'ordre 2 a l'init).
    # Le dt_staggered est ensuite remis a dt a la fin de la boucle en temps
    dt_crit = 2 * dt
    # ---------------------------------------------#
    #         MESH CREATION                        #
    # ---------------------------------------------#
    coord_mesh = np.loadtxt(meshfile, dtype=np.float64, skiprows=2, usecols=(1,))
    nodes_number = coord_mesh.shape[0]
    coord_init = np.zeros([nodes_number, 1], dtype=np.float64, order='C')
    coord_init[:, 0] = coord_mesh
    vit_init = np.zeros([nodes_number, 1], dtype=np.float64, order='C')
    enrichissement_registration = False
    if data.material_target.failure_model.type_of_enrichment == 'Hansbo':
        enrichissement_registration = True
    my_mesh = Mesh1dEnriched(
        initial_coordinates=coord_init, initial_velocities=vit_init,
        enrichment_type=data.material_target.failure_model.type_of_enrichment)

    # ---------------------------------------------#
    # Initialisation des vitesses des cibles et projectiles si elles sont dans le JDD
    # ---------------------------------------------#
    try:
        vitesse_projectile = data.material_projectile.initial_values.velocity_init
    except AttributeError:
        vitesse_projectile = 0
    try:
        vitesse_cible = data.material_target.initial_values.velocity_init
    except AttributeError:
    # le projectile n'est pas defini ou le champ vitesse initiale n'est pas defini
        vitesse_cible = 0

    print "Initilisation de la vitesse du projectile : {:} m/s".format(vitesse_projectile)
    my_mesh.nodes.upundemi[my_mesh.nodes.nodes_in_projectile, 0] = vitesse_projectile
    my_mesh.nodes.umundemi[my_mesh.nodes.nodes_in_projectile, 0] = vitesse_projectile
    print "Initilisation de la vitesse de la cible : {:} m/s".format(vitesse_cible)
    my_mesh.nodes.upundemi[my_mesh.nodes.nodes_in_target, 0] = vitesse_cible
    my_mesh.nodes.umundemi[my_mesh.nodes.nodes_in_target, 0] = vitesse_cible

    if data.geometric.initial_interface_position is not None:
        # Le noeud a l'interface entre les deux a comme vitesse une combinaison lineaire
        # entre les deux vitesses (Option ByVelocity plutot que ByMomentum)
        vitesse_interface = 0.5 * (vitesse_cible + vitesse_projectile)
        if (np.where(my_mesh.nodes.nodes_in_target)[0][0] ==
                np.where(my_mesh.nodes.nodes_in_projectile)[0][-1]):
            node_interface = np.where(my_mesh.nodes.nodes_in_target)[0][0]
        else:
            raise ValueError("""Impossible de trouver le noeud de l'interface """
                             """projectile / cible""")
        my_mesh.nodes.upundemi[node_interface, 0] = vitesse_interface
        my_mesh.nodes.umundemi[node_interface, 0] = vitesse_interface
        print ("Initilisation de la vitesse de l'interface : {:} m/s (noeud {:})"
               .format(vitesse_interface, node_interface))

    # ---------------------------------------------#
    #  FIGURES MANAGER SETUP                       #
    # ---------------------------------------------#
    the_figure_mng = FigureManager(my_mesh, dump=data.output.dump)
    if images_number != 0:
        the_figure_mng.set_time_controler(final_time / images_number)
        the_figure_mng.populate_figs()

    # ---------------------------------------------#
    #  OUTPUT MANAGER SETUP                        #
    # ---------------------------------------------#
    for db_el in data.output.databases:
        output_db = OutputDatabase(db_el.path)
        if db_el.iteration_period is not None:
            the_output_mng.register_database_iteration_ctrl(db_el.identifier, output_db,
                                                            db_el.iteration_period)
        else:
            the_output_mng.register_database_time_ctrl(db_el.identifier, output_db,
                                                       db_el.time_period)
        the_output_mng.register_all_fields(enrichissement_registration, my_mesh.cells,
                                           my_mesh.nodes, db_el.identifier,
                                           db_el.cell_indexes, db_el.node_indexes)
    # ---------------------------------------------#
    #         NODAL MASS COMPUTATION               #
    # ---------------------------------------------#
    my_mesh.compute_cells_sizes()
    my_mesh.compute_cells_masses()
    my_mesh.compute_nodes_masses()
    print "CALCULUS LAUNCHED!"
    compute_time = 0.

    # ************************************************* #
    #         DEBUT DE LA BOUCLE EN TEMPS               #
    # ************************************************* #
    while simulation_time < final_time:
        loop_begin_time = time.time()
        # if step % 1000 == 0:
        msg = ("""Iteration {:<4d} -- Time : {:15.9g} seconds with"""
               """ a time step of {:15.9g} seconds""".
               format(step, simulation_time, dt))
        print msg

        # ---------------------------------------------#
        #                OUTPUT MANAGEMENT             #
        # ---------------------------------------------#
        the_output_mng.update(simulation_time, step,
                              data.material_target.failure_model.type_of_enrichment,
                              data.material_target.failure_model.failure_treatment_value)
        #the_figure_mng.update(simulation_time, step)
        # commente en attendant de retrouver Qt => avec passage python3 + PyQt5

        # ---------------------------------------------#
        #         NODES VELOCITIES COMPUTATION         #
        # ---------------------------------------------#
        my_mesh.compute_new_nodes_velocities(dt_staggered)
        # ---------------------------------------------#
        #    APPLY VELOCITY BOUNDARY CONDITIONS        #
        # ---------------------------------------------#
        if left_boundary_condition.is_velocity():
            my_mesh.apply_velocity_boundary_condition(
                'left', left_boundary_condition.evaluate(simulation_time))
        if right_boundary_condition.is_velocity():
            my_mesh.apply_velocity_boundary_condition(
                'right', right_boundary_condition.evaluate(simulation_time))
        # ---------------------------------------------#
        #         CONTACT CORRECTION                   #
        # ---------------------------------------------#
        # my_mesh.compute_contact(dt)
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
        #    CELLS DEVIATOR STRESSES COMPUTATION       #
        # ---------------------------------------------#
        if data.material_projectile.constitutive_model.elasticity_model is not None:
            my_mesh.compute_deviator_elasticity(dt, my_mesh.cells.cell_in_projectile)
        if data.material_target.constitutive_model.elasticity_model is not None:
            my_mesh.compute_deviator_elasticity(dt, my_mesh.cells.cell_in_target)
        # ---------------------------------------------#
        #           PLASTICITY COMPUTATION             #
        # ---------------------------------------------#
        if data.material_projectile.constitutive_model.plasticity_model is not None:
            my_mesh.get_plastic_cells(projectile_plasticity_criterion,
                                      my_mesh.cells.cell_in_projectile)
        if data.material_target.constitutive_model.plasticity_model is not None:
            my_mesh.get_plastic_cells(target_plasticity_criterion,
                                      my_mesh.cells.cell_in_target)
        my_mesh.apply_plasticity_treatment(dt)
        # ---------------------------------------------#
        #         PSEUDOVISCOSITY COMPUTATION          #
        # ---------------------------------------------#
        my_mesh.compute_new_cells_pseudo_viscosity(dt)
        # ---------------------------------------------#
        #         CELLS PRESSURES COMPUTATION          #
        # ---------------------------------------------#
        my_mesh.compute_new_cells_pressures(dt)
        # ---------------------------------------------#
        #         STRESS TENSOR COMPUTATION            #
        # ---------------------------------------------#
        my_mesh.assemble_complete_stress_tensor()
        # ---------------------------------------------#
        #              RUPTURE                         #
        # ---------------------------------------------#
        if rupture_treatment is not None:
            my_mesh.get_ruptured_cells(rupture_criterion)
            my_mesh.apply_rupture_treatment(rupture_treatment, simulation_time)
        # ---------------------------------------------#
        #         NODES FORCES COMPUTATION             #
        # ---------------------------------------------#
        my_mesh.compute_new_nodes_forces()
        my_mesh.compute_new_cohesive_forces(time)
        # ---------------------------------------------#
        #         LOADING                              #
        # ---------------------------------------------#
        if left_boundary_condition.is_pressure():
            my_mesh.apply_pressure('left', left_boundary_condition.evaluate(simulation_time))
        if right_boundary_condition.is_pressure():
            my_mesh.apply_pressure('right', right_boundary_condition.evaluate(simulation_time))
        # ---------------------------------------------#
        #         TIME STEP COMPUTATION                #
        # ---------------------------------------------#
        dt_crit = my_mesh.compute_new_time_step()
        if dt != dt_crit:
            print "Reduction of the time step after failure. New time step is " + str(dt_crit)
        # reduction du pas de temps apres rupture
        dt = min(dt, dt_crit)  # pylint: disable=invalid-name

        # ---------------------------------------------#
        #                INCREMENTATION                #
        # ---------------------------------------------#
        my_mesh.increment()
        simulation_time += dt
        if not data.time.is_time_step_constant:
            dt_staggered = 0.5 * (dt_crit + dt)
            dt = dt_crit  # pylint: disable=invalid-name
        else:
            dt_staggered = dt
        step += 1
        loop_end_time = time.time()
        compute_time += loop_end_time - loop_begin_time

    print "Total time spent in compute operation is : {:15.9g} seconds".format(compute_time)
    plt.show(block=False)

    for item_time_data in history_list:
        item_time_data.write_fields_history()
        item_time_data.close_file()

    print 'Done !'

    the_output_mng.finalize()


if __name__ == "__main__":
    main()
