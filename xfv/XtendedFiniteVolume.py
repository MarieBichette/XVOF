# XtendedFiniteVolume doesn't respect the style but we don't care pylint: disable=invalid-name
#!/usr/bin/env python2.7
"""
todo: to complete
"""
import os.path
import sys
import time
import matplotlib  # pylint: disable=unused-import
matplotlib.use('Qt4Agg')


import matplotlib.pyplot as plt
import numpy as np

from xfv.src.figure_manager.figure_manager      import FigureManager
from xfv.src.data.data_container                import DataContainer
from xfv.src.mesh.mesh1denriched                import Mesh1dEnriched
from xfv.src.data.save_time_data                import CellTimeData, NodeTimeData
from xfv.src.output_manager.outputmanager       import OutputManager
from xfv.src.output_manager.outputdatabase      import OutputDatabase
from xfv.src.rupturetreatment.enrichelement     import EnrichElement
from xfv.src.rupturetreatment.imposedpressure   import ImposedPressure


def __create_mesh(meshfile, enrichment_type):
    """
    Create a Mesh1D object from meshfile

    :param meshfile: path to the mesh file
    :param enrichment_type: type of enrichment desired
    """
    coord_mesh = np.loadtxt(meshfile, dtype=np.float64, skiprows=2, usecols=(1,))
    nodes_number = coord_mesh.shape[0]
    coord_init = np.zeros([nodes_number, 1], dtype=np.float64, order='C')
    coord_init[:, 0] = coord_mesh
    vit_init = np.zeros([nodes_number, 1], dtype=np.float64, order='C')
    return Mesh1dEnriched(initial_coordinates=coord_init, initial_velocities=vit_init,
                          enrichment_type=enrichment_type)


def __init_velocity(nodes, data):
    """
    Initialize the nodes velocity

    :param nodes: the nodes
    :type nodes: Node
    :param data: Simulation data
    :type data: DataContainer
    """
    try:
        vitesse_projectile = data.material_projectile.initial_values.velocity_init
    except AttributeError:
        vitesse_projectile = 0
    try:
        vitesse_cible = data.material_target.initial_values.velocity_init
    except AttributeError:
        vitesse_cible = 0

    print("Projectile velocity initialized to : {:} m/s".format(vitesse_projectile))
    nodes.upundemi[nodes.nodes_in_projectile, 0] = vitesse_projectile
    nodes.umundemi[nodes.nodes_in_projectile, 0] = vitesse_projectile
    print("Target velocity initialized to : {:} m/s".format(vitesse_cible))
    nodes.upundemi[nodes.nodes_in_target, 0] = vitesse_cible
    nodes.umundemi[nodes.nodes_in_target, 0] = vitesse_cible

    if data.geometric.initial_interface_position is not None:
        # Need to find the node at the interface between projectile and target
        if (np.where(nodes.nodes_in_target)[0][0] ==
                np.where(nodes.nodes_in_projectile)[0][-1]):
            node_interface = np.where(nodes.nodes_in_target)[0][0]
        else:
            raise ValueError("""Unable to find the node at the interface"""
                             """projectile / target""")
        # Node at the interface between projectile and target is initialized
        # with a linear combination of both velocities
        vitesse_interface = 0.5 * (vitesse_cible + vitesse_projectile)
        nodes.upundemi[node_interface, 0] = vitesse_interface
        nodes.umundemi[node_interface, 0] = vitesse_interface
        print(("Interface velocity initialized to : {:} m/s (node {:})"
               .format(vitesse_interface, node_interface)))


def __init_time_figure_plot(output_data, path):
    """
    Initialize the time figure plots

    :param output_data: output data
    :type output_data: data_container.output_props
    :param path: path toward the directory holding data file
    :type path: str
    :return: the list of nodes or cells that must be followed
    :rtype: List[TimeData]
    """
    cells_for_time_figure = []
    nodes_for_time_figure = []
    try:
        for cell_number in output_data.cells_numbers:
            cells_for_time_figure.append(int(cell_number))
    except TypeError:
        pass

    try:
        for node_number in output_data.nodes_numbers:
            nodes_for_time_figure.append(int(node_number))
    except TypeError:
        pass
    h_list = [CellTimeData(cell_id, path) for cell_id in cells_for_time_figure]
    h_list += [NodeTimeData(node_id, path) for node_id in nodes_for_time_figure]
    return h_list


def __init_output(output_data, type_of_enrichment, mesh):
    """
    Returns the OutputManager initialized

    :param output_data: output data
    :type output_data: data_container.output_props
    :param type_of_enrichment: type of enrichment
    :type type_of_enrichment: str
    :param mesh: the mesh
    :type mesh: Mesh1DEnriched
    """
    enrichment_registration = type_of_enrichment == 'Hansbo'
    np.set_printoptions(formatter={'float': '{: 25.23g}'.format})
    the_output_mng = OutputManager()
    for db_el in output_data.databases:
        output_db = OutputDatabase(db_el.path)
        if db_el.iteration_period is not None:
            the_output_mng.register_database_iteration_ctrl(db_el.identifier, output_db,
                                                            db_el.iteration_period)
        else:
            the_output_mng.register_database_time_ctrl(db_el.identifier, output_db,
                                                       db_el.time_period)
        the_output_mng.register_all_fields(enrichment_registration, mesh.cells,
                                           mesh.nodes, db_el.identifier,
                                           db_el.cell_indexes, db_el.node_indexes)
    return the_output_mng


def main():  #pylint: disable=too-many-locals, too-many-branches, too-many-statements
    """
    Launch the program
    """
    # ------------------------------------------------------------------
    #             PARAMETERS INITIALIZATION
    # ------------------------------------------------------------------
    # ---- # DATA FILES
    path = sys.argv[1]
    data = DataContainer(os.path.join(path, "XDATA.xml"))
    meshfile = os.path.join(path, "mesh.txt")
    print("Running simulation for {:s}".format(os.path.normpath(os.path.abspath(path))))

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
        print("Enrichment method : {}".format(type_of_enrichment))
    # ---------------------------------------------#
    #         MESH CREATION                        #
    # ---------------------------------------------#
    type_of_enr = data.material_target.failure_model.type_of_enrichment
    my_mesh = __create_mesh(meshfile, type_of_enr)

    # ---------------------------------------------#
    # TARGET AND PROJECTILE VELOCITIES INITIALIZATION
    # ---------------------------------------------#
    __init_velocity(my_mesh.nodes, data)

    # ---------------------------------------------#
    #  FIGURES MANAGER SETUP                       #
    # ---------------------------------------------#
    the_figure_mng = FigureManager(my_mesh, dump=data.output.dump)
    images_number = data.output.number_of_images
    if images_number != 0:
        the_figure_mng.set_time_controler(final_time / images_number)
        the_figure_mng.populate_figs()

    # ---------------------------------------------#
    #  OUTPUT MANAGER SETUP                        #
    # ---------------------------------------------#
    history_list = __init_time_figure_plot(data.output, path)
    the_output_mng = __init_output(data.output, type_of_enr, my_mesh)
    # ---------------------------------------------#
    #         NODAL MASS COMPUTATION               #
    # ---------------------------------------------#
    my_mesh.compute_cells_sizes()
    my_mesh.compute_cells_masses()
    my_mesh.compute_nodes_masses()
    print("CALCULUS LAUNCHED!")
    compute_time = 0.

    # ************************************************* #
    #         DEBUT DE LA BOUCLE EN TEMPS               #
    # ************************************************* #
    simulation_time = 0.
    step = 0
    dt = initial_time_step  # pylint: disable=invalid-name
    dt_staggered = dt / 2
    # Le premier increment pour la vitesse a un pas de temps de dt/2 pour
    # tenir compte du fait que les vitesses sont
    # init a t=0 et pas t = -1/2 dt (pour garder l'ordre 2 a l'init).
    # Le dt_staggered est ensuite remis a dt a la fin de la boucle en temps
    dt_crit = 2 * dt
    while simulation_time < final_time:
        loop_begin_time = time.time()
        # if step % 1000 == 0:
        msg = ("""Iteration {:<4d} -- Time : {:15.9g} seconds with"""
               """ a time step of {:15.9g} seconds""".
               format(step, simulation_time, dt))
        print(msg)

        # ---------------------------------------------#
        #                OUTPUT MANAGEMENT             #
        # ---------------------------------------------#
        the_output_mng.update(simulation_time, step,
                              data.material_target.failure_model.type_of_enrichment,
                              data.material_target.failure_model.failure_treatment_value)
        the_figure_mng.update(simulation_time, step)

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
        target_model = data.material_target.constitutive_model
        projectile_model = None
        if data.data_contains_a_projectile:
            projectile_model = data.material_projectile.constitutive_model
        if projectile_model and projectile_model.elasticity_model is not None:
            my_mesh.compute_deviator_elasticity(dt, my_mesh.cells.cell_in_projectile)
        if target_model.elasticity_model is not None:
            my_mesh.compute_deviator_elasticity(dt, my_mesh.cells.cell_in_target)
        # ---------------------------------------------#
        #           PLASTICITY COMPUTATION             #
        # ---------------------------------------------#
        if projectile_model and projectile_model.plasticity_model is not None:
            my_mesh.get_plastic_cells(projectile_model.plasticity_criterion,
                                      my_mesh.cells.cell_in_projectile)
        if target_model.plasticity_model is not None:
            my_mesh.get_plastic_cells(target_model.plasticity_criterion,
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
        my_mesh.compute_new_cohesive_forces()
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
        dt = min(dt, dt_crit)  # pylint: disable=invalid-name
        if dt != dt_crit:
            print("Reduction of the time step after failure. New time step is " + str(dt_crit))

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

    print("Total time spent in compute operation is : {:15.9g} seconds".format(compute_time))
    plt.show(block=False)

    for item_time_data in history_list:
        item_time_data.write_fields_history()
        item_time_data.close_file()

    print('Done !')

    the_output_mng.finalize()


if __name__ == "__main__":
    main()
