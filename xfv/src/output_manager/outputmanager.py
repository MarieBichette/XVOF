"""
Implementing the OutputManager class
"""
from collections import namedtuple
import os
import numpy as np
from xfv.src.utilities.singleton import Singleton
from xfv.src.output_manager.outputtimecontroler import OutputTimeControler
from xfv.src.discontinuity.discontinuity import Discontinuity


DatabaseBuildInfos = namedtuple("DatabaseBuildInfos", ["database_object", "fields", "time_controler"])
Field = namedtuple("Field", ["name", "owner", "attr_name", "indexes"])


class OutputManager(object):
    """
    The manager of all outputs
    """
    __metaclass__ = Singleton

    def __init__(self):
        self.__db_build_infos = {}

    def register_database_time_ctrl(self, database_name, database_obj, delta_t):
        """
        Add a database to the manager. The database will be updated every deltat_t seconds
        """
        time_ctrl = OutputTimeControler(database_name, time_period=delta_t)
        dbinfos = DatabaseBuildInfos(database_obj, fields=[], time_controler=time_ctrl)
        self.__db_build_infos[database_name] = dbinfos

    def register_database_iteration_ctrl(self, database_name, database_obj, delta_it):
        """
        Add a database to the manager. The database will be updated every delta_it iterations
        """
        time_ctrl = OutputTimeControler(database_name, iteration_period=delta_it)
        dbinfos = DatabaseBuildInfos(database_obj, fields=[], time_controler=time_ctrl)
        self.__db_build_infos[database_name] = dbinfos

    def register_field(self, field_name, field_support, field_attr_name, indexes=None, database_names=None):
        """
        Add a field to the manager. Each field will be stored in every database in arguments
        if specified else in every database registered
        """
        if database_names is not None:
            for db in database_names:
                self.__db_build_infos[db].fields.append(
                    Field(name=field_name, owner=field_support, attr_name=field_attr_name, indexes=indexes))
        else:
            for build_infos in self.__db_build_infos.values():
                build_infos.fields.append(Field(name=field_name, owner=field_support,
                                                attr_name=field_attr_name, indexes=indexes))

    def register_all_fields(self, enrichment_registration, cells, nodes,
                            database_id, cell_indexes=None, node_indexes=None, disc_indexes=None):
        """
        Add all fields to the manager.
        :param enrichment_registration : bool to control if the
        :param cells: cells from which fields must be printed
        :param nodes: nodes from which fields must be printed
        :param database_id: identifier of the database
        :param cell_indexes: indexes of the cells to be printed
        :param node_indexes: indexes of the nodes to be printed
        :param disc_indexes : indexes of the discontinuity to be printed
        """
        if cell_indexes is None:
            cell_indexes = slice(0, cells.number_of_cells)
        if node_indexes is None:
            node_indexes = slice(0, nodes.number_of_nodes)
        self.register_field("NodeStatus", nodes, ("enriched",),
                            database_names=[database_id], indexes=node_indexes)
        self.register_field("NodeCoordinates", nodes, ("xt",),
                            database_names=[database_id], indexes=node_indexes)
        self.register_field("ClassicalNodeVelocity", nodes, ("umundemi",),
                            database_names=[database_id], indexes=node_indexes)

        self.register_field("CellStatus", cells, ("enriched",),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("CellSize", cells, ("size_t_plus_dt",),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalPressure", cells, ("pressure", "current_value"),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalDensity", cells, ("density", "current_value"),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalInternalEnergy", cells, ("energy", "current_value"),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalSoundVelocity", cells, ("sound_velocity", "current_value"),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalArtificialViscosity", cells, ("pseudo", "current_value"),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalStress", cells, ("stress", ),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalDeviatoricStress", cells, ("deviatoric_stress_current",),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalEquivalentPlasticStrainRate", cells, ("equivalent_plastic_strain_rate",),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalPlasticStrainRate", cells, ("plastic_strain_rate",),
                            database_names=[database_id], indexes=cell_indexes)

        if enrichment_registration:
            self.register_field("AdditionalPressure", None, ("additional_dof_pressure", "current_value"),
                                database_names=[database_id])
            self.register_field("AdditionalDensity", None, ("additional_dof_density", "current_value"),
                                database_names=[database_id])
            self.register_field("AdditionalInternalEnergy", None, ("additional_dof_energy", "current_value"),
                                database_names=[database_id])
            self.register_field("AdditionalArtificialViscosity", None,
                                ("additional_dof_artificial_viscosity", "current_value"),
                                database_names=[database_id])
            self.register_field("AdditionalSoundVelocity", None, ("additional_dof_sound_velocity", "current_value"),
                                database_names=[database_id])
            self.register_field("AdditionalLeftSize", None, ("left_part_size", "current_value"),
                                database_names=[database_id])
            self.register_field("AdditionalRightSize", None, ("right_part_size", "current_value"),
                                database_names=[database_id])
            self.register_field("AdditionalNodeVelocity", None, ("additional_dof_velocity_current",),
                                database_names=[database_id])
            self.register_field("AdditionalStress", None, ("additional_dof_stress",),
                                database_names=[database_id])
            self.register_field("AdditionalDeviatoricStress", None, ("additional_dof_deviatoric_stress_current",),
                                database_names=[database_id])
            self.register_field("AdditionalEquivalentPlasticStrainRate", None, ("additional_dof_equivalent_plastic_strain_rate",),
                                database_names=[database_id])
            self.register_field("AdditionalPlasticStrainRate", None, ("additional_dof_plastic_strain_rate",),
                                database_names=[database_id])
            self.register_field("AdditionalCohesiveForce", None, ("cohesive_force", "current_value"),
                                database_names=[database_id])
            self.register_field("AdditionalDiscontinuityOpening", None, ("discontinuity_opening", "current_value"),
                                database_names=[database_id])

    def update(self, time, iteration, type_of_enrichment, eps):
        """
        If the current time given in argument is above the time of next output then
        the manager asks each of its database to save fields. It's the same for
        a current iteration above the iteration of next output

        Additional_dof_fields are created when a new discontinuity is created.
        Need to treat them in a different way from classical fields.
        Based on this remark, additional_dof_fields have standard name "Additional..."
        Differentiation is made with test startswith(Additional)
        """
        for build_infos in self.__db_build_infos.values():
            if build_infos.time_controler.db_has_to_be_updated(time, iteration):
                build_infos.database_object.add_time(time)
                for field in build_infos.fields:
                    if not field.name.startswith("Additional") and field.indexes is not None:
                        value = getattr(field.owner, field.attr_name[0])
                        for attr_name in field.attr_name[1:]:
                            value = getattr(value, attr_name)
                        build_infos.database_object.add_field(
                            field.name, value.__getitem__(field.indexes), support=field.owner.__class__.__name__)
                    elif field.name.startswith("Additional"):
                        # Permet d'identifier les champs enrichis qui doivent se rapporter a un support disc
                        additional_field = []
                        for disc in Discontinuity.discontinuity_list():
                            cell_id = np.where(disc.mask_ruptured_cell)[0][0]
                            value = getattr(disc, field.attr_name[0])
                            for attr_name in field.attr_name[1:]:
                                value = getattr(value, attr_name).flatten()
                            try:  # enregistre un tenseur diagonal (3 valeurs)
                                additional_field.append((cell_id, value[:,0], value[:,1], value[:,2]))
                            except IndexError:
                                try:  # enrgistre un vecteur vitesse (2 composantes)
                                    additional_field.append((cell_id, value[0], value[1]))
                                except IndexError:  # enregistre un scalaire (champs thermo)
                                    additional_field.append((cell_id, value))
                        if len(Discontinuity.discontinuity_list()) > 0:
                            # sortir le build info de la boucle for des discontinuites permet d'enregistrer plusieurs
                            # discontinuites a la fois. Le seul "probleme" est que ces disc. ont toutes le meme support
                            # de type "Discontinuity" mais pas genant pour le moment a priori car on fait reference au
                            # nom de la classe et non a l'objet lui meme
                            build_infos.database_object.add_field(field.name, np.array(additional_field),
                                                                  support="Discontinuity",
                                                                  enrichment=type_of_enrichment,
                                                                  discontinuity_position=eps)
                            # todo : faire passer la position de la disc. dans le tableau de donnees car peut varier
                            # todo : d'une discontinuite a l'autre


    def finalize(self):
        """
        Close all the database
        """
        print("Flushing and writing database outputs!")
        for build_infos in self.__db_build_infos.values():
            build_infos.database_object.close()

    def __str__(self):
        msg = self.__class__.__name__
        msg += " manages following database: " + os.linesep
        for db_name in self.__db_build_infos:
            msg += "|_> " + db_name + os.linesep
            msg += "    Fields : {:s}".format(
                ",".join([field.name for field in self.__db_build_infos[db_name].fields])) + os.linesep
        return msg
