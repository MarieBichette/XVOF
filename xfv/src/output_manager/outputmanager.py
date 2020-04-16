"""
Implementing the OutputManager class
"""
from collections import namedtuple
import os
import numpy as np
from xfv.src.utilities.singleton import Singleton
from xfv.src.output_manager.outputtimecontroler import OutputTimeControler
from xfv.src.discontinuity.discontinuity import Discontinuity


DatabaseBuildInfos = namedtuple("DatabaseBuildInfos", ["database_object", "fields",
                                                       "time_controler"])
Field = namedtuple("Field", ["name", "owner", "attr_name", "indexes"])


class OutputManager(object, metaclass=Singleton):
    """
    The manager of all outputs
    """

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

    def register_field(self, field_name, field_support, field_attr_name, indexes=None,
                       database_names=None):
        """
        Add a field to the manager. Each field will be stored in every database in arguments
        if specified else in every database registered
        """
        if database_names is not None:
            for db in database_names:
                self.__db_build_infos[db].fields.append(
                    Field(name=field_name, owner=field_support, attr_name=field_attr_name,
                          indexes=indexes))
        else:
            for build_infos in list(self.__db_build_infos.values()):
                build_infos.fields.append(Field(name=field_name, owner=field_support,
                                                attr_name=field_attr_name, indexes=indexes))

    def register_all_fields(self, enrichment_registration, cells, nodes,
                            database_id):
        """
        Add all fields to the manager.
        :param enrichment_registration : bool to control if the
        :param cells: cells from which fields must be printed
        :param nodes: nodes from which fields must be printed
        :param database_id: identifier of the database
        """
        node_indexes = slice(0, nodes.number_of_nodes)
        cell_indexes = slice(0, cells.number_of_cells)
        enriched_cells = np.where(cells.enriched)[0]
        # Node field
        self.register_field("NodeStatus", nodes, ("enriched",),
                            database_names=[database_id], indexes=node_indexes)
        self.register_field("NodeCoordinates", nodes, ("xt",),
                            database_names=[database_id], indexes=node_indexes)
        self.register_field("ClassicalNodeVelocity", nodes, ("umundemi",),
                            database_names=[database_id], indexes=node_indexes)
        # Cell field
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
        self.register_field("ClassicalEquivalentPlasticStrainRate", cells,
                            ("equivalent_plastic_strain_rate",),
                            database_names=[database_id], indexes=cell_indexes)
        self.register_field("ClassicalPlasticStrainRate", cells, ("plastic_strain_rate",),
                            database_names=[database_id], indexes=cell_indexes)

        if enrichment_registration:
            # Enriched cell field -> cell support
            self.register_field("AdditionalPressure", cells,
                                ("additional_dof_pressure", "current_value"),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalDensity", cells,
                                ("additional_dof_density", "current_value"),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalInternalEnergy", cells,
                                ("additional_dof_energy", "current_value"),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalArtificialViscosity", cells,
                                ("additional_dof_artificial_viscosity", "current_value"),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalSoundVelocity", cells,
                                ("additional_dof_sound_velocity", "current_value"),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalLeftSize", cells, ("left_part_size", "current_value"),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalRightSize", cells, ("right_part_size", "current_value"),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalStress", cells, ("additional_dof_stress",),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalDeviatoricStress", cells,
                                ("additional_dof_deviatoric_stress_current",),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalEquivalentPlasticStrainRate", cells,
                                ("additional_dof_equivalent_plastic_strain_rate",),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalPlasticStrainRate", cells,
                                ("additional_dof_plastic_strain_rate",),
                                database_names=[database_id], indexes=enriched_cells)
            # Enriched node fields -> disc support
            self.register_field("AdditionalNodeVelocity", None,
                                ("additional_dof_velocity_current",),
                                database_names=[database_id])
            # Cohesive fields -> disc support
            self.register_field("AdditionalCohesiveForce", None,
                                ("cohesive_force", "current_value"),
                                database_names=[database_id])
            self.register_field("AdditionalDiscontinuityOpening", None,
                                ("discontinuity_opening", "current_value"),
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
        for build_infos in list(self.__db_build_infos.values()):
            if build_infos.time_controler.db_has_to_be_updated(time, iteration):
                build_infos.database_object.add_time(time)
                for field in build_infos.fields:
                    # Case : classical fields with cell or node support ------------------------
                    if not field.name.startswith("Additional") and field.indexes is not None:
                        value = self.get_value_of_field(field, field.owner)
                        build_infos.database_object.add_field(
                            field.name, value.__getitem__(field.indexes),
                            support=field.owner.__class__.__name__)
                    # Case : enriched cell field -----------------------------------------------
                    # Field registration in this case is a mix between cell classical fields and
                    # disc fields because of the refactoring of the enriched fields in Cell class

                    elif field.name.startswith("Additional") and field.owner is not None:
                        mask_enr = field.owner.enriched
                        if mask_enr.any():
                            cell_ids = np.where(mask_enr)[0]
                            value = self.get_value_of_field(field, field.owner)
                            if len(value.shape) == 1:
                                # Register cell scalar field
                                enr_field = np.array([cell_ids.tolist(),
                                                      value[mask_enr].tolist()]).transpose()

                            elif value.shape[1] == 3:
                                # Register cell tensor field
                                enr_field = np.array([cell_ids.tolist(),
                                                      value[mask_enr, 0].tolist(),
                                                      value[mask_enr, 1].tolist(),
                                                      value[mask_enr, 2].tolist()]).transpose()
                            else:
                                enr_field = []
                                raise ValueError(
                                    "Shape of the enriched cell field is not recognized")

                            build_infos.database_object.add_field(
                                field.name, enr_field,
                                support="Discontinuity", enrichment=type_of_enrichment,
                                discontinuity_position=eps)

                    # Case : enriched disc field (node and cohesive) -------------------------
                    elif field.name.startswith("Additional") and field.owner is None:
                        # Permet d'identifier les champs enrichis qui doivent se rapporter
                        # a un support disc
                        disc_field_collec = []
                        for disc in Discontinuity.discontinuity_list():
                            cell_id = np.where(disc.mask_ruptured_cell)[0][0]
                            value = self.get_value_of_field(field, disc)

                            if value.shape == (1, ):
                                # Register cell scalar field
                                disc_field_collec.append((cell_id, value[0]))

                            elif value.shape == (2, 1):
                                # Register node vector field
                                disc_field_collec.append((cell_id, value[0][0], value[1][0]))

                            elif value.shape == (1, 3):
                                # Register cell tensor field
                                disc_field_collec.append((cell_id, value[0, 0], value[0, 1],
                                                         value[0, 2]))
                            else:
                                raise ValueError("Unknown shape to register in database")

                        if len(Discontinuity.discontinuity_list()) > 0:
                            # sortir le build info de la boucle for des discontinuites permet
                            # d'enregistrer plusieurs discontinuites a la fois. Le seul "probleme"
                            # est que ces disc. ont toutes le meme support de type "Discontinuity"
                            # mais pas genant pour le moment a priori car on fait reference au
                            # nom de la classe et non a l'objet lui meme
                            build_infos.database_object.add_field(field.name,
                                                                  np.array(disc_field_collec),
                                                                  support="Discontinuity",
                                                                  enrichment=type_of_enrichment,
                                                                  discontinuity_position=eps)
                            # todo : faire passer la position de la disc. dans le tableau de
                            #  donnees car peut varier
                            # todo : d'une discontinuite a l'autre
                    # end enriched disc field -----------------------

    def get_value_of_field(self, field: Field, owner) -> np.array:
        """
        Get the np.array associated to the field following all the attribute names list
        """
        value = getattr(owner, field.attr_name[0])
        for attr_name in field.attr_name[1:]:
            value = getattr(value, attr_name).flatten()
        return value

    def finalize(self):
        """
        Close all the database
        """
        print("Flushing and writing database outputs!")
        for build_infos in list(self.__db_build_infos.values()):
            build_infos.database_object.close()

    def __str__(self):
        msg = self.__class__.__name__
        msg += " manages following database: " + os.linesep
        for db_name in self.__db_build_infos:
            msg += "|_> " + db_name + os.linesep
            msg += "    Fields : {:s}".format(",".join([field.name for field in
                                                        self.__db_build_infos[db_name].fields])) \
                   + os.linesep
        return msg
