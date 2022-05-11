"""
Implementing the OutputManager class
"""
from collections import namedtuple
import os
import numpy as np
from xfv.src.utilities.singleton import Singleton
from xfv.src.output_manager.outputtimecontroler import OutputTimeControler
from xfv.src.data.data_container import DataContainer


DatabaseBuildInfos = namedtuple("DatabaseBuildInfos", ["database_object", "fields",
                                                       "time_controler"])
Field = namedtuple("Field", ["name", "owner", "attr_name", "indexes"])

FieldConstruction = namedtuple("FieldConstruction", ["name", "support", "attr_name"])

field_list = dict()
field_list["NodeVelocity"] = FieldConstruction("ClassicalNodeVelocity", "nodes", ("umundemi",))
field_list["NodeCoordinates"] = FieldConstruction("NodeCoordinates", "nodes", ("xt",))
field_list["CellSize"] = FieldConstruction("CellSize", "cells", ("size_t_plus_dt",))
# TODO corriger ça en size t
field_list["Pressure"] = FieldConstruction("ClassicalPressure",
                                           "cells", ("pressure", "current_value"))
field_list["Density"] = FieldConstruction("ClassicalDensity", "cells", ("density", "current_value"))
field_list["InternalEnergy"] = FieldConstruction("ClassicalInternalEnergy",
                                                 "cells", ("energy", "current_value"))
field_list["SoundVelocity"] = FieldConstruction("ClassicalSoundVelocity",
                                                "cells", ("sound_velocity", "current_value"))
field_list["ArtificialViscosity"] = FieldConstruction("ClassicalArtificialViscosity",
                                                      "cells", ("pseudo", "current_value"))
field_list["Stress"] = FieldConstruction("ClassicalStress", "cells", ("stress", ))
field_list["DeviatoricStress"] = FieldConstruction("ClassicalDeviatoricStress",
                                                   "cells", ("deviatoric_stress_current", ))
field_list["EquivalentPlasticStrainRate"] = FieldConstruction(
    "ClassicalEquivalentPlasticStrainRate", "cells", ("equivalent_plastic_strain_rate", ))
field_list["PlasticStrainRate"] = FieldConstruction("ClassicalPlasticStrainRate",
                                                    "cells", ("plastic_strain_rate", ))
field_list["Porosity"] = FieldConstruction("ClassicalPorosity",
                                           "cells", ("porosity", "current_value"))
field_list["ShearModulus"] = FieldConstruction("ClassicalShearModulus",
                                               "cells", ("shear_modulus", "current_value"))
field_list["YieldStress"] = FieldConstruction("ClassicalYieldStress",
                                              "cells", ("yield_stress", "current_value"))

enr_field_list = dict()
enr_field_list["NodeVelocity"] = FieldConstruction("AdditionalNodeVelocity", None,
                                                   ("enr_velocity_current", ))
enr_field_list["NodeCoordinates"] = FieldConstruction("AdditionalNodeCoordinates", None,
                                                      ("enr_coordinates_current", ))
enr_field_list["CohesiveForce"] = FieldConstruction("AdditionalCohesiveForce", None,
                                                    ("cohesive_force", "current_value"))
enr_field_list["DiscontinuityOpening"] = FieldConstruction("AdditionalDiscontinuityOpening",
                                                           None, ("discontinuity_opening",
                                                                  "current_value"))
enr_field_list["DissipatedEnergy"] = FieldConstruction("AdditionalDissipatedEnergy",
                                                           None, ("dissipated_energy",
                                                                  "current_value"))
enr_field_list["Pressure"] = FieldConstruction("AdditionalPressure", "cells",
                                               ("enr_pressure", "current_value"))
enr_field_list["Density"] = FieldConstruction("AdditionalDensity", "cells",
                                              ("enr_density", "current_value"))
enr_field_list["InternalEnergy"] = FieldConstruction("AdditionalInternalEnergy", "cells",
                                                     ("enr_energy", "current_value"))
enr_field_list["ArtificialViscosity"] = FieldConstruction("AdditionalArtificialViscosity", "cells",
                                                          ("enr_artificial_viscosity",
                                                           "current_value"))
enr_field_list["SoundVelocity"] = FieldConstruction("AdditionalSoundVelocity", "cells",
                                                    ("enr_sound_velocity", "current_value"))
enr_field_list["Stress"] = FieldConstruction("AdditionalStress", "cells", ("enr_stress", ))
enr_field_list["DeviatoricStress"] = FieldConstruction("AdditionalDeviatoricStress", "cells",
                                                       ("enr_deviatoric_stress_current", ))
enr_field_list["EquivalentPlasticStrainRate"] = FieldConstruction(
    "AdditionalEquivalentPlasticStrainRate", "cells", ("enr_equivalent_plastic_strain_rate", ))
enr_field_list["PlasticStrainRate"] = FieldConstruction(
    "AdditionalPlasticStrainRate", "cells", ("enr_plastic_strain_rate", ))
enr_field_list["Porosity"] = FieldConstruction("AdditionalPorosity", "cells",
                                               ("enr_porosity", "current_value"))
field_list["ShearModulus"] = FieldConstruction("AdditionalShearModulus",
                                               "cells", ("enr_shear_modulus", "current_value"))
field_list["YieldStress"] = FieldConstruction("AdditionalYieldStress",
                                              "cells", ("enr_yield_stress", "current_value"))

class OutputManager(metaclass=Singleton):
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

        :param enrichment_registration: bool to control if the enriched fields should be registered
        :param cells: cells from which fields must be printed
        :param nodes: nodes from which fields must be printed
        :param database_id: identifier of the database
        """
        node_indexes = slice(0, nodes.number_of_nodes)
        cell_indexes = slice(0, cells.number_of_cells)
        enriched_cells = np.where(cells.enriched)[0]

        # Node and cell status should always been registered
        self.register_field("NodeStatus", nodes, ("enriched",),  # should always be registered
                            database_names=[database_id], indexes=node_indexes)
        self.register_field("CellStatus", cells, ("enriched",),  # should always be registered
                            database_names=[database_id], indexes=cell_indexes)

        # Classical field registration
        for key in field_list:
            if key in DataContainer().output.variables:  # registration if field is in the dataset
                field_infos = field_list[key]
                # Node field
                if field_infos.support == "nodes":
                    self.register_field(field_infos.name, nodes, field_infos.attr_name,
                                        database_names=[database_id], indexes=node_indexes)
                # Cell field
                if field_infos.support == "cells":
                    self.register_field(field_infos.name, cells, field_infos.attr_name,
                                        database_names=[database_id], indexes=cell_indexes)

        if enrichment_registration:
            # => Register the enriched field also
            # Left and right size should always been registered
            self.register_field("AdditionalLeftSize", cells, ("left_part_size", "current_value"),
                                database_names=[database_id], indexes=enriched_cells)
            self.register_field("AdditionalRightSize", cells, ("right_part_size", "current_value"),
                                database_names=[database_id], indexes=enriched_cells)

            for key in enr_field_list:  # registration if field is in the dataset
                if key in DataContainer().output.variables:
                    field_infos = enr_field_list[key]
                    # Enr Node field
                    if field_infos.support is None:
                        self.register_field(field_infos.name, None, field_infos.attr_name,
                                            database_names=[database_id])
                    # Enr Cell field
                    if field_infos.support == "cells":
                        self.register_field(field_infos.name, cells, field_infos.attr_name,
                                            database_names=[database_id], indexes=cell_indexes)

    def update(self, time, iteration, eps, discontinuity_list):
        """
        If the current time given in argument is above the time of next output then
        the manager asks each of its database to save fields. It's the same for
        a current iteration above the iteration of next output

        Additional_dof_fields are created when a new discontinuity is created.
        Need to treat them in a different way from classical fields.
        Based on this remark, enr_fields have standard name "Additional..."
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
                                support="Discontinuity", enrichment="Hansbo",
                                discontinuity_position=eps)

                    # Case : enriched disc field (node and cohesive) -------------------------
                    elif field.name.startswith("Additional") and field.owner is None:
                        # Permet d'identifier les champs enrichis qui doivent se rapporter
                        # a un support disc
                        disc_field_collec = []
                        for disc in discontinuity_list:
                            cell_id = disc.get_ruptured_cell_id
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

                        if len(discontinuity_list) > 0:
                            # sortir le build info de la boucle for des discontinuites permet
                            # d'enregistrer plusieurs discontinuites a la fois. Le seul "probleme"
                            # est que ces disc. ont toutes le meme support de type "Discontinuity"
                            # mais pas genant pour le moment a priori car on fait reference au
                            # nom de la classe et non a l'objet lui meme
                            build_infos.database_object.add_field(field.name,
                                                                  np.array(disc_field_collec),
                                                                  support="Discontinuity",
                                                                  enrichment="Hansbo",
                                                                  discontinuity_position=eps)
                            # todo : faire passer la position de la disc. dans le tableau de
                            #  donnees car peut varier
                            # todo : d'une discontinuite a l'autre
                    # end enriched disc field -----------------------

    def get_value_of_field(self, field: Field, owner) -> np.array:
        """
        Get the np.array associated to the field following all the attribute names list

        :param field: field to be extracted
        :param owner: object who supports the fields
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
