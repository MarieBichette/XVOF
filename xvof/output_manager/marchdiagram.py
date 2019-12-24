#!/usr/bin/env python2.7
"""
A script writing a march diagram after exploitation of the output database
"""
from __future__ import print_function
import numpy as np
from xvof.output_manager.outputdatabaseexploit import OutputDatabaseExploit


if __name__ == "__main__":
        my_hd = OutputDatabaseExploit("all_fields.hdf5")
        with open("/tmp/titi.txt", "w") as fo:
            for t in my_hd.saved_times:
                values = my_hd.extract_true_field_at_time("EquivalentPlasticStrainRate", t)
                t_arr = np.ndarray((values.shape[0]), dtype=np.float64, order="C")
                t_arr[:] = t
                v = np.ndarray((values.shape[0], 3), dtype=np.float64, order="C")
                v[:, 0] = values[:, 0]  # position
                v[:, 1] = t_arr  # temps
                v[:, 2] = values[:, 1]  # variable
                np.savetxt(fo, v)
                print("", file=fo)
