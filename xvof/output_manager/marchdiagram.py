#!/usr/bin/env python2.7
"""
A script writing a march diagram after exploitation of the output database
"""
from __future__ import print_function
import numpy as np
from xvof.output_manager.outputdatabaseexploit import OutputDatabaseExploit


if __name__ == "__main__":
    my_hd = OutputDatabaseExploit("all_fields.hdf5")
    with open("/tmp/toto.txt", "w") as fo:
        for t in my_hd.saved_times:
            values = my_hd.extract_true_field_at_time("Pressure", t)
            t_arr = np.ndarray((values.shape[0], 1), dtype=np.float64, order="C")
            t_arr[:] = t
            v = np.concatenate((t_arr, values), axis=1)
            np.savetxt(fo, v)
            print("\n", file=fo)
