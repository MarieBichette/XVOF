#!/usr/bin/env python

from distutils.core import setup

setup(name='XFV',
      version='1.0',
      description='A 1D hydrocode with discontinuous fields abilities',
      author='Guillaume Peillex',
      author_email='guillaume.peillex@gmail.com',
      maintainer='Guillaume Peillex',
      maintainer_email='guillaume.peillex@gmail.com',
      url='https://github.com/hippo91/XVOF',
      packages=['xvof.src', 'xvof.src.boundary_condition', 'xvof.src.cell', 'xvof.src.cohesive_model',
                'xvof.src.contact', 'xvof.src.data', 'xvof.src.discontinuity', 'xvof.src.equationsofstate',
                'xvof.executables', 'xvof.src.fields', 'xvof.src.figure_manager', 'xvof.src.mass_matrix',
                'xvof.src.mesh', 'xvof.src.node', 'xvof.src.output_figure', 'xvof.src.output_manager',
                'xvof.src.plasticitycriterion', 'xvof.src.pressurelaw', 'xvof.src.rheology', 'xvof.src.rupturecriterion',
                'xvof.src.rupturetreatment', 'xvof.src.solver', 'xvof.src.utilities'],
        scripts=['xvof/executables/VNR_1D_EnrichedElement.py']
     )