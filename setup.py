#!/usr/bin/env python
"""
Setup script for XtendedFiniteVolume project
"""
from setuptools import setup, find_packages


setup(name='XtendedFiniteVolume',
      version='1.0',
      description='A 1D hydrocode with discontinuous fields abilities',
      author='Guillaume Peillex',
      author_email='guillaume.peillex@gmail.com',
      maintainer='Guillaume Peillex',
      maintainer_email='guillaume.peillex@gmail.com',
      url='https://github.com/hippo91/XVOF',
      packages=find_packages(),
      scripts=['xfv/XtendedFiniteVolume.py'],
      install_requires=['h5py',
                        'lxml',
                        'matplotlib',
                        'numpy']
     )
