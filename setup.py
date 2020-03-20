#!/usr/bin/env python
"""
Setup script for XtendedFiniteVolume project
"""

from distutils.core import setup
import os
import os.path


def find_packages():
    """
    Return the list of all packages
    """
    packages = []
    for (dirpath, _, _) in os.walk('xfv'):
        if os.path.exists(os.path.join(dirpath, '__init__.py')):
            packages.append(dirpath.replace('/', '.'))
    return packages


setup(name='XtendedFiniteVolume',
      version='1.0',
      description='A 1D hydrocode with discontinuous fields abilities',
      author='Guillaume Peillex',
      author_email='guillaume.peillex@gmail.com',
      maintainer='Guillaume Peillex',
      maintainer_email='guillaume.peillex@gmail.com',
      url='https://github.com/hippo91/XVOF',
      packages=find_packages(),
      scripts=['xfv/XtendedFiniteVolume.py']
     )
