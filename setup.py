#!/usr/bin/env python

import aempy3

# Check for setuptools package:

try:
    from setuptools import setup
except ImportError:
    setuptools = False
    from distutils.core import setup
else:
    setuptools = True

LONG_DESC = """
aempy3 is an open source Python package for Airborne Electromagnetics (AEM)
data processing, analysis, modelling, visualization and interpretation.
"""

setup_kwargs = {}

# The advantage of setuptools is that EXE wrappers are created on Windows,
# which allows Tab-completion for the script names from the system Scripts
# folder.

# Add names of scripts here. You can also specify the function to call
# by adding :func_name after the module name, and the name of the script
# can be customized before the equals sign.

# setup_kwargs['entry_points'] = {'console_scripts':
#                     ['ws2vtk = aempy3.utils.ws2vtk:main',
#                      'modem_pyqt = aempy3.gui.modem_pyqt:main',
#                      'modem_plot_response = aempy3.gui.modem_plot_response:main',
#                      'modem_plot_pt_maps = aempy3.gui.modem_plot_pt_maps:main',
#                      'modem_mesh_builder = aempy3.gui.modem_mesh_builder:main',
#                      'modem2vtk = aempy3.utils.modem2vtk:main',
#                      'occam1d_gui = aempy3.gui.occam1d_gui:main',
#                      'edi_editor = aempy3.gui.edi_editor:main']}

# But many people will not have setuptools installed, so we need to handle
# the default Python installation, which only has Distutils:

# if setuptools is False:
    # Different script specification style for ordinary Distutils:

    setup_kwargs['scripts'] = [
        s.split(' = ')[1].replace('.', '/').split(':')[0] + '.py' for s in
        setup_kwargs['entry_points']['console_scripts']]
    del setup_kwargs['entry_points']

    # "You must explicitly list all packages in packages: the Distutils will not
    # recursively scan your source tree looking for any directory with an
    # __init__.py file"

setup_kwargs['packages'] = [
                            'aempy3',
                            'aempy3.core',
                            'aempy3.modules',
                            'aempy3.scripts',
                            'aempy3.docs',
                            'aempy3.examples'
]

setup_kwargs['install_requires'] = ['numpy>=1.17',
                                     'scipy>=1.4',
                                     'matplotlib>=3,
                                     'pyyaml>=5',
                                     'geopandas>=6']

setup(
	name="aempy3",
	version=aempy3.__version__,
	author="Duygu Kiyan, Volker Rath",
    author_email="volker.rath.dublin@gmail.com",
	description="Python toolkit for airborne EM data processing and inversion.",
	long_description=LONG_DESC,
    url="https://github.com/volkerrath/AEMPY33",
	include_package_data=True,
	license="GNU GENERAL PUBLIC LICENSE v3",
	classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        ],
	**setup_kwargs)
