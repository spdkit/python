.. spdkit documentation master file, created by
   sphinx-quickstart on Wed Feb 15 20:45:23 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to spdkit's documentation!
==================================
.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: spdkit
   :members:

Chemical objects
==================================
.. autoclass:: spdkit.Molecule
   :members:

.. autoclass:: spdkit.Atom
    :members:

.. autoclass:: spdkit.Lattice
    :members:

Gosh
==================================
.. autoclass:: spdkit.gosh.BlackBoxModel
   :members:

.. autoclass:: spdkit.gosh.Computed
   :members:

.. autoclass:: spdkit.gosh.JobHub
   :members:

.. autoclass:: spdkit.gosh.DbConnection
   :members:

.. autofunction:: spdkit.gosh.optimize


Read/write
==================
.. autoclass:: spdkit.io.GrepReader
   :members:

.. autoclass:: spdkit.io.TextViewer
   :members:

.. autoclass:: spdkit.io.OptimizationTrajactory
   :members:

.. autoclass:: spdkit.io.Template
   :members:

.. autofunction:: spdkit.io.read

.. autofunction:: spdkit.io.write

.. autofunction:: spdkit.io.guess_format_from_path

Utils
==================
.. autofunction:: spdkit.utils.jmol_selection_commands

.. autofunction:: spdkit.utils.set_verbosity

.. autofunction:: spdkit.utils.abbreviate_numbers_human_readable

.. autofunction:: spdkit.utils.parse_numbers_human_readable

Surface
==================
.. autofunction:: spdkit.surface.probe_surface_atoms
.. autofunction:: spdkit.surface.probe_adsorption_sites


Reaction Path
==================
.. autoclass:: spdkit.apps.ReactionPreview
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
