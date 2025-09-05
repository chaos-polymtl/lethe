==============================
Command-line Interface Options
==============================

--------------------------------
Print deal.II and Lethe Versions
--------------------------------
The flag ``-V`` allos the user to print the deal.II and Lethe versions used to
run an application. In an example of usage, ``lethe-fluid -V`` yields:

.. code-block:: shell
   
   Running: lethe-fluid -V
   lethe/1.0-77-g9fe541d99 deal.II/9.6.2-2966-g7a74e69540

.. hint::

   The ``-V`` tag can be used with all Lethe applications.


----------------------------
Delete Previous Output Files
----------------------------
The flag ``-R`` allows the user to delete previous output files before 
new ones are written. This is useful, for instance, when launching the simulation
multiple times with different numbers of steps, so that you make sure only the 
most recent output files are stored. One example is as follows:

.. code-block:: shell

   lethe-fluid -R tgv-matrix-based.prm

The previous (if any) ``.vtu``, ``.pvtu``, and ``.pvd`` files, stored in the output path 
determined in the ``prm`` file, will be removed before the files corresponding
to the current simulation are saved.

.. hint::

   The ``-R`` tag can be used in the applications ``lethe-fluid``, ``lethe-fluid-block``,
   ``lethe-fluid-matrix-free``, ``lethe-fluid-nitsche``, ``lethe-fluid-particles``,
   ``lethe-fluid-sharp``, ``lethe-fluid-vans``, and ``lethe-particles``.