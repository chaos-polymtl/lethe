========================================
Massive Data Visualization with Paraview
========================================

Lethe is developed with special care for scalability; it is intended to be used for solving extremely large problems by using massively parallel computing. As such, the output files that are obtained can sometimes be so heavy that the post-processing requires special considerations. This section covers key points to take into account in such cases.


------------------------------
Client-server visualization
------------------------------

The first limitation that is encountered when FEM problems contain :math:`\sim 15`M cells (when using a 32GB-RAM machine) is lack of memory. Even though the VTU files weigh 2GB, they are decompressed by Paraview when opened.


Clusters such as those of the Digital Research Alliance of Canada allow to launch interactive jobs in which a Paraview-server can be used. Detailed instructions on how to setup and use such a server are available on this `page <https://docs.alliancecan.ca/wiki/ParaView>`_. However, here is an overview of the required steps to use Client-server visualization:



1. Open 3 terminals: #1 will be used to run the server, #2 will be used to link the local-remote ports and #3 will be used to run the client.

2. Ensure that the files to visualize are on the cluster. ``scp`` can be used to move them there.
3. In #1:
    1. Log in to the cluster.
    2. Request an interactive job, specifying the number of tasks per node, the number of nodes, the memory per node, the required time and the account : ``salloc --ntasks-per-node=5 --mem=200G --time=2:00:00 --account=rrg-blaisbru --nodes=2``. The selection of the parameters is covered in the next subsection.
    3. Once the job is running, load the Paraview module : ``module load gcc/9.3.0 paraview-offscreen/5.11.0``. The Paraview version must be the same as the client. The cluster used in this example is Narval.
    4. Start the server : ``srun pvserver --force-offscreen-rendering``. ``srun`` is used to run the server with all the available cores. Make a note of the assigned node.
4. In #2: start the port tunneling by executing ``ssh user@narval.computecanada.ca -L 11111:nc10133:11111``. Use the assigned node (in this case ``nc10133``) and the port declared by the server (here ``11111``).

5. In #3:
    1. Launch Paraview: ``paraview``.
    2. Add the server in File -> Connect -> Add Server.
    3. Input the proper port, and assign a name to the connection.
    4. Then, click Connect and open files as in local mode: the shown files should be the ones on the cluster. Notice that the available memory shown in the lower right is now much higher than the one that is shown when Paraview is used locally.



------------------------------------
Splitting the output files in groups
------------------------------------

The ``group files`` parameter in the ``simulation control`` `subsection <../../parameters/cfd/simulation_control.html>`_ is used to split the obtained VTU files in multiple parts. Not only does this make them easier to manage and copy, but it also helps for post-processing.

Paraview, when used in parallel as it is described in this section, aims to balance the load on the available cores. This implies that if an output contains only one VTU, only one core will be able to open it. If this core belongs to a node that doesn't have enough memory, Paraview will not be able to open the file.

For this reason, we have to take into account post-processing in advance and select the proper number of parts. The number of parts should:

* Allow the individual parts to have a size equal or lower than 4GB, for performance.
* Be high enough so that Paraview can benefit from parallelization and share the load on its available cores.

As an example, a 300M cells problem can be split in about 10 parts of 30M cells each. These 10 parts need about 400GB of memory to be opened in Paraview. When using a cluster where nodes have 200GB of memory, we need 2 nodes to have access to enough memory. To ensure that each node has to open only 5 parts, we request 5 cores per node. We then have to use the line ``salloc –ntasks-per-node=5 –mem=200G –time=2:00:00 –account=rrg-blaisbru –nodes=2`` to request the appropriate resources. A higher number of parts can be used so that post-processing can be done using the actual available cores on each node. More information on the resources available on the Alliance clusters can be found `here <https://docs.alliancecan.ca/wiki/>`_.

