=====================================
Remote Data Visualization on ParaView
=====================================

Lethe is developed with special care for scalability; it is intended to be used for solving extremely large problems by using massively parallel computing. When simulations are carried out on computing clusters, it can be impossible or unpractical to download all results. As such, the output files that are obtained can sometimes be so heavy that the post-processing requires special considerations. This section covers key points to take into account in such cases.


------------------------------
Client-server visualization
------------------------------

The first limitation that is encountered when FEM problems contain :math:`\sim 15` M cells is lack of memory (when using a 32GB-RAM machine). Even though the VTU files weigh 2GB, they are decompressed by ParaView when opened.

Clusters such as those of the Digital Research Alliance of Canada allow to launch interactive jobs in which a ParaView-server can be used. Detailed instructions on how to set up and use such a server are available on this `page from the Alliance <https://docs.alliancecan.ca/wiki/ParaView>`_ and this `page from ParaView <https://www.paraview.org/Wiki/Setting_up_a_ParaView_Server>`_. However, here is an overview of the required steps to use Client-server visualization, in an example using the Narval cluster:

1. Open 3 terminals: #1 will be used to run the server, #2 will be used to link the local-remote ports and #3 will be used to run the client.

2. In #1:

    1. Log in to the cluster where the data is located.
    2. Request an interactive job, specifying the number of tasks per node, the number of nodes, the memory per node, the required time and the account: ``salloc --ntasks-per-node=5 --mem=200G --time=2:00:00 --account=def-someprof --nodes=2``. The selection of the parameters is covered in the next subsection.

    3. Once the job is running, load the ParaView module : ``module load gcc/9.3.0 paraview-offscreen/5.11.0``. The ParaView version must be the same as the client.

    4. Start the server: ``srun pvserver --force-offscreen-rendering``. ``srun`` is used to run the server with all the available cores. Make a note of the assigned node.

3. In #2: start the port tunneling by executing ``ssh user@narval.computecanada.ca -L 11111:nc10133:11111``. Use the assigned node (in this case ``nc10133``) and the port declared by the server (here ``11111``).

4. In #3:
    1. Launch ParaView: ``paraview``.

    2. Add the server in File -> Connect -> Add Server.
    3. Input the proper port, and assign a name to the connection.
    4. Then, click Connect and open files as in local mode: the shown files should be the ones on the cluster. Notice that the available memory shown in the lower right is now much higher than the one that is shown when ParaView is used locally.


------------------------------------
Splitting the output files in groups
------------------------------------

The ``group files`` parameter in the ``simulation control`` `subsection <../../parameters/cfd/simulation_control.html>`_ is used to split the obtained VTU files into multiple parts. Not only does this make them easier to manage and copy, but it also helps with post-processing.


ParaView, when used in parallel, aims to balance the load on the available cores. This implies that if an output contains only one VTU, only one core will be able to open it. If this core belongs to a node that does not have enough memory, ParaView will not be able to open the file.



For this reason, we have to take into account post-processing in advance and select the proper number of parts. The number of parts should:

* Allow the individual parts to have a size equal to or lower than 4 GB, for performance.

* Be high enough so that ParaView can benefit from parallelization and share the load on its available cores.


As an example, a 300 M cells problem can be split into about 10 parts of 30M cells each. These 10 parts need about 400 GB of memory to be opened in ParaView. When using a cluster where nodes have 200 GB of memory, we need 2 nodes to have access to enough memory. To ensure that each node has to open only 5 parts, we request 5 cores per node. We then have to use the line ``salloc --ntasks-per-node=5 --mem=200G --time=2:00:00 --account=def-someprof --nodes=2`` to request the appropriate resources. A higher number of parts can be used so that post-processing can be done using the available cores on each node. More information on the resources available on the Alliance clusters can be found `here <https://docs.alliancecan.ca/wiki/>`_.


