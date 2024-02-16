=====================================
Introduction on How to Use Pointwise
=====================================

--------------------------
Installation
--------------------------

Fidelity Pointwise is a  mesh generation proprietary software. To download it, please refer to the
official instructions provided by the Fidelity Pointwise documentation or to the designated person in your group.

For Linux users, start by extracting the .tgz file that holds the software. Open a terminal and move to the directory where you extracted the file. Make the file executable by typing the following:

.. code-block:: text
    
    chmod +x your_file.sh

Execute the installation file with:

.. code-block:: text

    ./your_file.sh

At the end of the installation, add the the following to your ``.bashrc`` system file:

.. code-block:: text

    #Pointwise license
    export CDS_LIC_FILE=license_server_ip_adress
    export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6

    #Launch pointwise
    alias pointwise=~/Fidelity/Pointwise/Pointwise2023.1.1/pointwise

Make sure to modify the ``license_server_ip_adress`` with the actual server IP adress so it can be accessed by the software. After these modifications, you can launch pointwise from the terminal by entering ``pointwise``.

----------------------------
Make your first 2D Geometry
----------------------------

This section will give you a good idea of the tools that can be used to create a 2D geometry of the elbow of a pipe.


1. Select the 2 Point Curve in the shortcut bar at the top of the pointwise.

2. Make sure the entity type is set to connector (the green curved line). Always remember that database entities (the pink rectangle) are your worst enemy if you are not using imported CAD geometries. A majority of the meshes can be made by only using connectors . 

3. In the XYZ box, enter the coordinates ``0 0 0`` of the first point of the 2 Point Curve and press enter. Reselect the XYZ box. Enter the coordinates ``0 1 0`` of the second point of the 2 Point Curve and press enter again. This will create the left wall of our geometry.

.. image:: images/two_point_curve.png
    :align: center



.. tip::

    In a hurry? You can enter a coordinate by clicking wherever you want on the geometry.

3. To make the top wall of our geometry, click on the top point of the first connector we just created. This will select.

dont forget to erase the copy button on sphinx