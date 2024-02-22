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

2. Make sure the entity type is set to connector (the green curved line). Always remember that a majority of the meshes can be assembled by only using connectors. The database entity type (the pink squigly rectangle) is your worst ennemy if you are not using imported CAD geometries. 

3. In the XYZ box, enter the coordinates ``0 0 0`` of the first point of the 2 Point Curve and press enter. Reselect the XYZ box. Enter the coordinates ``0 1 0`` of the second point of the 2 Point Curve and press enter again. This will create the left wall of our geometry.


At this point, you should have a vertical connector that starts from the origin (0,0,0) to the point (0,1,0) as is displayed in the image below.


.. image:: images/two_point_curve.png
    :align: center

4. Now that your first connector is created lets finish our geometry. Add the five following connectors. 

- (0, 1, 0) to (1, 1, 0)
- (1, 1, 0) to (1, 0.5, 0)
- (1, 0.5, 0) to (0.5, 0.5, 0)
- (0.5, 0.5, 0) to (0.5, 0, 0)
- (0.5, 0, 0) to (0, 0, 0)
- (0.5, 0.5, 0) to (0, 0.5, 0)
- (0.5, 0.5, 0) to (0.5, 1, 0)

.. tip::

    In a hurry? You can enter a coordinate by clicking wherever you want on the geometry. A small target lets you snap the point directly at the extremity of a connector. 

The addition of these connectors delimits the edges of the geometry. By default, these connectors have no dimension. Remember that a structured mesh will be created by selecting a closed quadrilateral of connectors. The only criteria is the sum of the size of one side of the quad needs to be equal to the size of its opposing side. Therefore, to make a 10x12 mesh, two opposing connectors will need to have a size of 12 and two others a size 10.

Sadly, our geometries are rarely nice convex quadrilaterals. It is therefore difficult to produce nice meshes. The goal is therefore to separate the geometry in sections that will facilitate the meshing by creating sections nice trapezoidal or rectangular sections. Exemples of meshed special geometries are presented at the end to help you figure out how to discritize some weird shapes.  

5. To size the connectors, select them and enter 20 in the dimension box located at the top of your window beside the green hashtag.

6. It is now time to mesh. In the create option. Select assemble special and domain. A new window with a square should appear on the left side of your screen. Select 





This criteria makes it easier for us to separate our geometries in structures that ressemble quads



3. To make the top wall of our geometry, click on the top point of the first connector we just created. This will select.



dont forget to erase the copy button on sphinx