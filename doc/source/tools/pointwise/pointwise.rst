=====================================
Introduction on How to Use Pointwise
=====================================

--------------------------
Installation
--------------------------

Fidelity Pointwise is a  mesh generation proprietary software. To download it, please refer to the
official instructions provided by the Fidelity Pointwise documentation or to the designated person in your group.

For Linux users, start by extracting the .tgz file that holds the software. To finalize the installation, open a terminal and move to the directory where you extracted the file. You will need to make it executable by typing the following:

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

dont forget to erase the copy button on sphinx