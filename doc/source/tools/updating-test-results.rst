=====================
Updating Test Results
=====================
Lethe's test suite compares individual test results against reference
results, stored in so-called golden files, which are known, to some
extent, to be good.
However, these golden files can turn out to be wrong, or can change
slightly, either in format or numerical precision.

Updating the golden files manually is tedious, and so the
`update-golden` build target exists to automate this process.
The `update-golden` target launches the
file:`contrib/utilities/update-golden.tl` script with appropriate
arguments.

The file:`contrib/utilities/update-golden.tl` script's only dependency
is TXR_.

.. _TXR: https://www.nongnu.org/txr/

-----
Setup
-----
TXR can be installed to file:`/usr/local/bin` as follows:

.. code-block:: shell

   curl 'https://www.kylheku.com/cgit/txr/snapshot/txr-285.tar.gz' | tar -xzf -
   cd txr-285
   ./configure
   make
   sudo make install

Once that TXR is installed, the `update-golden.tl` script can be run in two ways which will be discussed in :ref:`Usage`.
The first way requires adding the script as a target in your CMake build directory.
If Lethe's build was configured before TXR was installed, then the `update-golden` target will not exist.
In this situation, the build has to be reconfigured, for example with:

.. code-block:: shell

   cmake ../lethe

assuming the current directory is your `build` directory and that you have followed the installation instruction.

.. _Usage:
-----
Usage
-----
If the script was added as a target in your CMake build directory and assuming that you are in your `build` directory,
the following command updates the golden files:

.. code-block:: shell
   
   cmake --build . --target update-golden

The second option to run the `update-golden.tl` script is to run it manually using the following command:

.. code-block:: shell
    txr ../lethe/contrib/utilities/update-golden.tl  -v . ../lethe

assuming that TXR is installed in `usr/bin` on your machine, that you are in your `build` directory and that you have
followed the installation instructions.