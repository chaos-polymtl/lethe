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
   make
   sudo make install

If Lethe's build was configured before TXR was installed, then the
`update-golden` target will not exist.
In this situation, the build has to be reconfigured, for example with:

.. code-block:: shell

   cmake -S. -Bbuild

assuming the current directory is Lethe's root directory, and the build
directory file:`build`.

-----
Usage
-----
Assuming the build directory is file:`build`, the following command
updates the golden files:

.. code-block:: shell
   
   cmake --build build --target update-golden
