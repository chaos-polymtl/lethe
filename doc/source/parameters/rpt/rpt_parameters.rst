==============
RPT Parameters
==============

This subsection contains the general information required for the photon count calculation using a Monte Carlo technique

.. code-block:: text

  subsection rpt parameters
    set particle positions file          = positions_diagonal.particle
    set verbosity                        = verbose
    set export counts                    = true
    set counts file                      = counts_diagonal.csv
    set monte carlo iteration            = 100000
    set random number seed               = 0
    set reactor height                   = 0.3
    set reactor radius                   = 0.1
    set peak-to-total ratio              = 0.4
    set sampling time                    = 1
    set gamma-rays emitted               = 2
    set attenuation coefficient detector = 21.477
  end


- ``particle positions file``: Filename of the file with a set of particle positions inside the reactor.
- ``verbosity``: Enable with ``verbose`` to show photon counts numbers results in realtime in terminal
    Options: ``quiet`` *(by default)* or ``verbose``
- ``export counts``: Enable to export photon counts numbers in a file
    Options: ``true`` or ``false`` *(by default)*
- ``counts file``: Filename of export counts file (.csv or .dat)
    Options: Any text file with .csv or .dat extension *(default value:* ``counts.csv`` *)*
- ``monte carlo iteration``: Defines the number of traced gamma-rays from each particle position to the detector within the defined solid angle
    Options: Any positive integer *(default value: 1)*
- ``random number seed``: Seed number for the random number generator, using a particular number allows to run the same series of numbers
    Options: Any positive integer *(default value: std::time(nullptr))*
- ``reactor height``: Height of the cylindrical reactor vessel [m]
    Options: Any positive float *(default value: 0.1)*
- ``reactor radius``: Radius of the cylindrical reactor vessel [m]
    Options: Any positive float *(default value: 0.1)*


The following parameters are variables in the gamma-ray Monte-Carlo model from Beam *et al.* (1978) `[1] <https://doi.org/10.1016/0029-554X(78)90081-2>`_:

- ``peak-to-total ratio``: The proportion of the events appearing in the full energy peak to the total number of events []
    Options: Any positive float *(default value: 1)*
- ``sampling time``: The amount of time for which the RPT hardware records the photon count at each position [s]
    Options: Any positive float *(default value: 1)*
- ``gamma-rays emitted``: Number of gamma-rays emitted by each disintegration []
    Options: Any positive float *(default value: 1)*
- ``attenuation coefficient detector``: Total linear attenuation of light coefficient of the detector []
    Options: Any positive float *(default value: 1)*


References
~~~~~~~~~~~

`[1] <https://doi.org/10.1016/0029-554X(78)90081-2>`_ G. B. Beam, L. Wielopolski, R. P. Gardner, and K. Verghese, “Monte Carlo calculation of efficiencies of right-circular cylindrical NaI detectors for arbitrarily located point sources,” *Nucl. Instrum. Methods*, vol. 154, no. 3, pp. 501–508, Sep. 1978, doi: 10.1016/0029-554X(78)90081-2.

