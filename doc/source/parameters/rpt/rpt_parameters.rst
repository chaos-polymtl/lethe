RPT parameters
-------------------
.. role:: raw-html(raw)
    :format: html

This subsection contains the general information required for the photon count calculation.

.. code-block:: text

    # --------------------------------------------------
    # RPT Monte Carlo technique
    #---------------------------------------------------
    subsection rpt parameters
        set particle positions file           = positions_diagonal.particle
        set verbosity                         = verbose
        set export counts                     = true
        set counts file                       = counts_diagonal.csv
        set monte carlo iteration             = 100000
        set random number seed                = 0
        set reactor height                    = 0.3
        set reactor radius                    = 0.1
        set peak-to-total ratio               = 0.4
        set sampling time                     = 1
        set gamma-rays emitted                = 2
        set attenuation coefficient detector  = 21.477
    end

- ``particle positions file``: Filename of the text file with a set of positions inside the reactor.
- ``verbosity``: Enable with ``verbose`` to show photon counts numbers results in realtime in terminal
    Options : ``quiet`` (by default) or ``verbose``
- ``export counts``: Enable to export photon counts numbers in a file
    Options : ``true`` or ``false`` (by default)
- ``counts file``: Filename of export counts file (.csv or .dat)
    Options : Any text file with .csv or .dat extension *(default value:* ``counts.csv`` *)*
- ``monte carlo iteration``: Defines the number of traced gamma-rays from each particle position to the detector within the defined solid angle
    Options : Any positive integer *(default value: 1)*
- ``random number seed``: Seed number for the random number generator, using a particular number allows to run the same series of numbers
    Options : Any positive integer *(default value: 0)*
- ``reactor height``: Height of the cylindrical reactor vessel [m]
    Options : Any positive float *(default value: 0.1)*
- ``reactor radius``: Radius of the cylindrical reactor vessel [m]
    Options : Any positive float *(default value: 0.1)*

:raw-html:`<br />`

The following parameters are variables in the gamma-ray Monte-Carlo model from `Beam <https://www.sciencedirect.com/science/article/abs/pii/0029554X78900812?via%3Dihub>`_:

- ``peak-to-total ratio``: The proportion of the events appearing in the full energy peak to the total number of events []
    Options : Any positive float *(default value: 1)*
- ``sampling time``: The amount of time for which the RPT hardware records the photon count at each position [s]
    Options : Any positive float *(default value: 1)*
- ``dead time``: Dead time of the detector per accepted pulse [s]
    Options : Any positive float *(default value: 1)*
- ``gamma-rays emitted``: Number of gamma-rays emitted by each disintegration []
    Options : Any positive float *(default value: 1)*
- ``attenuation coefficient detector``: Total linear attenuation of light coefficient of the detector []
    Options : Any positive float *(default value: 1)*

