Parameter tuning
-------------------

This subsection contains information regarding the tuning parameters with NOMAD. Enable tuning parameters requires the ``verbosity parameter`` in the subsection ``rpt parameters`` to be disabled by setting it as ``quiet`` otherwise it will interact with `NOMAD <https://www.gerad.ca/en/software/nomad/>`_ since it needs the cost function value only. So far there are 3 types of cost functions implemented, one from `Larachi et al. (1994) <https://www.sciencedirect.com/science/article/abs/pii/0168900294913439?via%3Dihub>`_, the L1 function, and the L2 function. To tune parameters, the cost function compares the calculated counts with the Monte Carlo technique and the measured counts that are provided in the ``.experimental`` file. The three parameters ``dead time``, ``activity`` and ``attenuation coefficient reactor`` seen in the `detector parameter subsection <./detector_parameters.html>`_ of the ``.prm`` file are obtained using NOMAD. The second example `Tuning Parameters with NOMAD <../../examples/rpt/tuning-parameters-with-nomad/tuning-parameters-with-nomad.html>`_ explains how we can obtain the values of these parameters using NOMAD.

.. code-block:: text

    # --------------------------------------------------
    # Tuning with NOMAD
    #---------------------------------------------------
    subsection parameter tuning
        set tuning                           = true
        set cost function type               = larachi
        set experimental data file           = real_counts.experimental
    end




- ``tuning``: Enable to tune parameters with NOMAD by showing the cost function in the terminal
    Options: ``true`` or ``false`` *(by default)*
- ``cost function type``: Type of cost function to evaluate
    Options: ``larachi`` *(by default)*, ``l1`` or ``l2``

    - Larachi cost function :
        .. math::

            f=\sum_{i=n}^{N}\left(\frac{C_i - M_i}{C_i + M_i}\right)^2

    - L1 cost function :
        .. math::

            f=\frac{1}{N}\sum_{i=n}^{N}\left|C_i - M_i\right|

    - L2 cost function :
        .. math::

            f=\frac{1}{N}\sum_{i=n}^{N}\left(C_i - M_i\right)^2

- ``experimental data file``: Filename of the text file with experimental/artificial counts
    Options: Any text file with ``.experimental`` extension with the
    **required header**: *experimental_counts*


References
~~~~~~~~~~~

[1] Larachi, F., Kennedy, G., & Chaouki, J. (1994). A γ-ray detection system for 3-D particle tracking in multiphase reactors. *Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment*. 338(2), 568-576. https://doi.org/10.1016/0168-9002(94)91343-9

