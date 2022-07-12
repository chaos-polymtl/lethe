Detector parameters
--------------------

This subsection contains the specific information of the detector. ``Detector position file`` defines the files that Lethe uses to read the detector positions. This file includes the position of the detector face's center and the position of a point inside the detector on its axis.

.. code-block:: text

    #---------------------------------------------------
    # Detector parameters
    #---------------------------------------------------
    subsection detector parameters
        set detector positions file          = positions.detector
        set radius                           = 0.0381
        set length                           = 0.0762
        set dead time                        = 1e-5
        set activity                         = 2e6
        set attenuation coefficient reactor  = 10
    end


- ``detector positions file``: Filename of the text file with positions of detector(s)
    Options: Any text file with ``.detector`` extension with the **required header**:
    *face_positions_x face_positions_y face_positions_z middle_positions_x middle_positions_y middle_positions_z*
- ``radius``: Radius of detector(s) (all detectors must have the same dimensions)
    Options: Any positive float *(default value: 1)*
- ``length``: Length of detector(s) (all detectors must have the same dimensions)
    Options: Any positive float *(default value: 1)*


The following parameters are variables in the gamma-ray Monte-Carlo model from `Beam et al. (1978) <https://www.sciencedirect.com/science/article/abs/pii/0029554X78900812?via%3Dihub>`_:

- ``activity``: Radioactive source activity of the tracer [Beq]
    Options: Any positive float *(default value: 1)*
- ``dead time``: Dead time of the detector per accepted pulse [s]
    Options : Any positive float *(default value: 1)*
- ``attenuation coefficient reactor``: Total linear attenuation of light coefficient of the medium and reactor wall []
    Options: Any positive float *(default value: 1)*

References
~~~~~~~~~~~

[1] Beam, G.B., Wielopolski, L., Gardner,  R.P., & Verghese, K. (1978). Monte Carlo calculation of efficiencies of right-circular cylindrical NaI detectors for arbitrarily located point sources. *Nuclear Instruments and Methods*. 154(3), 501-508. https://doi.org/10.1016/0029-554X(78)90081-2

