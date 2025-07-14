================
Animation Videos
================

Here we cover the steps for saving a video using Paraview visualization tool; if you are using a different tool, keep in mind that the video formatting rules still apply.

Making the video
^^^^^^^^^^^^^^^^
``Save Animation`` on Paraview saves the animation in images. They can be combined in a video using ``ffmpeg``, with an appropriate framerate and from inside the images folder:

.. code-block:: text
 :class: copy-button
 
  ffmpeg -framerate 20 -i <imageprefix>.%04d.png -c:v libx264 -pix_fmt yuv420p ../out.mp4 

Formatting
^^^^^^^^^^
- Check if the legend contains the corresponding units of the represented field (e.g., [length/time] unit when plotting velocity).
- Make sure that the video contains the Lethe logo within it (add on Paraview with Sources/Annotations/Logo). The file can be found in the `repository <https://github.com/chaos-polymtl/lethe/tree/master/logo>`_.

Video Upload
^^^^^^^^^^^^
Animation videos are uploaded on the `Lethe Youtube page <https://www.youtube.com/@lethecfd6431>`_.

When adjusting the video configurations:

- Add a descriptive title to the video. This usually includes the type of solver (e.g., CFD, DEM, CFD-DEM) and the type of problem solved (e.g. `Transient flow around a cylinder <https://www.youtube.com/watch?v=NbN2kBdakH4>`_).
- When starting the upload process, a template for the video description will appear, which relates the video to its documentation page. 
- Most of the animations uploaded are part of an example listed in the Lethe documentation. If that is your case, remember to add the correspondent website link in the video description.
- Add the video to its correspondent playlist: CFD, CFD-DEM, or DEM.
- Ensure that the video audience is set as "No, it's not Made for Kids".

.. tip::

    In the Video Elements section, you can add cards to promote related content during the video. A good practice is to add a card to reference a related playlist or video to the current upload, so that people watching it know there's more to explore!

- In the Visibility section, choose the "Public" configuration for who's allowed to see the video.
- After checking if all your input information is correct, just click on "Publish" and your animation is ready to be accessed.
