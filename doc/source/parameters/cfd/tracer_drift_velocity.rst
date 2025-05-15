=====================
Tracer Drift Velocity
=====================

This subsection allows you to define a tracer drift velocity. This drift velocity is an additional velocity which is added to the fluid velocity when advecting the tracer. This enables user to model, with some strong hypothesis, the dynamics of a dispersed phase (such as bubbles or particles) within a fluid flow. 

.. code-block:: text

  subsection tracer drift velocity
    subsection drift velocity
      # Default values in 2D
      set Function expression = 0; 0
      # in 3D: set Function expression = 0; 0; 0
    end
  end
