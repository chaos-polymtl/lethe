=============================================================
Finite Element Method for Time-Harmonic Maxwell's Equations
=============================================================

This section describes the FEM formulation of the time-harmonic Maxwell's equations. For further details on the subject, the reader is referred to the following references: `Larson <https://books.google.ca/books/about/The_Finite_Element_Method_Theory_Impleme.html?id=Vek_AAAAQBAJ&redir_esc=y>`_ and `Jian-Ming <https://www.wiley.com/en-br/The+Finite+Element+Method+in+Electromagnetics%2C+3rd+Edition-p-9781118571361>`_ . 

Starting from the strong form of the electromagnetic field showed in the :doc:`Time Harmonic Maxwell's Equations <time_harmonic>` section:

.. math::
    \nabla \times \left( \mu_{em}^{-1} \nabla \times \mathbf{\hat{E}} \right) -\omega^2 \varepsilon_{em_{eff}} \mathbf{\hat{E}} = -i \omega \mathbf{\hat{J}}_{ext}, \\

we consider a domain :math:`\Omega` with boundary :math:`\Gamma`.

.. note::
    In general, it is not necessary to solve both the magnetic and the electric field wave equation since one results from the other. One can choose the field to solve that suits the best the problem under consideration.

The choice of boundary conditions can vary greatly depending on the problem to be solved (e.g., Sommerfeld condition, absorbing boundary condition, impedance condition, etc.), but the two most common types are:

- Dirichlet boundary conditions: :math:`\mathbf{n} \times \mathbf{\hat{E}} = \mathbf{\hat{E}}_{Dirichlet}`;
- Neumann boundary conditions: :math:`\mu_{em}^{-1}(\nabla \times \mathbf{\hat{E}}) \times \mathbf{n} = -i \omega \mathbf{J}_{ext_{Neumann}}`.

The :math:`\mathbf{n}` represent the unit normal vector. For simplicity, in the current derivation of the weak form, a perfect electric conductor (PEC) is considered, which implies that :math:`\mathbf{E}_{Dirichlet} = 0`. Now, multiplying the strong form by a complex test function :math:`\mathbf{v}` that satisfies :math:`\mathbf{v} \times \mathbf{n}=0` and integrating over the domain :math:`\Omega`:

.. math::
    \begin{align*}
     \int_{\Omega} -i \omega \mathbf{J}_{ext} \cdot \mathbf{v^*} \mathrm{d}\Omega = &
     \int_{\Omega}  \mu_{em}^{-1} (\nabla \times \mathbf{\hat{E}}) \cdot (\nabla \times \mathbf{v^*}) \mathrm{d}\Omega - \int_{\Omega} \omega^2 \varepsilon_{em_{eff}} \mathbf{\hat{E}} \cdot \mathbf{v^*} \mathrm{d}\Omega \\ & + \int_{\Gamma} \mu_{em}^{-1} (\nabla \times \mathbf{\hat{E}}) \cdot (\mathbf{v^*} \times \mathbf{n}) \mathrm{d}\Gamma .
    \end{align*} \\


Note that the terms above involve complex values, so the usual :math:`L^2` inner product is replaced by :math:`\int_{\Omega} \mathbf{u}\mathbf{v^*} \mathrm{d}\Omega`, where :math:`\mathbf{v^*}` is the complex conjugate of :math:`\mathbf{v}`. Finally, for the above integral to be well defined, :math:`\mathbf{E}` and :math:`\mathbf{v}` must belong to the Sobolev space :math:`H(curl, \Omega)` defined by:

.. math::
    H(curl, \Omega) = \{ \mathbf{v} \in [L^2(\Omega)]^3 : \nabla \times \mathbf{v} \in [L^2(\Omega)]^3 \},

and the Dirichlet boundary condition imposed the homogeneous curl space:

.. math::
    H_0(curl, \Omega) = \{ \mathbf{v} \in H(curl, \Omega) : \mathbf{n} \times \mathbf{v}|_{\Gamma} = 0 \}.
    
Thus the boundary term vanishes (:math:`\int_{\Gamma} \mu_{em}^{-1} (\nabla \times \mathbf{\hat{E}}) \cdot (\mathbf{v^*} \times \mathbf{n}) \mathrm{d}\Gamma = 0`) and the weak form of the time-harmonic Maxwell's equations is:

.. math::
    B(\mathbf{E}, \mathbf{v}) = &\int_{\Omega}  \mu_{em}^{-1} (\nabla \times \mathbf{\hat{E}}) \cdot (\nabla \times \mathbf{v^*}) \mathrm{d}\Omega + \int_{\Omega} \omega^2 \varepsilon_{em_{eff}} \mathbf{\hat{E}} \cdot \mathbf{v^*} \mathrm{d}\Omega  \\
    L(\mathbf{v}) = &\int_{\Omega} -i \omega \mathbf{J}_{ext} \cdot \mathbf{v^*} \mathrm{d}\Omega .

Formally, :math:`\mathbf{E}` should also satisfy Gauss's law (:math:`\nabla \cdot \mathbf{D} = \rho_f`), but it is implicitly taken into account by the electromagnetic wave equation and holds in the weak form presented above.