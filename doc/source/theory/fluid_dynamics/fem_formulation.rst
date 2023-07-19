========================================
Origin of the Finite Element Formulation
========================================

This section describes the FEM formulation used within Lethe. Starting from the strong form of the equations, we obtain the weak-form. We then briefly discuss the challenges associated with solving the Navier-Stokes equations before we introduce the two approaches that are available in Lethe to solve them.


Starting from :doc:`The Incompressible Navier-Stokes <navier-stokes>` equations:

.. math::
    \partial_l u_l &= 0 

    \partial_t u_k + u_l \partial_l u_k &= -\frac{1}{\rho} \partial_k p^* + \nu \partial_l \partial_l u_k + f_k

We consider a domain :math:`\Omega` of contour :math:`\Gamma`. Without loss of generality, we assume Dirichlet boundary conditions or zero stress  conditions 
on :math:`\Gamma`. We multiply by two test functions :math:`q` and :math:`\mathbf{v}=v_k` for pressure and velocity respectively and integrate over the domain :math:`\Omega`. The resulting set of equation is:

.. math::

  &\int_{\Omega}  q  \partial_l u_l d\Omega =0 
  \\
  &\int_{\Omega}  v_k \left(\partial_t u_k+ u_l \partial_l u_k + \partial_k p - \nu \partial_l \partial_l u_k - f_k \right) d\Omega =0


Because we want the pressure to be in :math:`\mathcal{L}^2` and the velocity to be in :math:`\mathcal{H}^1`, we integrate by parts the viscous stress and the pressure gradient terms. We thus obtain the weak form:


.. math::

  &\int_{\Omega}  q  \partial_l u_l \mathrm{d}\Omega =0 
  \\
  &\int_{\Omega}  v_k \left(\partial_t u_k+ u_l \partial_l u_k - f_k \right) \mathrm{d}\Omega 
  \\
  &  - \int_{\Omega} \left( \partial_k \right) v_k p \mathrm{d}\Omega  
 + \nu \int_{\Omega} \left( \partial_l v_k \right) \left( \partial_l u_k  \right) \mathrm{d}\Omega  
 \\
  &  + \int_{\Gamma} \left( v_k \right) \left( \partial_l u_k  +\delta_{lk} p \right) n_l \mathrm{d}\Gamma
   =0

where :math:`\delta_{lk}` is the Kronecker delta and :math:`n_l` is the outward pointing normal vector to a surface. Since we assume Dirichlet boundary conditions or zero stress  conditions 
on :math:`\Gamma`, this term may be discarded. Thus, when no boundary condition is applied, the boundary condition applied is:

.. math::

    \int_{\Gamma} \left( v_k \right) \left( \partial_l u_k  +\delta_{lk} p \right) n_l \mathrm{d}\Gamma=0

which can be seen as an outlet boundary condition where the normal stress is zero. In essence, this can be used to approximately impose an outlet boundary condition with a zero average pressure.
This weak form is non-linear because of :math:`u_l \partial_l u_k` term. 


----------------------------------
Solving the Non-linear Problem
----------------------------------

To solve non-linear problem, Lethe uses the `Newton-Raphson method <https://en.wikipedia.org/wiki/Newton%27s_method>`_. This method proceeds by solving recurrently for the correction vector :math:`\mathbf{\delta x}` which is obtained by solving the following system:

.. math::

    \mathbf{\mathcal{J}} \mathbf{\delta x} = - \mathbf{\mathcal{R}}

For the incompressible Navier-Stokes equation, this leads to a saddle point problem of the form:

.. math::
    
  \mathbf{\mathcal{J}} &= \left[ \begin{matrix} 	A & B^T  \\[0.3em]	B & 0 \end{matrix} \right] \\
  \mathbf{\mathcal{R}} &=  \left[ \begin{matrix} \mathbf{\mathcal{R}}_v   \\[0.3em]		\mathbf{\mathcal{R}}_q  \end{matrix} \right]
  
  
The residual is:
  
.. math::

    \mathbf{\mathcal{R}} &=    \int_{\Omega}  q  \partial_l u_l \mathrm{d}\Omega 
    +   \int_{\Omega}  v_k \left(\partial_t u_k+ u_l \partial_l u_k - f_k \right) \mathrm{d}\Omega \\
    &  - \int_{\Omega} p\left( \partial_k   v_k\right) \mathrm{d}\Omega  
    + \nu \int_{\Omega} \left( \partial_l v_k \right) \left( \partial_l u_k  \right) \mathrm{d}\Omega  
  
  
We recall that in FEM, the pressure :math:`p` and the velocity :math:`u_k` are obtained from the discrete nodal values from the following:

.. math::
   p &= \sum_j p_j \psi_{j}   \\
   u_k &= \sum_j u_{k,j} \phi_{k,j}   \\

where  :math:`\psi_j` is the :math:`j` interpolation function for pressure and  :math:`\phi_{k,j}` is the :math:`j` interpolation function for the :math:`k` component of the velocity.  :math:`p_j` is the nodal value of pressure and :math:`u_{k,j}` is the nodal value of the  :math:`k` component of the velocity.


The Jacobian is:
  
.. math::

    \mathbf{\mathcal{J}} &=    \int_{\Omega}  q  \partial_l \phi_{l,j} \mathrm{d}\Omega 
    +   \int_{\Omega}  v_k \left(\partial_t \phi_{k,j}+ \phi_{l,j} \partial_l u_k + u_l \partial_l \phi_{k,j}  \right) \mathrm{d}\Omega \\
    &  - \int_{\Omega} \psi_j  \left( \partial_k   v_k \right)\mathrm{d}\Omega  
    + \nu \int_{\Omega} \left( \partial_l v_k \right) \left( \partial_l \phi_{k,j}  \right) \mathrm{d}\Omega  


