Notation
=========

We first summarize the notation that we use in the Lethe theory guide.
This is intended as a reference for the rest of the documentation.

Conventions
-----------
In this documentation, we use both traditional boldface vector/tensor
notation and Einstein (index) notation. All expressions are written in
three-dimensional Cartesian coordinates with an orthonormal basis :math:`\{\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3\}`.

The following conventions are used:

* Latin indices :math:`i,j,k,l` run from 1 to 3.
* Repeated indices in a single term imply summation over that index
  (Einstein summation convention).
* An index that appears exactly once in each term of an equation is a
  *free* index. Free indices must match on both sides of an equation.
* An index that appears twice in a single term is a *dummy* (summed)
  index and is implicitly summed over.


Simple contraction
~~~~~~~~~~~~~~~~~~

The simple contraction of two vectors :math:`\mathbf{a}` and
:math:`\mathbf{b}` with components :math:`a_i` and :math:`b_i` is the
scalar:

.. math::

      a_i b_i \equiv  \mathbf{a} \cdot \mathbf{b}.

A simple contraction sums over a single repeated index and reduces the rank of the quantity by one.

Double contraction
~~~~~~~~~~~~~~~~~~

The double contraction of two second-order tensors :math:`\mathbf{A}`
and :math:`\mathbf{B}` with components :math:`A_{ij}` and
:math:`B_{ij}` is the scalar:

.. math::

   A_{ij} B_{ij} \equiv \mathbf{A} : \mathbf{B}.

A double contraction sums over two repeated indices and reduces the
rank of the quantity by two.

Dyadic (tensor) product
~~~~~~~~~~~~~~~~~~~~~~~

The dyadic product (or tensor product) of two vectors :math:`\mathbf{a}` and :math:`\mathbf{b}` with components :math:`a_i` and :math:`b_j` is the second-order tensor:

.. math::

   a_i b_j \equiv \mathbf{a} \otimes \mathbf{b} \equiv \mathbf{a} \mathbf{b}.

More generally, the tensor product combines two tensors by multiplying their components without summation over repeated indices, thereby increasing the total rank by the sum of the individual ranks.


Weak formulations
------------------
When writing weak form, we use the standard inner product notation from functional analysis when expressing products between trial and test spaces. This helps distinguish between algebraic operations and the weak formulation of the problem. To avoid ambiguity, we retain a subscript to indicate the domain over which each integration is performed. It can be either the domain :math:`\Omega`, its boundary :math:`\partial \Omega`, or the sum of the elements :math:`\Omega_e`. The corresponding notation convention for the integrals in vector notation is:

.. math::
   :nowrap:

   \begin{aligned}
   (a, b)_{\Omega} &= \int a b \,\mathrm{d}\Omega \\
   (\mathbf{a}, \mathbf{b})_{\Omega} &= \int \mathbf{a} \cdot \mathbf{b} \,\mathrm{d}\Omega \\
   (\nabla \mathbf{a}, \nabla \mathbf{b})_{\Omega} &= \int \nabla \mathbf{a} : \nabla \mathbf{b} \,\mathrm{d}\Omega \\
   (a, b)_{\Omega_e} &= \sum_e \int a b \,\mathrm{d}\Omega_e
   \end{aligned}

And in Einstein notation:


.. math::
   :nowrap:

   \begin{aligned}
   (a, b)_{\Omega} &= \int a \, b \,\mathrm{d}\Omega \\
   (a_i, b_i)_{\Omega} &= \int a_i \, b_i \,\mathrm{d}\Omega \\
   (\partial_j a_i, \partial_j b_i)_{\Omega} &= \int \partial_j a_i \, \partial_j b_i \,\mathrm{d}\Omega \\
   (a, b)_{\Omega_e} &= \sum_e \int a \, b \,\mathrm{d}\Omega_e
   \end{aligned}
