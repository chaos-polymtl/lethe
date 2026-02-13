// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @mainpage
 *
 * An outline of the main classes in Lethe and how they interact is given by the following
 * clickable graph:
 *
 * @dot
 *    digraph main_classes {
      graph [bgcolor="transparent", align=true, ranksep=1, nodesep=0.5];
      node [fontname="FreeSans",fontsize=15,
        shape=record,height=0.2,width=0.4,
        color="royalblue", fillcolor="white", style="filled"];
      edge [color="royalblue", weight=10];
      rankdir="LR";
      size = "16,10";


      //---------------------------------
      // Particles Ray Tracing
      //---------------------------------

      particles_ray_tracing_solver [label=<<B>RayTracingSolver</B> <br/>(lethe-particles-ray-tracing)>, href="https://chaos-polymtl.github.io/lethe/doxygen/classRayTracingSolver.html", tooltip="RayTracingSolver", fillcolor="#8ba5d4"];


      //---------------------------------
      // DEM
      //---------------------------------

      dem_solver [label=<<B>DEMSolver</B> <br/>(lethe-particles)>, href="https://chaos-polymtl.github.io/lethe/doxygen/classDEMSolver.html", tooltip="DEMSolver", fillcolor="#8ba5d4"];


      //---------------------------------
      // Physics Solver
      //---------------------------------

      physics_solver [label="PhysicsSolver", href="https://chaos-polymtl.github.io/lethe/doxygen/classPhysicsSolver.html"];

      navier_stokes_base [label="NavierStokesBase",href="https://chaos-polymtl.github.io/lethe/doxygen/classNavierStokesBase.html"];

      auxiliary_physics [label="AuxiliaryPhysics",href="https://chaos-polymtl.github.io/lethe/doxygen/classAuxiliaryPhysics.html"];

      physics_solver:e -> navier_stokes_base:w [dir=back];
      physics_solver:e -> auxiliary_physics:w [dir=back];

      navier_stokes_base_1 [label=<<B>FluidDynamicsMatrixBased</B> <br/>(lethe-fluid)>,href="https://chaos-polymtl.github.io/lethe/doxygen/classFluidDynamicsMatrixBased.html", tooltip="FluidDynamicsMatrixBased", fillcolor="#8ba5d4"];
      navier_stokes_base_2 [label=<<B>FluidDynamicsBlock</B> <br/> (lethe-fluid-block)>,href="https://chaos-polymtl.github.io/lethe/doxygen/classFluidDynamicsBlock.html", tooltip="FluidDynamicsBlock", fillcolor="#8ba5d4"];
      navier_stokes_base_3 [label=<<B>FluidDynamicsMatrixFree</B> <br/> (lethe-fluid-matrix-free)>,href="https://chaos-polymtl.github.io/lethe/doxygen/classFluidDynamicsMatrixFree.html", tooltip="FluidDynamicsMatrixFree", fillcolor="#8ba5d4"];

      navier_stokes_base:e -> navier_stokes_base_1:w [dir=back];
      navier_stokes_base:e -> navier_stokes_base_2:w [dir=back];
      navier_stokes_base:e -> navier_stokes_base_3:w [dir=back];

      auxiliary_physics_1 [label="VolumeOfFluid",href="https://chaos-polymtl.github.io/lethe/doxygen/classVolumeOfFluid.html"];
      auxiliary_physics_2 [label="CahnHilliard",href="https://chaos-polymtl.github.io/lethe/doxygen/classCahnHilliard.html"];
      auxiliary_physics_3 [label="HeatTransfer",href="https://chaos-polymtl.github.io/lethe/doxygen/classHeatTransfer.html"];
      auxiliary_physics_4 [label="Tracer",href="https://chaos-polymtl.github.io/lethe/doxygen/classTracer.html"];
      auxiliary_physics_5 [label="TimeHarmonicMaxwell",href="https://chaos-polymtl.github.io/lethe/doxygen/classTimeHarmonicMaxwell.html"];

      auxiliary_physics:e -> auxiliary_physics_1:w [dir=back];
      auxiliary_physics:e -> auxiliary_physics_2:w [dir=back];
      auxiliary_physics:e -> auxiliary_physics_3:w [dir=back];
      auxiliary_physics:e -> auxiliary_physics_4:w [dir=back];
      auxiliary_physics:e -> auxiliary_physics_5:w [dir=back];

      navier_stokes_base_1_1 [label=<<B>FluidDynamicsVANS</B> <br/>(lethe-fluid-vans)>,href="https://chaos-polymtl.github.io/lethe/doxygen/classFluidDynamicsVANS.html", tooltip="FluidDynamicsVANS", fillcolor="#8ba5d4"];
      navier_stokes_base_1_2 [label=<<B>FluidDynamicsSharp</B> <br/>(lethe-fluid-sharp)>,href="https://chaos-polymtl.github.io/lethe/doxygen/classFluidDynamicsSharp.html", tooltip="FluidDynamicsSharp", fillcolor="#8ba5d4"];
      navier_stokes_base_1_3 [label=<<B>FluidDynamicsNitsche</B> <br/>(lethe-fluid-nitsche)>,href="https://chaos-polymtl.github.io/lethe/doxygen/classFluidDynamicsNitsche.html", tooltip="FluidDynamicsNitsche", fillcolor="#8ba5d4"];

      navier_stokes_base_1:e -> navier_stokes_base_1_1:w [dir=back];
      navier_stokes_base_1:e -> navier_stokes_base_1_2:w [dir=back];
      navier_stokes_base_1:e -> navier_stokes_base_1_3:w [dir=back];

      navier_stokes_base_1_1_1 [label=<<B>CFDDEMSolver</B> <br/>(lethe-fluid-particles)>,href="https://chaos-polymtl.github.io/lethe/doxygen/classCFDDEMSolver.html", tooltip="CFDDEMSolver", fillcolor="#8ba5d4"];

      navier_stokes_base_1_1:e -> navier_stokes_base_1_1_1:w [dir=back];

      navier_stokes_base_3_1 [label=<<B>FluidDynamicsVANSMatrixFree</B> <br/>(lethe-fluid-vans-matrix-free)>,href="https://chaos-polymtl.github.io/lethe/doxygen/classFluidDynamicsVANSMatrixFree.html", tooltip="FluidDynamicsVANSMatrixFree", fillcolor="#8ba5d4"];

      navier_stokes_base_3:e -> navier_stokes_base_3_1:w [dir=back];

      navier_stokes_base_3_1_1 [label=<<B>CFDDEMMatrixFree</B> <br/>(lethe-fluid-particles-matrix-free)>,href="https://chaos-polymtl.github.io/lethe/doxygen/classCFDDEMMatrixFree.html", tooltip="CFDDEMMatrixFree", fillcolor="#8ba5d4"];

      navier_stokes_base_3_1:e -> navier_stokes_base_3_1_1:w [dir=back];
    }
 * @enddot
 */
