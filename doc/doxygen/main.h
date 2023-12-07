/**
 * @mainpage
 *
 * An outline of the main classes in Lethe and how they interact is given by the following
 * clickable graph:
 *
 * @dot
 *    digraph main_classes {
      graph [bgcolor="transparent", align=true, ranksep=1.5];
      node [fontname="FreeSans",fontsize=15,
        shape=record,height=0.2,width=0.4,
        color="royalblue", fillcolor="white", style="filled"];
      edge [color="royalblue", weight=10];
      rankdir="LR";
      size = "16,10";

      physics_solver [label="PhysicsSolver", href="https://lethe-cfd.github.io/lethe/html_doxygen/classPhysicsSolver.html"];

      navier_stokes_base [label="NavierStokesBase",href="https://lethe-cfd.github.io/lethe/html_doxygen/classNavierStokesBase.html"];

      auxiliary_physics [label="AuxiliaryPhysics",href="https://lethe-cfd.github.io/lethe/html_doxygen/classAuxiliaryPhysics.html"];

      physics_solver:e -> navier_stokes_base:w [dir=back];
      physics_solver:e -> auxiliary_physics:w [dir=back];

      navier_stokes_base_1 [label=<<B>GLSNavierStokesSolver</B> <br/>(lethe-fluid)>,href="https://lethe-cfd.github.io/lethe/html_doxygen/classGLSNavierStokesSolver.html", tooltip="GLSNavierStokesSolver"];
      navier_stokes_base_2 [label=<<B>GDNavierStokesSolver</B> <br/> (lethe-fluid-block)>,href="https://lethe-cfd.github.io/lethe/html_doxygen/classGDNavierStokesSolver.html", tooltip="GDNavierStokesSolver"];
      navier_stokes_base_3 [label=<<B>MFNavierStokesSolver</B> <br/> (lethe-fluid-matrix-free)>,href="https://lethe-cfd.github.io/lethe/html_doxygen/classMFNavierStokesSolver.html", tooltip="MFNavierStokesSolver"];

      navier_stokes_base:e -> navier_stokes_base_1:w [dir=back];
      navier_stokes_base:e -> navier_stokes_base_2:w [dir=back];
      navier_stokes_base:e -> navier_stokes_base_3:w [dir=back];

      auxiliary_physics_1 [label="VolumeOfFluid",href="https://lethe-cfd.github.io/lethe/html_doxygen/classVolumeOfFluid.html"];
      auxiliary_physics_2 [label="CahnHilliard",href="https://lethe-cfd.github.io/lethe/html_doxygen/classCahnHilliard.html"];
      auxiliary_physics_3 [label="HeatTransfer",href="https://lethe-cfd.github.io/lethe/html_doxygen/classHeatTransfer.html"];
      auxiliary_physics_4 [label="Tracer",href="https://lethe-cfd.github.io/lethe/html_doxygen/classTracer.html"];

      auxiliary_physics:e -> auxiliary_physics_1:w [dir=back];
      auxiliary_physics:e -> auxiliary_physics_2:w [dir=back];
      auxiliary_physics:e -> auxiliary_physics_3:w [dir=back];
      auxiliary_physics:e -> auxiliary_physics_4:w [dir=back];

      navier_stokes_base_1_1 [label=<<B>GLSVANSSolver</B> <br/>(lethe-fluid-vans)>,href="https://lethe-cfd.github.io/lethe/html_doxygen/classGLSVANSSolver.html", tooltip="GLSVANSSolver"];
      navier_stokes_base_1_2 [label=<<B>GLSSharpNavierStokesSolver</B> <br/>(lethe-fluid-sharp)>,href="https://lethe-cfd.github.io/lethe/html_doxygen/classGLSSharpNavierStokesSolver.html", tooltip="GLSSharpNavierStokesSolver"];
      navier_stokes_base_1_3 [label=<<B>GLSNitscheNavierStokesSolver</B> <br/>(lethe-fluid-nitsche)>,href="https://lethe-cfd.github.io/lethe/html_doxygen/classGLSNitscheNavierStokesSolver.html", tooltip="GLSNitscheNavierStokesSolver"];

      navier_stokes_base_1:e -> navier_stokes_base_1_1:w [dir=back];
      navier_stokes_base_1:e -> navier_stokes_base_1_2:w [dir=back];
      navier_stokes_base_1:e -> navier_stokes_base_1_3:w [dir=back];

      navier_stokes_base_1_1_1 [label=<<B>CFDDEMSolver</B> <br/>(lethe-fluid-particles)>,href="https://lethe-cfd.github.io/lethe/html_doxygen/classNavierStokesBase.html", tooltip="CFDDEMSolver"];

      navier_stokes_base_1_1:e -> navier_stokes_base_1_1_1:w [dir=back];

      dem_solver [label=<<B>DEM</B> <br/>(lethe-particles)>, href="https://lethe-cfd.github.io/lethe/html_doxygen/classDEMSolver.html", tooltip="DEM"];

      rpt_solver_1 [label=<<B>RPT</B> <br/>(lethe-rpt-3d)>, href="https://lethe-cfd.github.io/lethe/html_doxygen/classRPT.html", tooltip="RPT"];
      rpt_solver_2 [label=<<B>RPTCellReconstruction</B> <br/>(lethe-rpt-cell-reconstruction-3d)>, href="https://lethe-cfd.github.io/lethe/html_doxygen/classRPTCellReconstruction.html", tooltip="RPTCellReconstruction"];
      rpt_solver_3 [label=<<B>RPTFEMReconstruction</B> <br/>(lethe-rpt-fem-reconstruction-3d)>, href="https://lethe-cfd.github.io/lethe/html_doxygen/classRPTFEMReconstruction.html", tooltip="RPTFEMReconstruction"];
      rpt_solver_4 [label=<<B>RPTL2Projection</B> <br/>(lethe-rpt-l2-projection-3d)>, href="https://lethe-cfd.github.io/lethe/html_doxygen/classRPTL2Projection.html", tooltip="RPTL2Projection"];
    }
 * @enddot
 */
