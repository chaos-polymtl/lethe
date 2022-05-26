Mesh.Algorithm = 6;
Mesh.Hexahedra = 1;
Mesh.HighOrderNumLayers = 6;
Mesh.MeshSizeFactor = 1.5;
Mesh.MeshSizeMin = 0.3;
Mesh.MeshSizeMax = 0.4;
Mesh.MetisAlgorithm = 1;
Mesh.MetisEdgeMatching = 2;
Mesh.MetisRefinementAlgorithm = 2;
Mesh.Optimize = 1;
Mesh.OptimizeThreshold = 0.3;
Mesh.RecombinationAlgorithm = 2;
Mesh.RecombineAll = 1;
Mesh.RefineSteps = 10;
Mesh.Smoothing = 10;
Mesh.SmoothRatio = 1.8;

//+
Point(1) = {-0.0625, -0.0056, 0.2, 1.0};
//+
Point(2) = {0.0625, -0.0056, 0.2, 1.0};
//+
Point(3) = {0.0625, -0.0056, 0.05, 1.0};
//+
Point(4) = {0.0125, -0.0056, 0, 1.0};
//+
Point(5) = {0.0625, -0.0056, -0.02, 1.0};
//+
Point(6) = {0.0625, -0.0056, -0.2, 1.0};
//+
Point(7) = {-0.0625, -0.0056, -0.2, 1.0};
//+
Point(8) = {-0.0625, -0.0056, -0.02, 1.0};
//+
Point(9) = {-0.0125, -0.0056, 0, 1.0};
//+
Point(10) = {-0.0625, -0.0056, 0.05, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0.0056, 0} {
  Surface{1}; Layers{1}; Recombine;
}
