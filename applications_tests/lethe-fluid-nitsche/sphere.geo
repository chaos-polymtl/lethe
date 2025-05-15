SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
Physical Surface(0) = {1};
Physical Volume(0) = {1};
Mesh.MeshSizeFactor = 2;
