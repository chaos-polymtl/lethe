

SetFactory("OpenCASCADE");

H =  17.4625;   // Tank height
T = 28.57;  // Tank Diameter
z0= 26-1;

x0=2;//-2.835941-1.10379368685;
y0=-0.5;//2.554151-1.662944;

C = T/4;   // Clearance

Mesh.CharacteristicLengthMin = 1e-3;
Mesh.CharacteristicLengthMax = 2.5e-3;
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;
Mesh.HighOrderOptimize = 1;

Cylinder(1) = {0, 0, z0, 0, 0, H, T/2, 	2*Pi}; 

Merge "db_helical_Tiff.step";
Translate {x0, y0, 0} { Volume{2:6}; }

BooleanDifference(7)= { Volume{1}; Delete; }{ Volume{2:6}; Delete; };
Dilate { {0 ,0 ,0} , 1e-2 }{  Volume{7}; }

Physical Surface(1) = {1,3};
Physical Surface(2) = {2};
Physical Surface(3) = {4:28};
Physical Volume(0) = {7};


// 		xcenter, ycenter, ycenter, 	xaxis, yaxis, zaxis,	radius, 	angular opening
