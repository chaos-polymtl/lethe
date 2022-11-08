// Define a variable

H=0.45; // height of bed
W=0.07; // half width of bed from end of channel to end of bed
Z = 0.015; // depth of the bed in z
iW=0.005; // half width of channel
B=0.04; // point at channel inlet
nw=8;  // number of points  in the half width of bed after channel
nwi=2; // number of points in width of channel
nb=5; // number of points in the height of the channel
nh=46; // number of points in the height of the bed without channel
zl = 1; // number of cells in the z direction (depth)

Point(2) = {-iW, -B, 0};
Point(3) = {+iW, -B, 0};
Point(11) = {-W-iW, 0, 0};
Point(12) = {-iW, 0, 0};
Point(13) = {iW,0, 0};
Point(14) = {iW+W, 0, 0};
Point(21) = {-W-iW, H, 0};
Point(22) = {-iW, H, 0};
Point(23) = {iW,H, 0};
Point(24) = {iW+W, H, 0};

Line(1)={11,12};
Line(2)={12,2};
Line(3)={2,3};
Line(4)={3,13};
Line(5)={13,14};
Line(6)={14,24};
Line(7)={24,23};
Line(8)={23,22};
Line(9)={22,21};
Line(10)={21,11};

Line(101)={12,13};
Line(102)={12,22};
Line(103)={13,23};

Line Loop(1) = {1,102,9,10};
Line Loop(2) = {101,103,8,-102};
Line Loop(3) = {5,6,7,-103};
Line Loop(4) = {2,3,4,-101};
Plane Surface(1) = {1} ;
Plane Surface(2) = {2} ;
Plane Surface(3) = {3} ;
Plane Surface(4) = {4} ;
Transfinite Surface {1,2,3,4};
Recombine Surface{1,2,3,4};

Extrude {0, 0, Z} {
  Surface{1}; Surface{3}; Surface{2}; Surface{4}; Layers{zl}; Recombine;
}

Transfinite Line {10,102,103,6} = Ceil(nh) Using Progression 1;
Transfinite Line {1,5,9,7} = Ceil(nw) Using Progression 1;
Transfinite Line {3,101,8} = Ceil(nwi) Using Progression 1;
Transfinite Line {2,4} = Ceil(nb) Using Progression 1;

Physical Volume(0) = {1:4};

Physical Surface(0) = {1,2,3,4,124,138,125,191,169,147,178,186};
Physical Surface(1) = {182};
Physical Surface(2) = {120,142,164};
Physical Surface(3) = {112,134};
