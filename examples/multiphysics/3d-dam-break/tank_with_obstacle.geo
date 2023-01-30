// Gmsh project created on Fri Nov 11 15:28:42 2022
SetFactory("Built-in");

// ==============================================
// Define variables
// ==============================================

lc = 1.0;
f = 1;

// Geometry variables
left = 0;
right = 3.22;
left_box = 0.6635;
right_box = 0.8245;
front = -0.5;
back = 0.5;
box_front = -0.2015;
box_back = 0.2015;
bottom = 0;
top = 1;
box_top = 0.161;
height_bottom = 0.161;
height_top = 0.839;

// Division variables
ny_back   = 1*f;
nx_left   = 4*f;
nx_box    = 1*f;
nx_right  = 14*f;
ny_front  = 1*f;
ny_center = 2*f;
nz_bottom = 1*f;
nz_top    = 5*f;


// ==============================================
// z planes
// ==============================================

// Bottom level
Point(1)  = {left, back, bottom,lc};
Point(2)  = {left_box, back, bottom,lc};
Point(3)  = {left_box, box_back, bottom,lc};
Point(4)  = {left, box_back, bottom,lc};
Point(5)  = {left, box_front, bottom, lc};
Point(6)  = {left_box, box_front, bottom, lc};
Point(7)  = {left_box, front, bottom, lc};
Point(8)  = {left, front, bottom,lc};
Point(9)  = {right_box, front, bottom, lc};
Point(10) = {right_box, box_front, bottom, lc};
Point(11) = {right, box_front, bottom, lc};
Point(12) = {right, front, bottom, lc};
Point(13) = {right, box_back, bottom, lc};
Point(14) = {right_box, box_back, bottom, lc};
Point(15) = {right_box, back, bottom, lc};
Point(16) = {right, back, bottom, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {4,5};
Line(6) = {5,6};
Line(7) = {6,3};
Line(8)  = {6,7};
Line(9)  = {7,8};
Line(10) = {8,5};
Line(11) = {7,9};
Line(12) = {9,10};
Line(13) = {10,6};
Line(14) = {10,11};
Line(15) = {11,12};
Line(16) = {12,9};
Line(17) = {11,13};
Line(18) = {13,14};
Line(19) = {14,10};
Line(20) = {14,15};
Line(21) = {15,16};
Line(22) = {16,13};
Line(23) = {15,2};
Line(24) = {3,14};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {-3,-7,-6,-5};
Line Loop(3) = {6,8,9,10};
Line Loop(4) = {-13,-12,-11,-8};
Line Loop(5) = {14,15,16,12};
Line Loop(6) = {-18,-17,-14,-19};
Line Loop(7) = {21,22,18,20};
Line Loop(8) = {-23,-20,-24,-2};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};

// Top
Point(17) = {left, back, top,lc};
Point(18) = {left_box, back, top,lc};
Point(19) = {left_box, box_back, top,lc};
Point(20) = {left, box_back, top,lc};
Point(21) = {left, box_front, top, lc};
Point(22) = {left_box, box_front, top, lc};
Point(23) = {left_box, front, top, lc};
Point(24) = {left, front, top,lc};
Point(25) = {right_box, front, top, lc};
Point(26) = {right_box, box_front, top, lc};
Point(27) = {right, box_front, top, lc};
Point(28) = {right, front, top, lc};
Point(29) = {right, box_back, top, lc};
Point(30) = {right_box, box_back, top, lc};
Point(31) = {right_box, back, top, lc};
Point(32) = {right, back, top, lc};

Line(25) = {17,18};
Line(26) = {18,19};
Line(27) = {19,20};
Line(28) = {20,17};
Line(29) = {20,21};
Line(30) = {21,22};
Line(31) = {22,19};
Line(32) = {22,23};
Line(33) = {23,24};
Line(34) = {24,21};
Line(35) = {23,25};
Line(36) = {25,26};
Line(37) = {26,22};
Line(38) = {26,27};
Line(39) = {27,28};
Line(40) = {28,25};
Line(41) = {27,29};
Line(42) = {29,30};
Line(43) = {30,26};
Line(44) = {30,31};
Line(45) = {31,32};
Line(46) = {32,29};
Line(47) = {31,18};
Line(48) = {19,30};

Line Loop(9)  = {25,26,27,28};
Line Loop(10) = {27,31,30,29};
Line Loop(11) = {30,32,33,34};
Line Loop(12) = {37,36,35,32};
Line Loop(13) = {38,39,40,36};
Line Loop(14) = {42,41,38,43};
Line Loop(15) = {45,46,42,44};
Line Loop(16) = {47,44,48,26};
Line Loop(17) = {48,43,37,31};

Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};
Plane Surface(14) = {14};
Plane Surface(15) = {15};
Plane Surface(16) = {16};
Plane Surface(17) = {17};

// Box_top
Point(33) = {left, back, box_top, lc};
Point(34) = {left_box, back, box_top,lc};
Point(35) = {left_box, box_back, box_top,lc};
Point(36) = {left, box_back, box_top,lc};
Point(37) = {left, box_front, box_top, lc};
Point(38) = {left_box, box_front, box_top, lc};
Point(39) = {left_box, front, box_top, lc};
Point(40) = {left, front, box_top, lc};
Point(41) = {right_box, front, box_top, lc};
Point(42) = {right_box, box_front, box_top, lc};
Point(43) = {right, box_front, box_top, lc};
Point(44) = {right, front, box_top, lc};
Point(45) = {right, box_back, box_top, lc};
Point(46) = {right_box, box_back, box_top, lc};
Point(47) = {right_box, back, box_top, lc};
Point(48) = {right, back, box_top, lc};

Line(49) = {33,34};
Line(50) = {47,34};
Line(51) = {47,48};
Line(52) = {48,45};
Line(53) = {43,45};
Line(54) = {43,44};
Line(55) = {44,41};
Line(56) = {39,41};
Line(57) = {39,40};
Line(58) = {40,37};
Line(59) = {36,37};
Line(60) = {36,33};
Line(61) = {35,46};
Line(62) = {46,42};
Line(63) = {42,38};
Line(64) = {38,35};
Line(73) = {35,36};
Line(74) = {37,38};
Line(75) = {45,46};
Line(76) = {42,43};
Line(77) = {34,35};
Line(78) = {46,47};
Line(79) = {38,39};
Line(80) = {41,42};

Line Loop(18) = {61,62,63,64}; // Box top surface
Line Loop(25) = {49,77,73,60};
Line Loop(26) = {50,78,61,77};
Line Loop(27) = {51,52,75,78};
Line Loop(28) = {75,53,76,62};
Line Loop(29) = {76,54,55,80};
Line Loop(30) = {63,80,56,79};
Line Loop(31) = {79,57,58,74};
Line Loop(32) = {73,64,74,59};

Plane Surface(18) = {18};
Plane Surface(25) = {25};
Plane Surface(26) = {26};
Plane Surface(27) = {27};
Plane Surface(28) = {28};
Plane Surface(29) = {29};
Plane Surface(30) = {30};
Plane Surface(31) = {31};
Plane Surface(32) = {32};

// ==============================================
// y planes
// ==============================================

// Front
Line(65) = {40,24};
Line(66) = {23,39};
Line(67) = {41,25};
Line(68) = {28,44};
Line(69) = {8,40};
Line(70) = {39,7};
Line(71) = {9,41};
Line(72) = {44,12};

Line Loop(19) = {-33,66,57,65};
Line Loop(20) = {35,-67,-56,-66};
Line Loop(21) = {-40,68,55,67};
Line Loop(22) = {-55,72,16,71};
Line Loop(23) = {56,-71,-11,-70};
Line Loop(24) = {-57,70,9,69};

Plane Surface(19) = {19};
Plane Surface(20) = {20};
Plane Surface(21) = {21};
Plane Surface(22) = {22};
Plane Surface(23) = {23};
Plane Surface(24) = {24};

// Box_front
Line(81) = {37,21};
Line(82) = {22,38};
Line(83) = {42,26};
Line(84) = {27,43};
Line(85) = {43,11};
Line(86) = {10,42};
Line(87) = {38,6};
Line(88) = {5,37};

Line Loop(33) = {30,82,-74,81};
Line Loop(34) = {37,82,-63,83};
Line Loop(35) = {38,84,-76,83};
Line Loop(36) = {76,85,-14,86};
Line Loop(37) = {63,87,-13,86};
Line Loop(38) = {74,87,-6,88};

Plane Surface(33) = {33};
Plane Surface(34) = {34};
Plane Surface(35) = {35};
Plane Surface(36) = {36};
Plane Surface(37) = {37};
Plane Surface(38) = {38};

// Box_back
Line(89) = {36,20};
Line(90) = {19,35};
Line(91) = {46,30};
Line(92) = {29,45};
Line(93) = {45,13};
Line(94) = {14,46};
Line(95) = {35,3};
Line(96) = {4,36};

Line Loop(39) = {-27,90,73,89};
Line Loop(40) = {-48,90,61,91};
Line Loop(41) = {-42,92,75,91};
Line Loop(42) = {-75,93,18,94};
Line Loop(43) = {-61,95,24,94};
Line Loop(44) = {-73,95,3,96};

Plane Surface(39) = {39};
Plane Surface(40) = {40};
Plane Surface(41) = {41};
Plane Surface(42) = {42};
Plane Surface(43) = {43};
Plane Surface(44) = {44};

// Back
Line(97)  = {33,17};
Line(98)  = {18,34};
Line(99)  = {47,31};
Line(100) = {32,48};
Line(101) = {48,16};
Line(102) = {15,47};
Line(103) = {34,2};
Line(104) = {1,33};

Line Loop(45) = {25,98,-49,97};
Line Loop(46) = {47,98,-50,99};
Line Loop(47) = {45,100,-51,99};
Line Loop(48) = {51,101,-21,102};
Line Loop(49) = {50,103,-23,102};
Line Loop(50) = {49,103,-1,104};

Plane Surface(45) = {45};
Plane Surface(46) = {46};
Plane Surface(47) = {47};
Plane Surface(48) = {48};
Plane Surface(49) = {49};
Plane Surface(50) = {50};

// ==============================================
// x planes
// ==============================================

// Left
Line Loop(51) = {-28,-89,60,97};
Line Loop(52) = {-29,-89,59,81};
Line Loop(53) = {-34,-65,58,81};
Line Loop(54) = {-58,-69,10,88};
Line Loop(55) = {-59,-96,5,88};
Line Loop(56) = {-60,-96,4,104};

Plane Surface(51) = {51};
Plane Surface(52) = {52};
Plane Surface(53) = {53};
Plane Surface(54) = {54};
Plane Surface(55) = {55};
Plane Surface(56) = {56};

// Left_box
Line Loop(57) = {26,90,-77,-98};
Line Loop(58) = {31,90,-64,-82};
Line Loop(59) = {32,66,-79,-82};
Line Loop(60) = {79,70,-8,-87};
Line Loop(61) = {64,95,-7,-87};
Line Loop(62) = {77,95,-2,-103};

Plane Surface(57) = {57};
Plane Surface(58) = {58};
Plane Surface(59) = {59};
Plane Surface(60) = {60};
Plane Surface(61) = {61};
Plane Surface(62) = {62};

// Right_box
Line Loop(63) = {-44,-91,78,99};
Line Loop(64) = {-43,-91,62,83};
Line Loop(65) = {-36,-67,80,83};
Line Loop(66) = {-80,-71,12,86};
Line Loop(67) = {-62,-94,19,86};
Line Loop(68) = {-78,-94,20,102};

Plane Surface(63) = {63};
Plane Surface(64) = {64};
Plane Surface(65) = {65};
Plane Surface(66) = {66};
Plane Surface(67) = {67};
Plane Surface(68) = {68};

// Right
Line Loop(69) = {46,92,-52,-100};
Line Loop(70) = {41,92,-53,-84};
Line Loop(71) = {39,68,-54,-84};
Line Loop(72) = {54,72,-15,-85};
Line Loop(73) = {53,93,-17,-85};
Line Loop(74) = {52,93,-22,-101};

Plane Surface(69) = {69};
Plane Surface(70) = {70};
Plane Surface(71) = {71};
Plane Surface(72) = {72};
Plane Surface(73) = {73};
Plane Surface(74) = {74};

// ==============================================
// Volumes
// ==============================================

// Bottom
Surface Loop(1) = {1,50,56,44,62,25};
Surface Loop(2) = {2,44,55,38,61,32};
Surface Loop(3) = {3,38,54,24,60,31};
Surface Loop(4) = {4,37,60,23,66,30};
Surface Loop(5) = {5,36,66,22,72,29};
Surface Loop(6) = {6,42,67,36,73,28};
Surface Loop(7) = {7,48,68,42,74,27};
Surface Loop(8) = {8,49,62,43,68,26};

Volume(1) = {1};
Volume(2) = {2};
Volume(3) = {3};
Volume(4) = {4};
Volume(5) = {5};
Volume(6) = {6};
Volume(7) = {7};
Volume(8) = {8};

// Top
Surface Loop(9)  = {25,45,51,39,57,9};
Surface Loop(10) = {32,39,52,33,58,10};
Surface Loop(11) = {31,33,53,19,59,11};
Surface Loop(12) = {30,34,59,20,65,12};
Surface Loop(13) = {29,35,65,21,71,13};
Surface Loop(14) = {28,41,64,35,70,14};
Surface Loop(15) = {27,47,63,41,69,15};
Surface Loop(16) = {26,46,57,40,63,16};
Surface Loop(17) = {18,40,58,34,64,17};

Volume(9)  = {9};
Volume(10) = {10};
Volume(11) = {11};
Volume(12) = {12};
Volume(13) = {13};
Volume(14) = {14};
Volume(15) = {15};
Volume(16) = {16};
Volume(17) = {17};

// ==============================================
// Hexahedral Elements and Refininement
// ==============================================

Transfinite Line {1,3,6,9,25,27,30,33,49,57,73,74} = Ceil(nx_left+1) Using Progression 1;
Transfinite Line {13,11,23,24,37,35,47,48,50,56,61,63} = Ceil(nx_box+1) Using Progression 1;
Transfinite Line {14,16,18,21,38,40,42,45,51,55,75,76} = Ceil(nx_right+1) Using Progression 1;

Transfinite Line {2,4,20,22,28,26,44,46,52,60,77,78} = Ceil(ny_back+1) Using Progression 1;
Transfinite Line {5,7,17,19,29,31,41,43,53,59,62,64} = Ceil(ny_center+1) Using Progression 1;
Transfinite Line {10,8,12,15,34,32,36,39,54,58,79,80} = Ceil(ny_front+1) Using Progression 1;

Transfinite Line {65:68,81:84,89:92,97:100,105:108,113:116,121:124,129:132} = Ceil(nz_top+1) Using Progression 1;
Transfinite Line {69:72,85:88,93:96,101:104,109:112,117:120,125:128,133:136} = Ceil(nz_bottom+1) Using Progression 1;

Transfinite Surface{1:74};
Recombine Surface{1:74};

Transfinite Volume{1:17};
Recombine Volume{1:17};

// ==============================================
// Domain Boundaries and Final Volume
// ==============================================

Physical Surface(0) = {1:8};  	        // Bottom
Physical Surface(1) = {9:17}; 	        // Top
Physical Surface(2) = {51:56}; 	        // Left
Physical Surface(3) = {69:74}; 	        // Right
Physical Surface(4) = {19:24};          // Front
Physical Surface(5) = {45:50};	        // Back
Physical Surface(6) = {18,37,43,61,67}; // Box (obstacle)

Physical Volume(0) = {1:17};
