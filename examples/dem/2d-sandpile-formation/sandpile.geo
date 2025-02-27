// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

SetFactory("Built-in");


// Geometry variables
upper_right    = 0.505;
upper_left     = - upper_right;
upper_height   = 1;
angle_height   = 0.57;
channel_right  = 0.035;
channel_left   = - channel_right;
upper_channel_height = 0;
bottom_channel_height = -0.22;
bottom_right   = 1;
bottom_left    = - bottom_right;
bottom_height  = -0.57;


//+
Point(1) = {upper_left, upper_height, 0, 1.0};
//+
Point(2) = {channel_left, upper_height, 0, 1.0};
//+
Point(3) = {channel_right, upper_height, 0, 1.0};
//+
Point(4) = {upper_right, upper_height, 0, 1.0};
//+
Point(5) = {upper_right, angle_height, 0, 1.0};
//+
Point(6) = {channel_right, upper_channel_height, 0, 1.0};
//+
Point(7) = {channel_right, bottom_channel_height, 0, 1.0};
//+
Point(8) = {bottom_right, bottom_channel_height, 0, 1.0};
//+
Point(9) = {bottom_right, bottom_height, 0, 1.0};
//+
Point(10) = {channel_right, bottom_height, 0, 1.0};
//+
Point(11) = {channel_left, bottom_height, 0, 1.0};
//+
Point(12) = {bottom_left, bottom_height, 0, 1.0};
//+
Point(13) = {bottom_left, bottom_channel_height, 0, 1.0};
//+
Point(14) = {channel_left, bottom_channel_height, 0,  1.0};
//+
Point(15) = {channel_left, upper_channel_height, 0, 1.0};
//+
Point(16) = {upper_left, angle_height, 0, 1.0};
//+

//0.57 au lieu de 0.6 sur Points 5 et 16 ie 50.49 degres au lieu de 51.92 (52 dans le papier) pour eviter les pertes de particules

// Upper left
Line(1) = {1, 2};
Line(2) = {2, 15};
Line(3) = {15, 16};
Line(4) = {16, 1};

Transfinite Line {1,3} = 5 Using Progression 1;
Transfinite Line {2,4} = 6 Using Progression 1;
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Transfinite Surface {1};
Physical Surface(1) = {1};

// Upper Right
Line(5) = {3, 4};
Line(6) = {4, 5};
Line(7) = {5, 6};
Line(8) = {6, 3};

Transfinite Line {5,7} = 5 Using Progression 1;
Transfinite Line {6,8} = 6 Using Progression 1;
Line Loop(2) = {5,6,7,8};
Plane Surface(2) = {2};

Transfinite Surface {2};
Physical Surface(2) = {2};

//Upper Middle
Line(9) = {2, 3};
Line(10) = {6, 15};

Transfinite Line {9,10} = 2 Using Progression 1;
Transfinite Line {2, 8} = 6 Using Progression 1;
Line Loop(3) = {9,-8,10,-2};
Plane Surface(3) = {3};

Transfinite Surface {3};
Physical Surface(3) = {3};

// Bottom Left
Line(11) = {13, 14};
Line(12) = {14, 11};
Line(13) = {11, 12};
Line(14) = {12, 13};

Transfinite Line {11,13} = 6 Using Progression 1;
Transfinite Line {12,14} = 5 Using Progression 1;
Line Loop(4) = {11,12,13,14};
Plane Surface(4) = {4};

Transfinite Surface {4};
Physical Surface(4) = {4};


// Bottom Right
Line(15) = {7, 8};
Line(16) = {8, 9};
Line(17) = {9, 10};
Line(18) = {10, 7};

Transfinite Line {15,17} = 6 Using Progression 1;
Transfinite Line {16,18} = 5 Using Progression 1;
Line Loop(5) = {15,16,17,18};
Plane Surface(5) = {5};

Transfinite Surface {5};
Physical Surface(5) = {5};

// Bottom Middle
Line(19) = {14, 7};
Line(20) = {10, 11};

Transfinite Line {19,20} = 2 Using Progression 1;
Transfinite Line {12,18} = 5 Using Progression 1;
Line Loop(6) = {19,-18,20,-12};
Plane Surface(6) = {6};

Transfinite Surface {6};
Physical Surface(6) = {6};

// Middle Middle
Line(21) = {6, 7};
Line(22) = {14, 15};

Transfinite Line {10,19} = 2 Using Progression 1;
Transfinite Line {21,22} = 4 Using Progression 1;
Line Loop(7) = {-10, 21, -19, 22};
Plane Surface(7) = {7};

Transfinite Surface {7};
Physical Surface(7) = {7};


// Reconstruction

Physical Line(1) = {1,9,5};
Physical Line(2) = {6,7};
Physical Line(3) = {21};
Physical Line(4) = {15};
Physical Line(5) = {16};
Physical Line(6) = {17,20,13};
Physical Line(7) = {14};
Physical Line(8) = {11};
Physical Line(9) = {22};
Physical Line(10) = {3,4};

Recombine Surface{1,2,3,4,5,6,7};