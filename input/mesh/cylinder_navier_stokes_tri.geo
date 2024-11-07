SetFactory("OpenCASCADE");

cl1 = 1.0;

//+
Point(1) = {0.0, -0.0, 0, cl1};
//+
Point(2) = {2.2, -0.0, 0, cl1};
//+
Point(3) = {2.2, 0.41, 0, cl1};
//+
Point(4) = {0.0, 0.41, 0, cl1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {0.2, 0.2, 0, 0.05, 0, 2*Pi};

//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};


//+
Physical Curve("bottom", 6) = {1};
//+
Physical Curve("right", 7) = {2};
//+
Physical Curve("top", 8) = {3};
//+
Physical Curve("left", 9) = {4};
//+
Physical Curve("cylinder", 10) = {5};
//+
Physical Surface("fluid", 11) = {1};

//+
Transfinite Curve {1, -3} = 50 Using Progression 1;
Transfinite Curve {2, -4} = 10 Using Progression 1;
Transfinite Curve {5} = 20 Using Progression 1;