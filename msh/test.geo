
n = 1;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {1.0, 0.0, 0.0};
Point(3) = {1.0, 1.0, 0.0};
Point(4) = {0.0, 1.0, 0.0};
Point(5) = {0.25, 0.0, 0.0};
Point(6) = {0.0, 0.25, 0.0};
Point(7) = {0.25, 0.25, 0.0};
// Point(8) = {0.75, 0.75, 0.0};

Line(1) = {1,5};
Line(2) = {5,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,6};
Line(6) = {6,1};
Line(7) = {5,7};
Line(8) = {7,6};
Line(9) = {7,3};

Curve Loop(1) = {1,7,8,6};
Curve Loop(2) = {-7,2,3,-9};
Curve Loop(3) = {-8,9,4,5};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
// Curve{9} In Surface{1};
// Curve{5,6,7,8} In Surface{1};
Transfinite Curve{2,3,4,5,9} = n+1;
Transfinite Curve{1,6,7,8} = 3*n+1;
Transfinite Surface{1};

Physical Curve("Γᵍ") = {1,2,3,4,5,6};
Physical Surface("Ω") = {1,2,3};

Mesh.Algorithm = 6;
Mesh.MshFileVersion = 2;
Mesh 2;