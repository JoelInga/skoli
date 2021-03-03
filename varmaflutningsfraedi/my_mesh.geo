Point(1) = {0, -0.005, 0.0, 0.0002};
Point(2) = {0.0, 0.005, 0.0, 0.0002};
Point(3) = {0.003, 0.005, 0.0, 0.0002};
Point(4) = {0.003, 0.0015, 0.0, 0.0002};
Point(5) = {0.013000000000000001, 0.0015, 0.0, 0.0002};
Point(6) = {0.013000000000000001, -0.0015, 0.0, 0.0002};
Point(7) = {0.003, -0.0015, 0.0, 0.0002};
Point(8) = {0.003, -0.005, 0.0, 0.0002};
Curve(1) = {1, 2};
Curve(2) = {2, 3};
Curve(3) = {3, 4};
Curve(4) = {4, 5};
Curve(5) = {5, 6};
Curve(6) = {6, 7};
Curve(7) = {7, 8};
Curve(8) = {8, 1};
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// Physical Groups
Physical Curve("innri", 1) = {1};
Physical Curve("ytri", 2) = {3, 4, 5, 6, 7};
Physical Surface("ribba", 3) = {1};
Physical Curve("einangrun", 4) = {2, 8};
Coherence;
