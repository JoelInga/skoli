Point(1) = {0.0, 0.0, 0.0, 0.0005};
Point(2) = {0.0, 0.005021996212849453, 0.0, 0.0005};
Point(3) = {0.003, 0.005021996212849453, 0.0, 0.0005};
Point(4) = {0.003, 0.0, 0.0, 0.0005};
Curve(1) = {1, 2};
Curve(2) = {2, 3};
Curve(3) = {3, 4};
Curve(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Physical Groups
Physical Curve("innri", 1) = {1};
Physical Curve("ytri", 2) = {3, 4};
Physical Surface("ribba", 3) = {1};
Physical Curve("einangrun", 4) = {2};
Coherence;
