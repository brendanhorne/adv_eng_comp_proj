nodes = [% node ID, x-coordinate, y-coordinate, z-rotation, x-force, y-force, z-moment
            1 0  0 0 0 0 0;
            2 0  1 0 0 0 0];


elements = [% element ID, start node, end node, x1-force, y1-force, z1-moment, x2-force, y2-force, z2-moment, E, I, A, rho
            1 1 2 0 0 0 0 0 0 200e9 8.33333333333333e-6 0.01, 3e3];

boundary_conditions = [% node ID, x-position-fixity, y-position-fixity, z-rotation-fixity
            1 1 1 1];
        
% springs = [1 1 2 0 4e9 0];        
springs = [];   