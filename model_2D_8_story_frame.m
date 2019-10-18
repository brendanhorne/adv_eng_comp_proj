total_nodes = 45;
total_members = 72;

%Nodes
%nodes = [node no. , x-coord, y-coord]

nodes = zeros(total_nodes,3);
    for i = 1:45
        %Node Number
        nodes(i,1) = i;
        
        %x coordinate
        nodes(i,2) = (rem((i-1),5)) * 5;
        
        %y coordinates
        nodes(i,3) = floor((i-1)/5)*3.5;
        
        %z rotation
        nodes(i,4) = 0;
        
        % forces
        nodes(i,5:7) = [0 0 0];
    end
    
%Elements
%elements = [ elements no. , starting node, ending node, element type]

epf = 9; %elements per floor
bpf = 4; %Beams per floor
cpf = 5; %Columns per floor
bsi = cpf + 1; % beam starting index
n = 0; %floor number

E = 35e9;
ro = 3000;
k = 4e9;
Ic = (0.4^4)/12;
yb = (0.25*0.45*0.225 + 0.6*0.15*0.525) /  (0.25*0.45 + 0.6*0.15);
Ib =  0.25*0.45*(yb - 0.45/2)^2 + (0.25*0.45^3)/12 + 0.6*0.15*(0.45 + 0.15/2 - yb)^2 + (0.6*0.15^3)/12;
Is = (2.5^3*0.4)/12;
Ac = 0.16;
Ab = 0.6*0.15+0.45*0.25;
As = 1;
    
elements = zeros(total_members,13);
    for i = 1:total_members
        elements(i,1) = i;
        n = i * 1/9;
        elements(i,2) = i - floor((i - 1)/9)*4;
        elements(i,3) = i + 5 - floor((i + 3)/9)*4;
        if  elements (i, 3) - elements(i,2) < 5 % beam
            elements(i,10) = E;
            elements(i,11) = Ib;
            elements(i,12) = Ab;
            elements(i,13) = ro;
        elseif rem((i-3),9) == 0 % shear wall
            elements(i,10) = E;
            elements(i,11) = Is;
            elements(i,12) = As;
            elements(i,13) = ro;
        else
            elements(i,10) = E;
            elements(i,11) = Ic;
            elements(i,12) = Ac;
            elements(i,13) = ro;
        end
    end

node46 = [46 10 0 0 0 0 0]; %add in node for spring
nodes = [nodes; node46];   

springs = [% spring ID, start node, end node, x-k, y-k, z-k (GLOBAL)
            1 3 46 0 4e9 0];
        
boundary_conditions = [% node ID, x-position-fixity, y-position-fixity, z-rotation-fixity
            1 1 1 1;
            2 1 1 1;
            3 1 0 1;
            4 1 1 1;
            5 1 1 1;
            46 1 1 1];