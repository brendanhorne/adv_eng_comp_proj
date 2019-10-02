clc, clear, close
% All units in SI kg, m, s, N, Pa etc.

% INPUT SPACE
nodes_out = fopen('data/single_column/node_displacements.csv','w');
fprintf(nodes_out,'delta_t,node_id,x,y,z\r\n');
nodes = [% node ID, x-coordinate, y-coordinate, z-rotation, x-force, y-force, z-moment
            1 0 0 0 0 0 0;
            2 0 1 0 0 0 0];

elements_out = fopen('data/single_column/element_forces.csv','w');
fprintf(elements_out,'delta_t,element_id,N1,V1,M1,N2,V2,M2\r\n');
elements = [% element ID, start node, end node, x1-force, y1-force, z1-moment, x2-force, y2-force, z2-moment, E, I, A, rho
            1 1 2 0 0 0 0 0 0 200e6 8.33333333333333e-6 0.01, 3e3];

boundary_conditions = [% node ID, x-position-fixity, y-position-fixity, z-rotation-fixity
            1 1 1 1];

acceleration = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
                1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 ...
                0 -0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 ...
                -1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 ...
                ];
        
% PROGRAM SPACE

% form column vectors of element properties
E = elements(1:end,10); 
I = elements(1:end,11);
A = elements(1:end,12);
rho = elements(1:end,13);
        
% determine the number of boundary conditions
boundary_conditions_only = boundary_conditions(1:end,2:end);
number_possible_fixed_dof = numel(boundary_conditions_only);
boundary_conditions_list = reshape(boundary_conditions_only,number_possible_fixed_dof,1);
number_of_fixed_dof = 0;
for i = 1:size(boundary_conditions_list)
    if boundary_conditions_list(i) == 1
       number_of_fixed_dof = number_of_fixed_dof + 1;
    end
end
fixed_dof = zeros(number_of_fixed_dof,1);
dof_per_node = 3;
counter = 1;

% build the list of fixed dof
for bc = 1:size(boundary_conditions,1)
  if boundary_conditions(bc,2) == 1
      fixed_dof(counter) = boundary_conditions(bc,1)*dof_per_node-2;
      counter = 1 + counter;
  end
  if boundary_conditions(bc,3) == 1
      fixed_dof(counter) = boundary_conditions(bc,1)*dof_per_node-1;
      counter = 1 + counter;
  end
  if boundary_conditions(bc,4) == 1
      fixed_dof(counter) = boundary_conditions(bc,1)*dof_per_node;
      counter = 1 + counter;
  end
end

% determine the free dofs
fixed_dof = sort(fixed_dof);
total_dof = dof_per_node * size(nodes,1);
all_dof = 1:total_dof;
free_dof = zeros((total_dof-size(fixed_dof,1)),1);
free_dof = reshape(all_dof(~ismember(all_dof,fixed_dof)),size(free_dof,1),1);

h = 0.01;
gamma = 0.5;
beta = 0.25;

time_and_motion_data = fopen('data/single_column/time_and_motion_data.csv','w');
fprintf(time_and_motion_data,'delta_t,t,u, udot, udotdot\r\n');
udotdot = 0;
udot = 0.0;
u = 0.0;

% BEGIN LOOP
for delta_t = 1:(size(acceleration,2)-1)
fprintf(time_and_motion_data,'%g, %g, %g, %g, %g \r\n',delta_t, ((delta_t-1)*h), u, udot, udotdot);

% create a blank global stiffness matrix
Kg = zeros(total_dof);

% assemble the element stiffness matricies and add to the global stiffness
% matrix
for e = 1:size(elements,1)
    [dof] = getElementDegreesOfFreedom(e,elements,dof_per_node);
    [L,theta] = getElementLengthAndAngle(e,elements,nodes);
    [T,Tt] = getTransformationMatrix(theta);
    Ke_dash = getElementStiffnessMatrix(e,A,E,I,L);
    Ke = Tt * Ke_dash * T;
    Kg(dof,dof) = Kg(dof, dof) + Ke;                
end

% load in mass matrix from elements
M = zeros(total_dof);
Me = zeros(dof_per_node *2);
for e = 1:size(elements)
    [dof] = getElementDegreesOfFreedom(e,elements,dof_per_node);
    [L, theta] = getElementLengthAndAngle(e,elements,nodes);
    [T,Tt] = getTransformationMatrix(theta);
    % element stiffness matrix
    Me(1,1) = A(e) * L * rho(e) * 0.5 / 1000; 
%     Me(2,2) = A(e) * L * rho(e) * 0.5 / 1000;
%     Me(3,3) = A(e) * L * rho(e) * 0.5 / 1000;
    Me(4,4) = A(e) * L * rho(e) * 0.5 / 1000;
%     Me(5,5) = A(e) * L * rho(e) * 0.5 / 1000;
%     Me(6,6) = A(e) * L * rho(e) * 0.5 / 1000;
    Me = Tt * Me * T;
    M(dof,dof) = M(dof,dof) + Me;
end

% Damping matrix
am = 0.05;
ak = 0.05;
C = M.*0.05 + Kg.*0.05;

% Effective Stiffness Matrix
Kge = Kg +(2/h)*C + (4/h^2)*M;

% remove rows and columns for fixed dof
Kge = Kge(~ismember(1:size(Kge,1),fixed_dof),~ismember(1:size(Kge,1),fixed_dof));

forces=zeros(size(nodes,1)*dof_per_node,1);

% load in forces from nodes
for force=1:size(nodes,1)
    forces(force*dof_per_node-2,1)=nodes(force,5);
    forces(force*dof_per_node-1,1)=nodes(force,6);
    forces(force*dof_per_node,1)=nodes(force,7);
end

% load in forces from elements
for e = 1:size(elements)
    [dof] = getElementDegreesOfFreedom(e,elements,dof_per_node);
    forces(dof) = forces(dof) + reshape(elements(e,4:9),2*dof_per_node,1);
end

% assemble the force matrix
P = forces + C.*(2/h*u + udot)+ M.*(4/(h^2)*u + 4/h*udot + udotdot);

% remove rows and columns for fixed dof
analysis_forces=P(~ismember(1:size(forces,1),fixed_dof),1);

% solve the displacements
displacements = Kge\analysis_forces;
node_displacements = zeros(total_dof,1);
node_displacements(free_dof) = node_displacements(free_dof) + displacements;

% transform displacements for elements into local
element_displacements = zeros(size(elements,1),2*dof_per_node); %x1, y1, z1, x2, y2, z2 displacments
element_local_displacements = zeros(size(elements,1),2*dof_per_node);
element_loads = zeros(size(elements,1),2*dof_per_node);

% add local displacements to elements
for e = 1:size(elements,1)
    [dof] = getElementDegreesOfFreedom(e,elements,dof_per_node);
    [L,theta] = getElementLengthAndAngle(e,elements,nodes);
    [T,Tt] = getTransformationMatrix(theta);
    Ke_dash = getElementStiffnessMatrix(e,A,E,I,L);
    ED_column = node_displacements(dof);
    ED_row = reshape(ED_column,1,size(ED_column,1));
    element_displacements(e,1:2*dof_per_node) = ED_row;
    ED_dash_column = T*ED_column;
    ED_dash_row = reshape(ED_dash_column,1,size(ED_dash_column,1));
    element_local_displacements(e,1:2*dof_per_node) = ED_dash_row;
    element_loads(e,1:2*dof_per_node) = Ke_dash*ED_dash_column - reshape(elements(e,4:9),2*dof_per_node,1);
    load_conversion_matrix = [1 1 -1 1 -1 1];
    element_loads(e,1:2*dof_per_node) = element_loads(e,1:2*dof_per_node).*load_conversion_matrix;
end

temp_nodes = nodes;
for n = 1:size(nodes,1)
    dof = (n*dof_per_node -2):(n*dof_per_node);
    temp_nodes(n,2:4) = nodes(n,2:4) + reshape(node_displacements(dof),1,dof_per_node);
end

for n = 1:size(nodes,1)
    output = temp_nodes(n,2:4);
    fprintf(nodes_out,'%g, %g, %g, %g, %g \r\n',delta_t, n, output);
    fprintf('\r\n');
end

for e = 1:size(element_loads)
    output = element_loads(e,1:end);
    fprintf(elements_out,'%g, %g, %g, %g, %g, %g, %g, %g \r\n',delta_t, e, output);
    fprintf('\r\n');
end

% update motion data
udotdotn = udotdot;
udotdot = acceleration(delta_t+1);
udotn = udot;
udot = udot + (1-gamma)*h*udotdotn + gamma * h * udotdot;
u = u + h*udotn + ((1-2*beta)/2)*h^2*udotdotn + beta * h^2*udotdot;

end
fclose('all');

% FUNCTION SPACE

% Returns the transformation matrix for a given angle
function [T,Tt] = getTransformationMatrix(theta)
    c = cos(theta);
    s = sin(theta);
    Tt = [  c -s 0 0 0 0;
            s c 0 0 0 0;
            0 0 1 0 0 0;
            0 0 0 c -s 0;
            0 0 0 s c 0;
            0 0 0 0 0 1];
    T = transpose(Tt);
end

% Returns the Length and Angle of a given element
function [L,theta] = getElementLengthAndAngle(i,elements,nodes)
    start_node = elements(i, 2);
    x1 = nodes(start_node,2);
    y1 = nodes(start_node,3);
 
    end_node = elements(i,3);
    x2 = nodes(end_node,2);
    y2 = nodes(end_node,3);
    
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    theta = atan((y2-y1)/(x2-x1));
end

% Returns the indicies of the degrees of freedom for a given element
function [dof] = getElementDegreesOfFreedom(i,elements,dof_per_node)
    start_node = elements(i, 2);
    start_node_dof = start_node * [dof_per_node,dof_per_node,dof_per_node] - [2, 1, 0];
    end_node = elements(i,3);
    end_node_dof = end_node * [dof_per_node,dof_per_node,dof_per_node] - [2, 1, 0];
    dof = [start_node_dof, end_node_dof];
end

% Returns the local element stiffness matrix for a given element
function [Ke_dash] = getElementStiffnessMatrix(e,A,E,I,L)
    Ke_dash = [
        A(e)*E(e)/L 0 0 -A(e)*E(e)/L 0 0;
        0 12*E(e)*I(e)/L^3 6*E(e)*I(e)/L^2 0 -12*E(e)*I(e)/L^3 6*E(e)*I(e)/L^2;
        0 6*E(e)*I(e)/L^2 4*E(e)*I(e)/L 0 -6*E(e)*I(e)/L^2 2*E(e)*I(e)/L;
        -A(e)*E(e)/L 0 0 A(e)*E(e)/L 0 0;
        0 -12*E(e)*I(e)/L^3 -6*E(e)*I(e)/L^2 0 12*E(e)*I(e)/L^3 -6*E(e)*I(e)/L^2;
        0 6*E(e)*I(e)/L^2 2*E(e)*I(e)/L 0 -6*E(e)*I(e)/L^2 4*E(e)*I(e)/L];
end