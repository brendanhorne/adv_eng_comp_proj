clc, clear, close
% All units in SI kg, m, s, N, Pa etc.

% INPUT SPACE
out = fopen('test.txt','w');
nodes = [% node ID, x-coordinate, y-coordinate, z-rotation, x-force, y-force, z-moment
            1 0 0 0 0 0 0;
            2 0 1 0 0 0 0];

fprintf('NODES \r\n');
fprintf('%-12s %-12s %-12s %-12s %-12s %-12s %-12s \r\n','Node id','x-cord(m)', 'y-cord(m)','z-rot(rad)','x-force(N)','y-force(N)','z-moment(Nm)');
fprintf('%-12d %-12d %-12d %-12d %-12d %-12d %-12d \r\n',transpose(nodes));
fprintf('\r\n');

elements = [% element ID, start node, end node, x1-force, y1-force, z1-moment, x2-force, y2-force, z2-moment, E, I, A, rho
            1 1 2 0 0 0 0 0 0 200e9 8.33e-6 0.01, 3e3];

fprintf('ELEMENTS \r\n');
fprintf('%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s \r\n','Element id','start node', 'end node','N1(N)','V1(N)','M1(Nm)','N2(N)','V2(N)','M2(Nm)','E(Pa)','I(m^4)','A(m^2)');
fprintf('%-12d %-12d %-12d %-12d %-12d %-12d %-12d %-12d %-12d %-12.3g %-12.3g %-12.3g \r\n',transpose(elements));
fprintf('\r\n');

boundary_conditions = [% node ID, x-position-fixity, y-position-fixity, z-rotation-fixity
            1 1 1 1];

fprintf('BOUNDARY CONDITIONS \r\n');
fprintf('0-free, 1-fixed \r\n');
fprintf('%-12s %-12s %-12s %-12s \r\n','Node id','x-disp', 'y-disp','z-rot');
fprintf('%-12d %-12d %-12d %-12d \r\n',transpose(boundary_conditions));
fprintf('\r\n');
        
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

% SELF WEIGHT MATRIX
sw = [0 -1 0];

% determine the free dofs
fixed_dof = sort(fixed_dof);
total_dof = dof_per_node * size(nodes,1);
all_dof = 1:total_dof;
free_dof = zeros((total_dof-size(fixed_dof,1)),1);
free_dof = reshape(all_dof(~ismember(all_dof,fixed_dof)),size(free_dof,1),1);

% create a blank global stiffness matrix
Kg = zeros(total_dof);

% assemble the element stiffness matricies and add to the global stiffness
% matrix
for e = 1:size(elements,1)
    [dof] = getElementDegreesOfFreedom(e,elements,dof_per_node);
    [L,theta] = getElementLengthAndAngle(e,elements,nodes);
    [T,Tt] = getTransformationMatrix(theta);
    Ke_dash = getElementStiffnessMatrix(e,A,E,I,L);
    fprintf(['ELEMENT #','%d',' STIFFNESS MATRIX: \r\n'],e);
    fprintf('%-8.3g %-8.3g %-8.3g %-8.3g %-8.3g %-8.3g \r\n',transpose(Ke_dash));
    fprintf('\r\n');
    Ke = Tt * Ke_dash * T;
    Kg(dof,dof) = Kg(dof, dof) + Ke;                
end


fprintf('GLOBAL STIFFNESS MATRIX: \r\n');
for i = 1:total_dof
   for j = 1:total_dof
      fprintf('%-8.3g ',Kg(j,i)); % rather than transpose, swap i&j
   end
   fprintf('\r\n');
end
fprintf('\r\n');

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

fprintf('FORCE VECTOR: \r\n');
fprintf('%-8s %-8s %-8s \r\n','x(N)','y(N)','z(Nm)');
output = reshape(transpose(forces),(total_dof/3),3);
fprintf('%-8.3g %-8.3g %-8.3g \r\n',output);
fprintf('\r\n');

% load in mass matrix from elements
SW = [sw, sw];
M = zeros(total_dof);
Me = zeros(dof_per_node *2);
for e = 1:size(elements)
    [dof] = getElementDegreesOfFreedom(e,elements,dof_per_node);
    [L, theta] = getElementLengthAndAngle(e,elements,nodes);
    for d = 1:(dof_per_node*2)
        Me(d,d) = A(e) * L * rho(e) * 0.5;
    end    
    M(dof,dof) = M(dof,dof) + Me;
end

h = 0.01;
udotdot = 1;
udot = udotdot*h;
u = udot * h;

fprintf('MASS MATRIX: \r\n');
output = transpose(M);
fprintf('%-8.3g %-8.3g %-8.3g %-8.3g %-8.3g %-8.3g \r\n',output);
fprintf('\r\n');

% Damping matrix
am = 0.05;
ak = 0.05;
C = M.*0.05 + Kg.*0.05

% assemble the force matrix
P = forces + C.*(2/h*u + udot)+ M.*(4/(h^2)*u + 4/h*udot + udotdot)

% remove rows and columns for fixed dof
Kg = Kg(~ismember(1:size(Kg,1),fixed_dof),~ismember(1:size(Kg,1),fixed_dof));

% remove rows and columns for fixed dof
analysis_forces=P(~ismember(1:size(forces,1),fixed_dof),1);

% solve the displacements
displacements = Kg\analysis_forces;
node_displacements = zeros(total_dof,1);
node_displacements(free_dof) = node_displacements(free_dof) + displacements;

fprintf('NODE DISPLACEMENTS: \r\n');
fprintf('%-8s %-8s %-8s \r\n','x(m)','y(m)','z(rad)');
output = reshape(transpose(node_displacements),(total_dof/3),3);
fprintf('%-8.3g %-8.3g %-8.3g \r\n',output);
fprintf('\r\n');

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
element_loads
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