clc, clear, close
nodes = [% node ID, x-coordinate, y-coordinate, z-rotation, x-force, y-force, z-moment
            1 0 0 0 0 0 0; 
            2 3 0 0 0 0 75;
            3 4.5 -2 0 0 0 0];

boundary_conditions = [% node ID, x-position-fixity, y-position-fixity, z-rotation-fixity
            1 1 1 1
            3 1 1 1];

elements = [% element ID, start node, end node, x-force, y-force, z-moment, E, I, A
            1, 1, 2 0 -3 -1.5 0 -3 1.5 290e6 0.055 0.16;
            2, 2, 3 0 0 0 0 0 0 290e6 0.028 0.12];

E = elements(1:end,10);
I = elements(1:end,11);
A = elements(1:end,12);

boundary_conditions_only = boundary_conditions(1:end,2:end);
number_possible_fixed_dof = numel(boundary_conditions_only);
boundary_conditions_list = reshape(boundary_conditions_only,number_possible_fixed_dof,1);
fixed_only_boundary_condition_list = ismember(ones(number_possible_fixed_dof,1),boundary_conditions_list);
number_of_fixed_dof = numel(fixed_only_boundary_condition_list);
fixed_dof = zeros(number_of_fixed_dof,1);
dof_per_node = 3;
counter = 1;

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

fixed_dof = sort(fixed_dof);
total_dof = dof_per_node * size(nodes,1);
all_dof = 1:total_dof;
free_dof = zeros((total_dof-size(fixed_dof,1)),1);
free_dof = reshape(all_dof(~ismember(all_dof,fixed_dof)),size(free_dof,1),1);

Kg = zeros(total_dof);

for e = 1:size(elements,1)
    [dof] = getElementDegreesOfFreedom(e,elements,dof_per_node);
    [L,theta] = getElementLengthAndAngle(e,elements,nodes);
    [T,Tt] = getTransformationMatrix(theta);
    Ke_dash = getElementStiffnessMatrix(e,A,E,I,L);
    Ke = Tt * Ke_dash * T;
    Kg(dof,dof) = Kg(dof, dof) + Ke;                
end

Kg = Kg(~ismember(1:size(Kg,1),fixed_dof),~ismember(1:size(Kg,1),fixed_dof));

forces=zeros(size(nodes,1)*dof_per_node,1);

for force=1:size(nodes,1)
    forces(force*dof_per_node-2,1)=nodes(force,5);
    forces(force*dof_per_node-1,1)=nodes(force,6);
    forces(force*dof_per_node,1)=nodes(force,7);
end

for e = 1:size(elements)
    [dof] = getElementDegreesOfFreedom(e,elements,dof_per_node);
    forces(dof) = forces(dof) + reshape(elements(e,4:9),2*dof_per_node,1);
end

analysis_forces=forces(~ismember(1:size(forces,1),fixed_dof),1);

displacements = Kg\analysis_forces;
node_displacements = zeros(total_dof,1);
node_displacements(free_dof) = node_displacements(free_dof) + displacements;

element_displacements = zeros(size(elements,1),2*dof_per_node); %x1, y1, z1, x2, y2, z2 displacments
element_local_displacements = zeros(size(elements,1),2*dof_per_node);
element_loads = zeros(size(elements,1),2*dof_per_node);

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
