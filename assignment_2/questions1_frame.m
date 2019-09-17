clc, clear, close
nodes = [% node ID, x-coordinate, y-coordinate, z-rotation, x-force, y-force, z-moment
            1 0 0 0 0 0 0; 
            2 2 0 0 0 0 0;
            3 4 0 0 10000 0 0;
            4 6 0 0 0 0 0;
            5 8 0 0 -20000 0 0;
            6 10 0 0 0 0 0];

boundary_conditions = [% node ID, x-position-fixity, y-position-fixity, z-rotation-fixity
            1 1 1 1];
        
elements = [% element ID, start node, end node, E, I, A
            1, 1, 2 200e9 0.1 0.01;
            2, 2, 3 200e9 0.1 0.01;
            3, 3, 4 200e9 0.1 0.01;
            4, 4, 5 200e9 0.1 0.01;
            5, 5, 6 200e9 0.1 0.01];

E = elements(1:end,4);
I = elements(1:end,5);
A = elements(1:end,6);
        
fixed_dofs = [];
dof_per_node = 3;

for bc = 1:size(boundary_conditions,1)
  if boundary_conditions(bc,2) == 1
      fixed_dofs(end+1,1) = boundary_conditions(bc,1)*dof_per_node-2;
  end
  if boundary_conditions(bc,3) == 1
      fixed_dofs(end+1,1) = boundary_conditions(bc,1)*dof_per_node-1;
  end
    if boundary_conditions(bc,4) == 1
      fixed_dofs(end+1,1) = boundary_conditions(bc,1)*dof_per_node;
  end
end
fixed_dofs = sort(fixed_dofs);
total_dof = dof_per_node * size(nodes,1);
Kg = zeros(total_dof);
for e = 1:size(elements,1)
    start_node = elements(e, 2);
    x1 = nodes(start_node,2);
    y1 = nodes(start_node,3);
    start_node_dofs = start_node * [dof_per_node,dof_per_node,dof_per_node] - [2, 1, 0];
    end_node = elements(e,3);
    x2 = nodes(end_node,2);
    y2 = nodes(end_node,3);
    end_node_dofs = end_node * [dof_per_node,dof_per_node,dof_per_node] - [2, 1, 0];
    dofs = [start_node_dofs, end_node_dofs];
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    theta = atan((y2-y1)/(x2-x1));
    c = cos(theta);
    s = sin(theta);
    Tt = [
            c -s 0 0 0 0;
            s c 0 0 0 0;
            0 0 1 0 0 0;
            0 0 0 c -s 0;
            0 0 0 s c 0;
            0 0 0 0 0 1];
    T = transpose(Tt);
    Ke_dash = [
            A(e)*E(e)/L 0 0 -A(e)*E(e)/L 0 0;
            0 12*E(e)*I(e)/L^3 6*E(e)*I(e)/L^2 0 -12*E(e)*I(e)/L^3 6*E(e)*I(e)/L^2;
            0 6*E(e)*I(e)/L^2 4*E(e)*I(e)/L 0 -6*E(e)*I(e)/L^2 2*E(e)*I(e)/L;
            -A(e)*E(e)/L 0 0 A(e)*E(e)/L 0 0;
            0 -12*E(e)*I(e)/L^3 -6*E(e)*I(e)/L^2 0 12*E(e)*I(e)/L^3 -6*E(e)*I(e)/L^2;
            0 6*E(e)*I(e)/L^2 2*E(e)*I(e)/L 0 -6*E(e)*I(e)/L^2 4*E(e)*I(e)/L];
    Ke = Tt * Ke_dash * T;
    Kg(dofs,dofs) = Kg(dofs, dofs) + Ke;                
end
Kg = Kg(~ismember(1:size(Kg,1),fixed_dofs),~ismember(1:size(Kg,1),fixed_dofs));

forces=zeros(size(nodes,1)*dof_per_node,1);
for force=1:size(nodes,1)
    forces(force*dof_per_node-2,1)=nodes(force,5);
    forces(force*dof_per_node-1,1)=nodes(force,6);
    forces(force*dof_per_node,1)=nodes(force,7);
end
forces=forces(~ismember(1:size(forces,1),fixed_dofs),1);

displacements=Kg\forces