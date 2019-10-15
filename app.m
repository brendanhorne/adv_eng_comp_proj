clc, clear, close
% All units in SI kg, m, s, N, Pa etc.

% READ INPUT

% 2D Frame
% model_2D_8_story_frame
% Northridge_earthquake
% dof_of_interest = 45*3-1;
% element_of_interest = 1;
% load_of_interest = 1;
% h = 0.02;
% Time = 0:h:((size(load,2)-2)*h);
% transient_direction = [0 1 0];

% Single Column
single_column
calfem_test_data
dof_of_interest = 4;
element_of_interest = 1;
h = 0.001;
Time = 0:h:((size(load,2)-2)*h);
transient_direction = [1 0 0];

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

% TODO add spring stiffness to DOF appropriate.

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
total_dof = dof_per_node * size(nodes,1);%add spring
all_dof = 1:total_dof;
free_dof = zeros((total_dof-size(fixed_dof,1)),1);
free_dof = reshape(all_dof(~ismember(all_dof,fixed_dof)),size(free_dof,1),1);

gamma = 0.5;
beta = 0.25;

% BETA VALUES
b1 = h*h*0.5*(1-2*beta);
b2 = (1-gamma)*h;
b3 = gamma*h;
b4 = beta*h*h;

% Transient Load Vector
transient_vector = zeros(total_dof,1);
for d = 1:total_dof 
   r = rem(d,3);
   if r == 1
       transient_vector(d) = transient_direction(1);
   end
   if r == 2
       transient_vector(d) = transient_direction(2);
   end
   if r == 0
       transient_vector(d) = transient_direction(3);
   end
end

% create a blank global stiffness matrix
Kg = zeros(total_dof);
M = zeros(total_dof);
C = zeros(total_dof);
dof_mass = zeros(total_dof,1);

% Damping matrix
am = 0.05;
ak = 0.05;

% assemble the element stiffness matricies and add to the global stiffness
for e = 1:size(elements,1)
    [dof] = getElementDegreesOfFreedom(e,elements,dof_per_node);
    [L,theta] = getElementLengthAndAngle(e,elements,nodes);
    [T,Tt] = getTransformationMatrix(theta);
    Ke_dash = getElementStiffnessMatrix(e,A,E,I,L);    
    Ke = Tt * Ke_dash * T;
    Kg(dof,dof) = Kg(dof, dof) + Ke;     
    Me_dash = getMassMatrix(e,A,L,rho);
    Me = Tt * Me_dash * T;
    M(dof,dof) = M(dof, dof) + Me;
    Ce_dash = am.*Me_dash + ak.*Ke_dash;
    Ce = Tt * Ce_dash * T;
    C(dof,dof) = C(dof, dof) + Ce;  
    dof_mass(dof,1) = dof_mass(dof,1) + ones(dof_per_node*2,1).*L*A(e)*rho(e)/2;
end

% spring element
for s = 1:size(springs)
    Kspring = [springs(s,4) 0 0 -springs(s,4) 0 0;
                0 springs(s,5) 0 0 -springs(s,5) 0;
                0 0 springs(s,6) 0 0  -springs(s,6);
                -springs(s,4) 0 0 springs(s,4) 0 0;
                0 -springs(s,5) 0 0 springs(s,5) 0;
                0 0 -springs(s,6) 0 0 springs(s,6)];
    [dof] = geSpringDegreesOfFreedom(s,springs,dof_per_node);
    Kg(dof,dof) = Kg(dof,dof) + Kspring;
    C(dof,dof) = C(dof,dof) + ak.*Kspring;
end

% Effective Stiffness Matrix
Kge = M+b3*C+b4*Kg;

% remove rows and columns for fixed dof
Kge = Kge(~ismember(1:size(Kge,1),fixed_dof),~ismember(1:size(Kge,1),fixed_dof));

% Assembe load vector
transient_vector = dof_mass.*transient_vector;

dnew = zeros(total_dof,1);
vnew = zeros(total_dof,1);
anew = zeros(total_dof,1);
dt = h;
dold=dnew;      vold=vnew;      aold=anew;

% Initialise the results variables
node_results = zeros(total_dof, size(load,2)-1);
element_results = cell(size(elements,1),1);
for e = 1:size(element_results,1)
   element_results{e,1} = zeros(2*dof_per_node,size(load,2)-1);
end

% BEGIN LOOP
for delta_t = 1:(size(load,2)-1)
time = delta_t * h;  
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
dold=dnew;      vold=vnew;      aold=anew;
dpred=dold+dt*vold+b1*aold;
vpred=vold+b2*aold;

% assemble the force matrix
P= transient_vector.*load(delta_t+1)-C*vpred-Kg*dpred;

% remove rows and columns for fixed dof
P=P(~ismember(1:total_dof,fixed_dof),1);
dpred=dpred(~ismember(1:total_dof,fixed_dof),1);
vpred=vpred(~ismember(1:total_dof,fixed_dof),1);

% solve the displacements
anew=Kge\P;
dnew=dpred+b4*anew;  
vnew=vpred+b3*anew;
dold=dnew;      vold=vnew;      aold=anew;
anew = zeros(total_dof,1);
anew(free_dof) = anew(free_dof) + aold;
vnew = zeros(total_dof,1);
vnew(free_dof) = vnew(free_dof) + vold;
displacements = dnew;
node_displacements = zeros(total_dof,1);
node_displacements(free_dof) = node_displacements(free_dof) + displacements;
dnew = zeros(total_dof,1);
dnew(free_dof) = dnew(free_dof) + dold;

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
    Ke = Tt*Ke_dash*T;
    element_loads(e,1:2*dof_per_node) = Ke_dash*node_displacements(dof,1);
%     ED_column = node_displacements(dof);
%     ED_row = reshape(ED_column,1,size(ED_column,1));
%     element_displacements(e,1:2*dof_per_node) = ED_row;
%     ED_dash_column = T*ED_column;
%     ED_dash_row = reshape(ED_dash_column,1,size(ED_dash_column,1));
%     element_local_displacements(e,1:2*dof_per_node) = ED_dash_row;
%     element_loads(e,1:2*dof_per_node) = Ke_dash*ED_dash_column - reshape(elements(e,4:9),2*dof_per_node,1);
%     load_conversion_matrix = [1 1 -1 1 -1 1];
%     element_loads(e,1:2*dof_per_node) = element_loads(e,1:2*dof_per_node).*load_conversion_matrix;
    element_results{e,1}(:,delta_t) = element_loads(e,:);
end
node_results(:,delta_t) = node_displacements;

end

% figure
% plot(Time,node_results(dof_of_interest,:))
% title(['Displacement at dof ',num2str(dof_of_interest)]);
% xlabel('Time (s)');
% ylabel('Displacement (m)');

figure
subplot(2,3,1);
plot(Time,element_results{element_of_interest,1}(1,:))
title('N1')
xlabel('t')
ylabel('N')

subplot(2,3,2);
plot(Time,element_results{element_of_interest,1}(2,:))
title('V1')
xlabel('t')
ylabel('N')

subplot(2,3,3);
plot(Time,element_results{element_of_interest,1}(3,:))
title('M1')
xlabel('t')
ylabel('Nm')

subplot(2,3,4);
plot(Time,element_results{element_of_interest,1}(4,:))
title('N2')
xlabel('t')
ylabel('N')

subplot(2,3,5);
plot(Time,element_results{element_of_interest,1}(5,:))
title('V2')
xlabel('t')
ylabel('N')

subplot(2,3,6);
plot(Time,element_results{element_of_interest,1}(6,:))
title('M2')
xlabel('t')
ylabel('Nm')

% FUNCTION SPACE

% Returns the transformation matrix for a given angle
function [T,Tt] = getTransformationMatrix(theta)
    c = cos(theta);
    s = sin(theta);
    T = [  c s 0  0 0 0;
           -s c 0  0 0 0;
            0 0 1  0 0 0;
            0 0 0  c s 0;
            0 0 0 -s c 0;
            0 0 0  0 0 1];
    Tt = transpose(T);
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

% Returns the indicies of the degrees of freedom for a given spring
function [dof] = geSpringDegreesOfFreedom(i,springs,dof_per_node)
    start_node = springs(i, 2);
    start_node_dof = start_node * [dof_per_node,dof_per_node,dof_per_node] - [2, 1, 0];
    end_node = springs(i,3);
    end_node_dof = end_node * [dof_per_node,dof_per_node,dof_per_node] - [2, 1, 0];
    dof = [start_node_dof, end_node_dof];
end

% Returns the local element stiffness matrix for a given element
function [Ke_dash] = getElementStiffnessMatrix(e,A,E,I,L)
    Ke_dash = [
        A(e)*E(e)/L  0                 0               -A(e)*E(e)/L  0                 0;
        0            12*E(e)*I(e)/L^3  6*E(e)*I(e)/L^2  0           -12*E(e)*I(e)/L^3  6*E(e)*I(e)/L^2;
        0            6*E(e)*I(e)/L^2   4*E(e)*I(e)/L    0           -6*E(e)*I(e)/L^2   2*E(e)*I(e)/L;
        -A(e)*E(e)/L 0                 0                A(e)*E(e)/L  0                 0;
        0           -12*E(e)*I(e)/L^3 -6*E(e)*I(e)/L^2  0            12*E(e)*I(e)/L^3 -6*E(e)*I(e)/L^2;
        0            6*E(e)*I(e)/L^2   2*E(e)*I(e)/L    0           -6*E(e)*I(e)/L^2   4*E(e)*I(e)/L];
end

% get the mass matrix
function [Me_dash] = getMassMatrix(e,A,L,rho)
    m = A(e) * L * rho(e);
    Me_dash = m/420.* [   140 0     0     70  0     0;
                          0   156   22*L  0   54   -13*L;
                          0   22*L  4*L^2 0   13*L -3*L^2;
                          70  0     0     140 0     0;
                          0   54    13*L  0   156  -22*L;
                          0  -13*L -3*L^2 0  -22*L  4*L^2];

end