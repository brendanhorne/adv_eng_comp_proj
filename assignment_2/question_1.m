% Assignment 1
clc
clearvars
%% Geometry Information
%Node information: Node No., x-coord, x-force,
node_infor=[1, 0, 0;
            2, 2, 0;
            3, 4, 10000;
            4, 6, 0;
            5, 8, -20000;
            6, 10, 0];
        
%Element information (Connectivity table)
%Element no, starting node, ending node
element_infor=[1, 1, 2;
               2, 2, 3;
               3, 3, 4;
               4, 4, 5;
               5, 5, 6];
           
%Boundary condition
%BC=[Node no, x-direction movement, y-direction movement] 1-fix, 0-free
BC=[1, 1];
fix_dof=[];

%size(matrix/array, dimension [row=1,column=2])
for pp=1:size(BC,1) %load up array with dof ids to remove before building stiffness matrix
  if BC(pp,2)==1
      fix_dof(pp,1)=BC(pp,1)*2-1;
  else
      fix_dof(pp,1)=0;
  end
  
end
fixx=reshape(fix_dof,[size(fix_dof,1)*size(fix_dof,2),1]); %reshape the array into a column
fix=sort(fixx(find(fixx~=0))); %find unique non zero elements and sort in ascending order

%% Material Information
E=200e9;
A=0.01;

%% Greate the global stiffness matrix 
node_dof=1; %DOF per node
DOF=size(node_infor, 1)*node_dof;
KK=zeros(DOF, DOF);

for ii=1:size(element_infor,1)
    node1=element_infor(ii,2);
    node2=element_infor(ii,3);
    LL(ii,1)=sqrt((node_infor(node1,2)-node_infor(node2,2))^2);
    ke = E*A/LL(ii,1);
    KK_e = ke*[1, -1; -1, 1];
%   loc_vec=[node1*2-1, node1*2, node2*2-1, node2*2];
    loc_vec=[node1, node2];
    KK(loc_vec, loc_vec)=KK(loc_vec, loc_vec)+KK_e;
    
end
%Global stiffness matrix without fixed DOF
Kt=KK(~ismember(1:size(KK,1),fix),~ismember(1:size(KK,1),fix));

%% Global force vector
FF=zeros(size(node_infor,1),1);
for hh=1:size(node_infor,1)
    FF(hh,1)=node_infor(hh,3);
end
Ft=FF(~ismember(1:size(FF,1),fix),1);

%% Displacement
Ut=Kt\Ft 
    