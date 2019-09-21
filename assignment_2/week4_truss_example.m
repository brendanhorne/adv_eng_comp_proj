clc
clearvars
%% Geometry Information
%Node information: Node No.  x-coord, y-coord, x-force, y-force
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
for pp=1:size(BC,1)
  if BC(pp,2)==1
      fix_dof(pp,1)=BC(pp,1)*2-1;
  else
      fix_dof(pp,1)=0;
  end
  
  if BC(pp,3)==1
      fix_dof(pp,2)=BC(pp,1)*2;
  else
      fix_dof(pp,2)=0;
  end
end
fixx=reshape(fix_dof,[size(fix_dof,1)*size(fix_dof,2),1]);
fix=sort(fixx(find(fixx~=0)));

%% Material Information
E=200e9;
A=0.01;

%% Greate the global stiffness matrix 
node_dof=2; %DOF per node
DOF=size(node_infor, 1)*node_dof;
KK=zeros(DOF, DOF);

for ii=1:size(element_infor,1)
    node1=element_infor(ii,2);
    node2=element_infor(ii,3);
    LL(ii,1)=sqrt((node_infor(node1,3)-node_infor(node2,3))^2+(node_infor(node1,2)-node_infor(node2,2))^2);
    
    theta(ii,1)=atan((node_infor(node1,3)-node_infor(node2,3))/(node_infor(node1,2)-node_infor(node2,2)));
    
    ll=cos(theta(ii,1));
    mm=sin(theta(ii,1));
    
    KK_e=(E*A/LL(ii,1))*[ll^2, ll*mm, -ll^2, -ll*mm;
                         ll*mm, mm^2, -ll*mm, -mm^2;
                         -ll^2, -ll*mm, ll^2, ll*mm;
                        -ll*mm, -mm^2, ll*mm, mm^2];
    loc_vec=[node1*2-1, node1*2, node2*2-1, node2*2];
    KK(loc_vec, loc_vec)=KK(loc_vec, loc_vec)+KK_e;
    
end
%Global stiffness matrix without fixed DOF
Kt=KK(~ismember(1:size(KK,1),fix),~ismember(1:size(KK,1),fix));

%% Global force vector
FF=zeros(size(node_infor,1)*2,1);
for hh=1:size(node_infor,1)
    FF(hh*2-1,1)=node_infor(hh,4);
    FF(hh*2,1)=node_infor(hh,5);
end
Ft=FF(~ismember(1:size(FF,1),fix),1);

%% Displacement
Ut=Kt\Ft
    
