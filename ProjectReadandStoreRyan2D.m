clear, clc, close

%% Read and Store for 2D frame

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
    end
    
%Elements
%elements = [ elements no. , starting node, ending node, element type]

epf = 9; %elements per floor
bpf = 4; %Beams per floor
cpf = 5; %Columns per floor
bsi = cpf + 1; % beam starting index
n = 0; %floor number

elements = zeros(total_members,4);
    for i = 1:total_members
        elements(i,1) = i;
        n = i * 1/9;
        elements(i,2) = i - floor((i - 1)/9)*4;
        elements(i,3) = i + 5 - floor((i + 3)/9)*4;
        if  elements (i, 3) - elements(i,2) < 5
        elements(i,4) = 2;
        elseif rem((i-3),9) == 0
            elements(i,4) = 3;
        else
            elements(i,4) = 1;
        end
    end
    
    E = 35e9;
    ro = 3000;
    k = 4e9;
    Ic = (0.4^4)/12;
    yb = (0.25*0.45*0.225 + 0.6*0.15*0.525) /  (0.25*0.45 + 0.6*0.15);
    Ib =  0.25*0.45*(yb - 0.45/2)^2 + (0.25*0.45^3)/12 + 0.6*0.15*(0.45 + 0.15/2 - yb)^2 + (0.6*0.15^3)/12;
    Is = (02.5*4.5^3)/12;
    
    %Section Properties = [ section number, E, I, A]; (1 = column, 2= beam, 3 = shear
    %wall) 
    
    sec_props = [1, 1, Ic, 0.16;
                 2, 1, Ib, 0.2025;
                 3, 1, Is, 1];
                

    csvwrite('nodedata2D.csv',nodes)
    csvwrite('elementdata2D.csv',elements)
    csvwrite('sectiondata2D.csv',sec_props)
    
    type('nodedata2D.csv')
    type('elementdata2D.csv')
    type('sectiondata2D.csv')


    