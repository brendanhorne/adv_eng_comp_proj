clear, clc, close

total_nodes = 135;
total_members = 296;

%Nodes
%nodes = [node no. , x-coord, y-coord, z-coord]

nodes = zeros(135,4);
    for i = 1:135
        %Node Number
        nodes(i,1) = i;
        
        %x coordinate
        nodes(i,2) = (rem((i-1),5)) * 5;
        
        %y coordinates
        nodes(i,3) = floor((i-1)/15)*3.5;
        
        %z coordinates
        nodes(i,4) = floor(rem(i-1,15)/5)*6;
        

    end
    
%Elements
%elements = [ elements no. , starting node, ending node, element type]

epf = 37; %elements per floor
bpf = 22; %Beams per floor
cpf = 15; %Columns per floor
bsi = cpf + 1; % beam starting index
n = 0; %floor number

elements = zeros(296,4);
    for i = 1:296
        elements(i,1) = i;
        if rem((i-1),epf) == 0
            for j = 1:epf
                if j <= cpf
                    elements((i+j-1),2) = j+n*cpf;
                    elements((i+j-1),3) = elements((i+j-1),2) + cpf;
                    if rem(j,5) == 3
                        elements((i+j-1),4) = 3;
                    else
                        elements((i+j-1),4) = 1;
                    end
                    
                elseif  j > cpf
                    k = (j + n*epf);
                    m = (k - (bsi + n*epf))/9;
                    elements(k,2) = (k-n*bpf)-floor(m+(5/9))*4;
                    elements(k,3) = (k-n*bpf)+1-floor(m)*4;
                    elements(k,4) = 2;
                end
            end
            n = n + 1;
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
                

    csvwrite('nodedata.csv',nodes)
    csvwrite('elementdata.csv',elements)
    csvwrite('sectiondata.csv',sec_props)
    
    type('nodedata.csv')
    type('elementdata.csv')
    type('sectiondata.csv')


    