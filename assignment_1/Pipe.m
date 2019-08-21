classdef Pipe
    %Pipe A pipe object
    %   This is a basic pipe, for whatever.
    
    properties
        length
        diameter
        thickness
    end
    
    methods
        function obj = Pipe(l,d,t)
            %Pipe Construct an instance of this class
            %   Add a length, diameter and thickness to the pipe
            obj.length = l;
            obj.diameter = d;
            obj.thickness = t;
        end
        function C = capacity(obj)
            %capacity Calculate the capacity of the pipe
            C = obj.length * pi * ((obj.diameter-2*obj.thickness)/4)^2;
        end
        function area = A(obj)
            %A Calculate the area of the cross-section
            area = pi * (((obj.diameter/4)^2 - ((obj.diameter - 2*obj.thickness)/4)^2));
        end
    end
end

