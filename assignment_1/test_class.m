classdef test_class
    %CLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = buildProperty(obj,inputArg1,inputArg2)
            %CLASS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function obj = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            temp = obj.Property1;
            obj.Property1 = inputArg + temp;
        end
    end
end

