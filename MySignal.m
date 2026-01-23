classdef MySignal
    properties
        start
        stop
        Tres
        values
    end

    methods
        function obj = MySignal(start, stop, Tres, values)
            arguments
                start (1,1) double
                stop  (1,1) double
                Tres  (1,1) double {mustBePositive}
                values (1,:) double % Validate that values is a 1xN signal vector
            end
            narginchk(4,4); % Validate correct number of inputs
            obj.start = start;
            obj.stop = stop;
            obj.Tres = Tres;
            obj.values = values;
        end

        function output = MyConvolve(obj, other)
            arguments
                obj
                other (1,1) MySignal
            end
            
        end 
    end
end
