classdef Cathode

    properties
        Voltage=0 % voltage of head of cathode must be less than anode voltage
        IsInMiddle=false % bool value for representing cathode is in middle of thruseter
        Radius=0.02 % radius of cathode
        Height=0.1 % height of cathode respect to surface of hall thruster
    end

    % set and get methods
    methods
        function obj = set.Radius(obj, newValue)
            if newValue <= 0
                error("raduis of cathode must be posetive number");
            end
            obj.Radius = newValue;
        end

        function obj = set.Height(obj, newValue)
            if newValue <= 0
                error("height of cathode respect to thruster surface must be posetive number");
            end
            obj.Height = newValue;
        end

        function obj = set.Voltage(obj, newValue)
            if newValue > 0
                error("voltage of cathode must be negative number or zero");
            end
            obj.Voltage = newValue;
        end
    end
end

