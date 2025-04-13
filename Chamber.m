classdef Chamber

    properties
        Width=0.1 % width of chamber in meter
        Depth=0.1 % depth of chamber in meter
        MagentRegionFrom=0.7 % start ratio for magnetic region respect to depth
        MagnetRegionTo=0.9 % end ratio for magnetic region respect to depth
        MagnetPower=100 % represent magnetic field of magnets 
        StartRadius=0.15 % distance of starting radius from previous node
        AnodeWidth=0.75 % width ratio of anode at the bottom in respect to width of chamber
        AnodeVoltage=1000 % voltage of anode that must be more than cathod voltage
    end
    
    methods
        % plot chamber dimentions 
        function plot2D(obj, ax)
            cla(ax)

            % main part
            x = [-0.25 * obj.Width - obj.StartRadius, ...
                0, ...
                0, ...
                obj.Width, ...
                obj.Width, ...
                obj.Width + 0.25 * obj.Width, ...
                obj.Width + 0.25 * obj.Width, ...
                -0.25 * obj.Width - obj.StartRadius];
            y = [obj.Depth, ...
                obj.Depth, ...
                0, ...
                0, ...
                obj.Depth, ...
                obj.Depth, ...
                -0.25 * obj.Depth, ...
                -0.25 * obj.Depth];
            patch(ax, "XData", x, "YData", y, "FaceColor", "#4DBEEE", "EdgeColor", "none");
            hold(ax, "on");
            % anotation
            doubleArow(ax, [-obj.StartRadius, obj.Depth], ...
                       [-obj.StartRadius, 0], ...
                       "Depth:"+obj.Depth, "right"); % depth

            doubleArow(ax, [-obj.StartRadius, 0], ...
                       [0, 0], ...
                       "SR:"+obj.StartRadius, "top"); % start radius

            doubleArow(ax, [0, obj.Depth - 0.01], ...
                       [obj.Width, obj.Depth - 0.01], ...
                       "Width:"+obj.Width, "down"); % width
 
            
            % magnet part
            magnet_width = 0.01;
            x_left = [-magnet_width, ...
                0, ...
                0, ...
                -magnet_width];
            y_left = obj.Depth * ...
                [obj.MagnetRegionTo, ...
                obj.MagnetRegionTo, ...
                obj.MagentRegionFrom, ...
                obj.MagentRegionFrom];
            patch(ax, "XData", x_left, "YData", y_left, "FaceColor", "#0072BD", "EdgeColor", "none");

            x_right = [obj.Width, ...
                obj.Width + magnet_width, ...
                obj.Width + magnet_width, ...
                obj.Width];
            y_right = obj.Depth * ...
                [obj.MagnetRegionTo, ...
                obj.MagnetRegionTo, ...
                obj.MagentRegionFrom, ...
                obj.MagentRegionFrom];
            patch(ax, "XData", x_right, "YData", y_right, "FaceColor", "#0072BD", "EdgeColor", "none");
            % anotation
            doubleArow(ax, [obj.Width, 0], ...
                       [obj.Width, obj.Depth * obj.MagentRegionFrom], ...
                       "MRF:"+obj.MagentRegionFrom, "right"); % magnet region from 

            doubleArow(ax, [obj.Width + magnet_width, 0], ...
                       [obj.Width + magnet_width, obj.Depth * obj.MagnetRegionTo], ...
                       "MRT:"+obj.MagnetRegionTo, "right"); % magnet region to 
            

            % anode part
            anode_depth = 0.01;
            x = obj.Width/2 + (obj.AnodeWidth * obj.Width / 2) * ...
                [-1, 1, 1, -1];
            y = [0, 0, -anode_depth, -anode_depth];
            patch(ax, "XData", x, "YData", y, "FaceColor", 	"#A2142F", "EdgeColor", "none");
            % anotation
            doubleArow(ax, [x(1), anode_depth], ...
                       [x(2), anode_depth], ...
                       "AW:"+obj.AnodeWidth, "top"); % anode width 
            
            % set limits
            ax.XLim = [-0.25 * obj.Width - obj.StartRadius,  obj.Width + 0.25 * obj.Width];
            ax.YLim = [-0.25 * obj.Depth,  obj.Depth];
            ax.XTick = [];
            ax.YTick = [];
        end
    
        % create chamber space
        function gm = createGeometry(obj, R0)
            % radius of chamber space
            R_in = R0 + obj.StartRadius;
            R_out = R_in + obj.Width;

            % create empty space of chamber
            R = [R_in, ...
                R_in + obj.Width/2 - (obj.AnodeWidth * obj.Width)/2, ...
                R_in + obj.Width/2 + (obj.AnodeWidth * obj.Width)/2, ...
                R_out];
            gm = multicylinder(R, (obj.Depth * obj.MagentRegionFrom), "Void", [true, false, false, false]);
            gm = extrude(gm,[2,6,9],obj.Depth * [obj.MagnetRegionTo - obj.MagentRegionFrom, 1 - obj.MagnetRegionTo]);
            while gm.NumCells > 1
                gm = mergeCells(gm,[1 2]);
            end

            % translate space to the serface of thurster
            gm = translate(gm, [0, 0, -obj.Depth]);
        end
    end

    % set and get methods
    methods
        function obj = set.Width(obj, newValue)
            MINIMUM = 0.05;
            obj.Width = max(newValue, MINIMUM);
        end

        function obj = set.Depth(obj, newValue)
            MINIMUM = 0.05;
            obj.Depth = max(newValue, MINIMUM);
        end

        function obj = set.MagentRegionFrom(obj, newValue)
            MINIMUM = 0.02;
            MAXIMUM = 0.98;
            MIN_WIDTH = 0.02;
            newValue = min(max(newValue, MINIMUM), MAXIMUM);
            obj.MagentRegionFrom = min(newValue, obj.MagnetRegionTo - MIN_WIDTH);
        end

        function obj = set.MagnetRegionTo(obj, newValue)
            MINIMUM = 0.02;
            MAXIMUM = 0.98;
            MIN_WIDTH = 0.02;
            newValue = min(max(newValue, MINIMUM), MAXIMUM);
            obj.MagnetRegionTo = max(newValue, obj.MagentRegionFrom + MIN_WIDTH);
        end

        function obj = set.MagnetPower(obj, newValue)
            MINIMUM = 1;
            obj.MagnetPower = max(newValue, MINIMUM);
        end

        function obj = set.StartRadius(obj, newValue)
            MINIMUM = 0.05;
            obj.StartRadius = max(newValue, MINIMUM);
        end

        function obj = set.AnodeWidth(obj, newValue)
            MINIMUM = 0.02;
            MAXIMUM = 0.98;
            obj.AnodeWidth = min(max(newValue, MINIMUM), MAXIMUM);
        end

        function obj = set.AnodeVoltage(obj, newValue)
            MINIMUM = 0;
            obj.AnodeVoltage = max(newValue, MINIMUM);
        end
    end
end


