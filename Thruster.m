classdef Thruster

    properties
        FinalWidth=0.1 % the width of last part after finishing last chamber
        Cathode=Cathode() % Cathode object
        Chambers={Chamber(), Chamber()} % cell array of chambers
    end

    properties (Dependent)
        MaxHeight % maximum height of thruster from bottom
        MaxRadius % maximum raduis of thruster from center
        MaxDepth % maximum depth of chambers 
        PanelRadius % raduis of panel containig chambers 
    end
    
    methods
        % plot 2D cut of thruster
        function plot2D(obj, ax)
            % main part
            surface_y = obj.MaxHeight - obj.Cathode.Height;
            maxRadius = obj.MaxRadius;
            if ~obj.Cathode.IsInMiddle % cathode is not in middle
                maxRadius = maxRadius - 2 * obj.Cathode.Radius;
            end
            x = [-maxRadius, -maxRadius + obj.FinalWidth];
            y = [0, 0];
            % left side 
            for i = length(obj.Chambers) : -1 : 1
                chamber = obj.Chambers{i};

                x(end+1) = x(end);
                y(end+1) = y(end) - chamber.Depth;

                x(end+1) = x(end) + chamber.Width;
                y(end+1) = y(end);

                x(end+1) = x(end);
                y(end+1) = y(end) + chamber.Depth;

                x(end+1) = x(end) + chamber.StartRadius;
                y(end+1) = y(end);
            end
            x(end+1) = -x(end);
            y(end+1) = y(end);
            % right side 
            for i = 1 : length(obj.Chambers)
                chamber = obj.Chambers{i};

                x(end+1) = x(end) + chamber.StartRadius;
                y(end+1) = y(end);

                x(end+1) = x(end);
                y(end+1) = y(end) - chamber.Depth;

                x(end+1) = x(end) + chamber.Width;
                y(end+1) = y(end);

                x(end+1) = x(end);
                y(end+1) = y(end) + chamber.Depth;
            end
            x(end+1) = x(end) + obj.FinalWidth;
            y(end+1) = y(end);
            % close part
            delta1 = 0.05;
            x = [x, x(end), x(1)];
            y = [y, -surface_y - delta1, -surface_y - delta1];
            patch(ax, "XData", x, "YData", y, "FaceColor", "#4DBEEE", "EdgeColor", "none");
            hold(ax, "on");
            % anotation
            doubleArow(ax, [x(end-3), y(end-3)], ...
                       [x(end-2), y(end-3)], ...
                       "FW:"+obj.FinalWidth, "top"); % final width

            doubleArow(ax, [0, 0], ...
                       [0, -obj.MaxDepth], ...
                       "MD:"+obj.MaxDepth, "left"); % maximum depth
            
            doubleArow(ax, [0, 0.01], ...
                       [obj.MaxRadius, 0.01], ...
                       "MR:"+obj.MaxRadius, "top"); % maximum radius

            % magnets 
            x = [];
            y = [];
            delta2 = 0.01;
            if ~obj.Cathode.IsInMiddle % cathode is not in middle
                R = 0;
            else
                R = obj.Cathode.Radius;
            end
            for i = 1 : length(obj.Chambers)
                chamber = obj.Chambers{i};

                % right left
                x = [x, R + chamber.StartRadius + delta2 * [0;-1;-1;0]];
                y = [y, chamber.Depth * [chamber.MagnetRegionTo;chamber.MagnetRegionTo;chamber.MagentRegionFrom;chamber.MagentRegionFrom] - chamber.Depth];
                % right right
                x = [x, R + chamber.StartRadius + chamber.Width + delta2 * [0;1;1;0]];
                y = [y, chamber.Depth * [chamber.MagnetRegionTo;chamber.MagnetRegionTo;chamber.MagentRegionFrom;chamber.MagentRegionFrom] - chamber.Depth];
                % left left
                x = [x, -(R + chamber.StartRadius + delta2 * [0;-1;-1;0])];
                y = [y, chamber.Depth * [chamber.MagnetRegionTo;chamber.MagnetRegionTo;chamber.MagentRegionFrom;chamber.MagentRegionFrom] - chamber.Depth];
                % left right
                x = [x, -(R + chamber.StartRadius + chamber.Width + delta2 * [0;1;1;0])];
                y = [y, chamber.Depth * [chamber.MagnetRegionTo;chamber.MagnetRegionTo;chamber.MagentRegionFrom;chamber.MagentRegionFrom] - chamber.Depth];
                    
                R = R + chamber.StartRadius + chamber.Width;
            end
            patch(ax, "XData", x, "YData", y, "FaceColor", "#0072BD", "EdgeColor", "none");
        
            % anode
            x = [];
            y = [];
            delta2 = 0.01;
            if ~obj.Cathode.IsInMiddle % cathode is not in middle
                R = 0;
            else
                R = obj.Cathode.Radius;
            end
            for i = 1 : length(obj.Chambers)
                chamber = obj.Chambers{i};
                
                % right
                x = [x, R + chamber.StartRadius + chamber.Width/2 + (chamber.AnodeWidth * chamber.Width / 2) * [-1; 1; 1; -1]];
                y = [y, [0; 0; -delta2; -delta2] - chamber.Depth];
                % left
                x = [x, -(R + chamber.StartRadius +  chamber.Width/2 + (chamber.AnodeWidth * chamber.Width / 2) * [-1; 1; 1; -1])];
                y = [y, [0; 0; -delta2; -delta2] - chamber.Depth];

                R = R + chamber.StartRadius + chamber.Width;
            end
            patch(ax, "XData", x, "YData", y, "FaceColor", 	"#A2142F", "EdgeColor", "none");

            % cathode
            x = (obj.Cathode.Radius / 2) * [-1,-1,1,1];
            y = [0, obj.Cathode.Height, obj.Cathode.Height, 0];
            if ~obj.Cathode.IsInMiddle % cathode is not in middle
                x = x - (obj.MaxRadius - obj.Cathode.Radius) + obj.Cathode.Radius/2;
            end 
            patch(ax, "XData", x, "YData", y, "FaceColor", 	"#EDB120", "EdgeColor", "none");
            % anotation
            doubleArow(ax, [-obj.MaxRadius, 0], ...
                       [-obj.MaxRadius, obj.Cathode.Height], ...
                       "CH:"+obj.Cathode.Height, "right"); % final width
        end

        % add new chamber to thruster
        function obj = addChamber(obj, newChamber)
            obj.Chambers{end+1} = newChamber;
        end

        % create geometry of trusters 
        function [gm_chambers, gm_cathode] = createGeometry(obj)
            % initialize R0
            if ~obj.Cathode.IsInMiddle % cathode is not in middle
                R0 = 0;
                T_cathode = [0, obj.MaxRadius + obj.Cathode.Radius, 0.01];
            else
                R0 = obj.Cathode.Radius;
                T_cathode = [0, 0, 0.01];
            end

            % create chamber spaces
            gm_chambers = cell(1, length(obj.Chambers));
            for i = 1 : length(obj.Chambers)
                chamber = obj.Chambers{i};
                gm_chambers{i} = chamber.createGeometry(R0);
                R0 = R0 + chamber.StartRadius + chamber.Width;
            end

            % create cathode geometry
            gm_cathode = multicylinder(obj.Cathode.Radius, obj.Cathode.Height);
            gm_cathode = translate(gm_cathode, T_cathode);
        end
    end

    % get and set methods
    methods
        function value = get.MaxRadius(obj)
            value = obj.FinalWidth;

            % sum of all occupied radius by chambers 
            for i = 1 : length(obj.Chambers) 
                chamber = obj.Chambers{i};
                value = value + chamber.StartRadius + chamber.Width;
            end

            % occupied radius with cathode
            if obj.Cathode.IsInMiddle
                value = value + obj.Cathode.Radius;
            else
                % multiply by 2 because the chathode is out of thruster and
                % it is tanganted to edge of thruseter
                value = value + 2 * obj.Cathode.Radius;
            end
        end
        
        function value = get.MaxHeight(obj)
            value = -inf;

            % maximum depth of chambers
            for i = 1 : length(obj.Chambers)
                chamber = obj.Chambers{i};
                value = max(value, chamber.Depth);
            end
            
            % consider cathode height
            value = value + obj.Cathode.Height;
        end
    
        function value = get.MaxDepth(obj)
            depthes = zeros(1, length(obj.Chambers));
            for i = 1 : length(obj.Chambers)
                chamber = obj.Chambers{i};
                depthes(i) = chamber.Depth;
            end
            value = max(depthes);
        end
        
        function value = get.PanelRadius(obj)
            value = obj.FinalWidth;

            % sum of all occupied radius by chambers 
            for i = 1 : length(obj.Chambers) 
                chamber = obj.Chambers{i};
                value = value + chamber.StartRadius + chamber.Width;
            end

            % occupied radius with cathode
            if obj.Cathode.IsInMiddle
                value = value + obj.Cathode.Radius;
            end
        end
    end
end

