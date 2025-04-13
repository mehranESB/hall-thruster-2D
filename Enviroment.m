classdef Enviroment

    properties (Constant = true)
        VoltageContactThikness = 0.01
        WallThikness = 0.02;
        CoreThikness = 0.03;
        CoilThikness = 0.02;
        CoilDistance = 0.01;
    end

    properties
        Height=0.6 % simulation height from bottom of thruster in meter
        Radius=0.75 % simulation raduis with center of thruster in meter
        Thruster=Thruster() % Thruster objeject
    end

    properties
        Geometry % 2D geometry of space
        Geometry3D % 3D geometry of space
        Points = struct() % structure of creatical points 2D
        Points3D = struct() % structure of creatical points in 3D
    end

    methods

        % plot enviroment region in 2D
        function plot2D(obj, ax)
            % enviroment
            x = [-obj.Radius, -obj.Radius, obj.Radius, obj.Radius];
            y = [-(obj.Thruster.MaxDepth + 0.05), obj.Height, obj.Height, -(obj.Thruster.MaxDepth + 0.05)];
            patch(ax, "XData", x, "YData", y, "FaceColor", "#77AC30", "EdgeColor", "none");
            hold on

            % thruster
            obj.Thruster.plot2D(ax);

            doubleArow(ax, [0, 0], ...
                       [0, obj.Height], ...
                       "Height:"+obj.Height, "left"); % height

            % set limits 
            ax.XLim = [-obj.Radius, obj.Radius];
            ax.YLim = [-(obj.Thruster.MaxDepth + 0.05), obj.Height];
        end
        
        % create wholde geometry
        function [gd, ns,points] = createGeometry(obj)
            % check all parameters to be currect
            obj = checkParam(obj);

            % geometrix matrix for all parts
            [obj, gd_body, ns_body] = createBody(obj);
            [obj, gd_contact, ns_contact] = createVoltageContacts(obj);
            [obj, gd_magnetCore, ns_magnetCore] = createMagneticCores(obj);
            [obj, gd_magnetCoil, ns_magnetCoil] = createMagneticCoils(obj);
            [obj, gd_space, ns_space] = createSpace(obj);

            % create a whole geomtry
            maxSize = max([size(gd_body,1), size(gd_contact,1), size(gd_magnetCore,1), size(gd_magnetCoil,1), size(gd_space,1)]);
            gd_body = [gd_body; zeros(maxSize - size(gd_body, 1), size(gd_body, 2))];
            gd_contact = [gd_contact; zeros(maxSize - size(gd_contact, 1), size(gd_contact, 2))];
            gd_magnetCore = [gd_magnetCore; zeros(maxSize - size(gd_magnetCore, 1), size(gd_magnetCore, 2))];
            gd_magnetCoil = [gd_magnetCoil; zeros(maxSize - size(gd_magnetCoil, 1), size(gd_magnetCoil, 2))];
            gd_space = [gd_space; zeros(maxSize - size(gd_space, 1), size(gd_space, 2))];
            gd = [gd_body, gd_contact, gd_magnetCore, gd_magnetCoil, gd_space];

            % create ns 
            ns = [ns_body, ns_contact, ns_magnetCore, ns_magnetCoil, ns_space];
            
            % coordinates point
            points = obj.Points;
        end
    
        % create enviroment geometry
        function obj = createGeometry3D(obj)

            [gm3D, points] = obj.createBase();
            gm3D = obj.mergeAll(gm3D);

            % get faceIds 
            points = obj.getFaceIDes(gm3D, points);

            % extrude chambers and free space
            gm3D = obj.extrudeChambersAndFree(gm3D, points);
            gm3D = obj.mergeAll(gm3D);

            % create cathode void
            [gm3D, points] = obj.voidCathode(gm3D, points);

            % get all face IDes
            points = obj.getAllFaces(gm3D, points);

            % rotate and translate
            gm3D = rotate(gm3D, 180, [0,0,0], [1,0,0]);
            gm3D = translate(gm3D, [0, 0, obj.Height]);

            % generate mesh
            gm3D = generateMesh(fegeometry(gm3D), Hmax=0.02);

            % assign properties
            obj.Geometry3D = gm3D;
            obj.Points3D = points;
        end
    
        % plot geometry
        function plot3D(obj)
            mesh = obj.Geometry3D.Mesh;
            if isempty(mesh)
                return;
            end

            points = obj.Points3D;
            colorMapData = ones(size(mesh.Nodes, 2), 1);
            cmap = [1, 1, 1; ... % envioroment 1
                1, 0, 0;... % anode 2
                0, 0, 1;... % magnet 3
                1, 1, 0;... % cathode 4
                0, 0, 0]; % body 5

            % anodes
            nodeIDes = findNodes(mesh,"region","Face",points.AnodeFaceIDes);
            colorMapData(nodeIDes) = 2;

            % magnets 
            nodeIDes = findNodes(mesh,"region","Face",[points.MagnetInFaceIDes, points.MagnetOutFaceIDes]);
            colorMapData(nodeIDes) = 3;

            % cathode
            nodeIDes = findNodes(mesh,"region","Face",points.CathodeFaceIDes);
            colorMapData(nodeIDes) = 4;

            % body
            nodeIDes = findNodes(mesh,"region","Face",points.ThrusterBodyFaceIDes);
            colorMapData(nodeIDes) = 5;

            pdeplot3D(mesh, "ColorMapData", colorMapData, "FaceAlpha", 0.7);
            colormap(cmap)
            clim([1 5])
        end
    end

    % get and set methods
    methods

        function obj = set.Height(obj, newValue)
            if newValue <= obj.Thruster.MaxHeight
                error("height of enviroment must be greater than thruster.MaxHeight");

            end
            obj.Height = newValue;
        end

        function obj = set.Radius(obj, newValue)
            if newValue <= obj.Thruster.MaxRadius
                error("radius of enviroment must be greater than thruster.MaxRadius");
            end
            obj.Radius = newValue;
        end
    end
    
    % helper methods for 2d geometry
    methods (Access = private)

        % create 2d body of thruster
        function [obj, gd, ns] = createBody(obj)

            bodyFacePoints_rights = [];
        
            xSegments_right = [];
            ySegments_right = [];

            % create right part of body 
            xSegments_right(end+1) = 0;
            ySegments_right(end+1) = 0;

            % chambers
            for i = 1 : length(obj.Thruster.Chambers)

                chamber = obj.Thruster.Chambers{i};

                % start chamber segment
                xSegments_right(end+1) = xSegments_right(end) + chamber.StartRadius;
                ySegments_right(end+1) = 0;
    
                % go to depth of chamber
                xSegments_right(end+1) = xSegments_right(end);
                ySegments_right(end+1) = -chamber.Depth;

                % width of chamber
                xSegments_right(end+1) = xSegments_right(end) + chamber.Width;
                ySegments_right(end+1) = -chamber.Depth;

                % add points 
                bodyFacePoints_rights = [bodyFacePoints_rights;[xSegments_right(end) - chamber.Width/2, ySegments_right(end)-1e-3]];

                % go back into surface
                xSegments_right(end+1) = xSegments_right(end);
                ySegments_right(end+1) = 0;
            end

            buttomThikness = 2*obj.CoilDistance + 2*obj.CoilThikness + obj.CoreThikness;

            % edge part of thruster
            xSegments_right(end+1) = xSegments_right(end) + obj.Thruster.FinalWidth;
            ySegments_right(end+1) = 0;

            xSegments_right(end+1) = xSegments_right(end);
            ySegments_right(end+1) = -obj.Thruster.MaxDepth - buttomThikness;

            % close segment
            xSegments_right(end+1) = 0;
            ySegments_right(end+1) = -obj.Thruster.MaxDepth - buttomThikness;

            % create left part 
            xSegments_left = -1 .* xSegments_right;
            ySegments_left = ySegments_right;
            bodyFacePoints_lefts = [-1 .* bodyFacePoints_rights(:,1), bodyFacePoints_rights(:,2)];

            % create cathod body
            xSegment_cathode = obj.Thruster.Cathode.Radius * [-1,1,1,-1];
            ySegment_cathode = [...
                obj.Thruster.Cathode.Height, ...
                obj.Thruster.Cathode.Height, ...
                -obj.Thruster.MaxDepth - buttomThikness, ...
                -obj.Thruster.MaxDepth - buttomThikness];

            % set the cathod in the right position
            if obj.Thruster.Cathode.IsInMiddle
                xSegments_right = xSegments_right + obj.Thruster.Cathode.Radius;
                bodyFacePoints_rights(:,1) = bodyFacePoints_rights(:,1) + obj.Thruster.Cathode.Radius;
                xSegments_left = xSegments_left - obj.Thruster.Cathode.Radius;
                bodyFacePoints_lefts(:,1) = bodyFacePoints_lefts(:,1) - obj.Thruster.Cathode.Radius;
            else
                xSegment_cathode = xSegment_cathode - (obj.Thruster.PanelRadius + obj.Thruster.Cathode.Radius);
            end

            % add points
            bodyFacePoints_rights = [bodyFacePoints_rights; [mean(xSegments_right([1,end]))+1e-3, mean(ySegments_right([1,end]))]];
            bodyFacePoints_lefts = [bodyFacePoints_lefts; [mean(xSegments_left([1,end]))-1e-3, mean(ySegments_left([1,end]))]];
            obj.Points.Face.body = [...
                    bodyFacePoints_rights;...
                    bodyFacePoints_lefts;...
                    [mean(xSegment_cathode), mean(ySegment_cathode)]];

            % create geometry description matrix
            gd_right = [2, length(xSegments_right), xSegments_right, ySegments_right];
            gd_left = [2, length(xSegments_left), xSegments_left, ySegments_left];
            gd_cathode = [3, 4, xSegment_cathode, ySegment_cathode, zeros(1, 2*length(xSegments_right)-8)];
            gd = [gd_right; gd_left; gd_cathode]';
            ns = {'rightBody', 'leftBody', 'cathodeBody'};
        end
    
        % create electric vlotage place 
        function [obj, gd, ns] = createVoltageContacts(obj)
        
            gd = [];
            ns = {};
            
            % catode tip contact
            xSegment_cathode = obj.Thruster.Cathode.Radius * 0.75 * [-1,1,1,-1];
            ySegment_cathode =  obj.Thruster.Cathode.Height + [0, 0, obj.VoltageContactThikness, obj.VoltageContactThikness];
            if ~obj.Thruster.Cathode.IsInMiddle
                xSegment_cathode = xSegment_cathode - (obj.Thruster.PanelRadius + obj.Thruster.Cathode.Radius);
            end
            newRect = [3, 4, xSegment_cathode, ySegment_cathode];
            gd(:, end+1) = newRect';
            ns{end+1} = 'cathodeContact';

            % add points
            obj.Points.Edge.cathodeContact = [mean(xSegment_cathode(3:4)), mean(ySegment_cathode(3:4))];
            obj.Points.Face.cathodeContact = [mean(xSegment_cathode), mean(ySegment_cathode)];
            
            % create anode contacts in chambers 
            R = 0;
            for i = 1 : length(obj.Thruster.Chambers)
                
                chamber = obj.Thruster.Chambers{i};
                width = chamber.Width;

                % middle point of chamber
                middleX = R + chamber.StartRadius + (width / 2);
                if obj.Thruster.Cathode.IsInMiddle
                    middleX = middleX + obj.Thruster.Cathode.Radius;
                end

                % x segment info
                xSegment_chamber = middleX + ((width * chamber.AnodeWidth) / 2) * [-1, 1, 1, -1];
                ySegment_chamber = -chamber.Depth + [0, 0, obj.VoltageContactThikness, obj.VoltageContactThikness];
                
                % add to geometry description matrix
                newRect = [3, 4, xSegment_chamber, ySegment_chamber];
                gd(:, end+1) = newRect';
                ns{end+1} = ['chamberContactR', num2str(i)];

                % add points
                obj.Points.Edge.("chamberContactR" + i) = [mean(xSegment_chamber(3:4)), mean(ySegment_chamber(3:4))];
                obj.Points.Face.("chamberContactR" + i) = [mean(xSegment_chamber), mean(ySegment_chamber)];

                % left side 
                newRect = [3, 4, -xSegment_chamber, ySegment_chamber];
                gd(:, end+1) = newRect';
                ns{end+1} = ['chamberContactL', num2str(i)];

                % add points
                obj.Points.Edge.("chamberContactL" + i) = [mean(-xSegment_chamber(3:4)), mean(ySegment_chamber(3:4))];
                obj.Points.Face.("chamberContactL" + i) = [mean(-xSegment_chamber), mean(ySegment_chamber)];
            
                % update R
                R = R + chamber.StartRadius + width;
            end
        end
    
        % create magnetic cores
        function [obj, gd, ns] = createMagneticCores(obj)
        
            gd = [];
            ns = {};

            R = 0;
            for i = 1 : length(obj.Thruster.Chambers)
                
                chamber = obj.Thruster.Chambers{i};
                
                xSegment_core_right = [];
                ySegment_code_right = [];
                
                % create core as polyline 
                xSegment_core_right(end+1) = chamber.StartRadius + R;
                ySegment_code_right(end+1) = chamber.Depth * (chamber.MagentRegionFrom-1);

                xSegment_core_right(end+1) = xSegment_core_right(end) - obj.WallThikness;
                ySegment_code_right(end+1) = ySegment_code_right(end);

                xSegment_core_right(end+1) = xSegment_core_right(end);
                ySegment_code_right(end+1) = -chamber.Depth - (obj.CoilDistance + obj.CoilThikness);

                xSegment_core_right(end+1) = xSegment_core_right(1) + chamber.Width + obj.WallThikness;
                ySegment_code_right(end+1) = ySegment_code_right(end);

                xSegment_core_right(end+1) = xSegment_core_right(end);
                ySegment_code_right(end+1) = chamber.Depth * (chamber.MagentRegionFrom-1);

                xSegment_core_right(end+1) = xSegment_core_right(end) - obj.WallThikness;
                ySegment_code_right(end+1) = ySegment_code_right(end);

                xSegment_core_right(end+1) = xSegment_core_right(end);
                ySegment_code_right(end+1) = chamber.Depth * (chamber.MagnetRegionTo-1);

                xSegment_core_right(end+1) = xSegment_core_right(end) + obj.WallThikness + obj.CoreThikness;
                ySegment_code_right(end+1) = ySegment_code_right(end);

                xSegment_core_right(end+1) = xSegment_core_right(end);
                ySegment_code_right(end+1) = -chamber.Depth - (obj.CoilDistance + obj.CoilThikness + obj.CoreThikness);

                xSegment_core_right(end+1) = xSegment_core_right(2) - obj.CoreThikness;
                ySegment_code_right(end+1) = ySegment_code_right(end);

                xSegment_core_right(end+1) = xSegment_core_right(end);
                ySegment_code_right(end+1) = chamber.Depth * (chamber.MagnetRegionTo-1);

                xSegment_core_right(end+1) = xSegment_core_right(1);
                ySegment_code_right(end+1) = ySegment_code_right(end);
            
                % modify core regard to core
                if obj.Thruster.Cathode.IsInMiddle
                    xSegment_core_right = xSegment_core_right + obj.Thruster.Cathode.Radius;
                end

                % add to geometry description matrix
                newRect = [2, length(xSegment_core_right), xSegment_core_right, ySegment_code_right];
                gd(:, end+1) = newRect';
                ns{end+1} = ['magnetCoreR', num2str(i)];

                % add points
                obj.Points.Face.("magnetCoreR" + i) = [mean(xSegment_core_right(1:2)), mean(ySegment_code_right([1,end]))];

                % left side 
                xSegment_core_left = -xSegment_core_right;
                ySegment_code_left = ySegment_code_right;
                newRect = [2, length(xSegment_core_left), xSegment_core_left, ySegment_code_left];
                gd(:, end+1) = newRect';
                ns{end+1} = ['magnetCoreL', num2str(i)];

                % add points
                obj.Points.Face.("magnetCoreL" + i) = [mean(xSegment_core_left(1:2)), mean(ySegment_code_left([1,end]))];

                % update R
                R = R + chamber.StartRadius + chamber.Width;
            end
        end
    
        % create magnetic coils
        function [obj, gd, ns] = createMagneticCoils(obj)
            gd = [];
            ns = {};

            R = 0;
            for i = 1 : length(obj.Thruster.Chambers)

                chamber = obj.Thruster.Chambers{i};
            
                xSegment_coil_top_right = [];
                ySegment_coil_top_right = [];

                xSegment_coil_down_right = [];
                ySegment_coil_down_right = [];

                % create coil as rectangle for top
                xSegment_coil_top_right(end+1) = chamber.StartRadius + R;
                ySegment_coil_top_right(end+1) = -chamber.Depth - obj.CoilDistance;
                
                xSegment_coil_top_right(end+1) = xSegment_coil_top_right(end) + chamber.Width;
                ySegment_coil_top_right(end+1) = ySegment_coil_top_right(end);

                xSegment_coil_top_right(end+1) = xSegment_coil_top_right(end);
                ySegment_coil_top_right(end+1) = ySegment_coil_top_right(end) - obj.CoilThikness;

                xSegment_coil_top_right(end+1) = xSegment_coil_top_right(1);
                ySegment_coil_top_right(end+1) = ySegment_coil_top_right(end);

                % create coil as rectangle for down
                xSegment_coil_down_right(end+1) = chamber.StartRadius + R;
                ySegment_coil_down_right(end+1) = -chamber.Depth - (obj.CoilDistance + obj.CoilThikness + obj.CoreThikness);
                
                xSegment_coil_down_right(end+1) = xSegment_coil_down_right(end) + chamber.Width;
                ySegment_coil_down_right(end+1) = ySegment_coil_down_right(end);

                xSegment_coil_down_right(end+1) = xSegment_coil_down_right(end);
                ySegment_coil_down_right(end+1) = ySegment_coil_down_right(end) - obj.CoilThikness;

                xSegment_coil_down_right(end+1) = xSegment_coil_down_right(1);
                ySegment_coil_down_right(end+1) = ySegment_coil_down_right(end);

                % modify core regard to core
                if obj.Thruster.Cathode.IsInMiddle
                    xSegment_coil_top_right = xSegment_coil_top_right + obj.Thruster.Cathode.Radius;
                    xSegment_coil_down_right = xSegment_coil_down_right + obj.Thruster.Cathode.Radius;
                end

                % add to geometry description matrix
                newRect = [3, 4, xSegment_coil_top_right, ySegment_coil_top_right];
                gd(:, end+1) = newRect';
                ns{end+1} = ['magnetCoilTopR', num2str(i)];

                newRect = [3, 4, xSegment_coil_down_right, ySegment_coil_down_right];
                gd(:, end+1) = newRect';
                ns{end+1} = ['magnetCoilDownR', num2str(i)];

                % add points
                obj.Points.Face.("magnetCoilTopR" + i) = [mean(xSegment_coil_top_right), mean(ySegment_coil_top_right)];
                obj.Points.Rect.("magnetCoilTopR" + i) = [xSegment_coil_top_right', ySegment_coil_top_right'];

                obj.Points.Face.("magnetCoilDownR" + i) = [mean(xSegment_coil_down_right), mean(ySegment_coil_down_right)];
                obj.Points.Rect.("magnetCoilDownR" + i) = [xSegment_coil_down_right', ySegment_coil_down_right'];

                % left side 
                xSegment_coil_top_left = -xSegment_coil_top_right;
                ySegment_coil_top_left = ySegment_coil_top_right;
                newRect = [3, 4, xSegment_coil_top_left, ySegment_coil_top_left];
                gd(:, end+1) = newRect';
                ns{end+1} = ['magnetCoilTopL', num2str(i)];

                xSegment_coil_down_left = -xSegment_coil_down_right;
                ySegment_coil_down_left = ySegment_coil_down_right;
                newRect = [3, 4, xSegment_coil_down_left, ySegment_coil_down_left];
                gd(:, end+1) = newRect';
                ns{end+1} = ['magnetCoilDownL', num2str(i)];

                % add points
                obj.Points.Face.("magnetCoilTopL" + i) = [mean(xSegment_coil_top_left), mean(ySegment_coil_top_left)];
                obj.Points.Rect.("magnetCoilTopL" + i) = [xSegment_coil_top_left', ySegment_coil_top_left'];

                obj.Points.Face.("magnetCoilDownL" + i) = [mean(xSegment_coil_down_left), mean(ySegment_coil_down_left)];
                obj.Points.Rect.("magnetCoilDownL" + i) = [xSegment_coil_down_left', ySegment_coil_down_left'];

                % update R
                R = R + chamber.StartRadius + chamber.Width;
            end

        end
    
        % create space 
        function [obj, gd, ns] = createSpace(obj)
        
            buttomThikness = 2*obj.CoilDistance + 2*obj.CoilThikness + obj.CoreThikness;

            xSegment_space = obj.Radius * [-1, 1, 1, -1];
            ySegment_space = [-obj.Thruster.MaxDepth - buttomThikness-0.05,...
                -obj.Thruster.MaxDepth - buttomThikness - 0.05,...
                obj.Height, ...
                obj.Height];

            % add to geometry description matrix
            newRect = [3, 4, xSegment_space, ySegment_space];
            gd = newRect';
            ns = {'space'};

            % add points
            obj.Points.Face.("space") = [0, obj.Height];
            obj.Points.Rect.("space") = [xSegment_space', ySegment_space'];
        end
    
        % check parameters 
        function obj = checkParam(obj)
        
            % check start radious for closeness of chambers
            minimumValue = obj.CoreThikness + obj.WallThikness + 0.01;
            for i = 1 : length(obj.Thruster.Chambers)
                chamber = obj.Thruster.Chambers{i};
                if chamber.StartRadius < minimumValue
                    chamber.StartRadius = minimumValue;
                    obj.Thruster.Chambers{i} = chamber;
                end
            end

            % check for enviromental radius 
            obj.Radius = max(obj.Radius, obj.Thruster.MaxRadius + 0.05);
        end
    end

    % helper methods for 3d geometry
    methods (Access = private)

        % create 2D geometry for surface of thruster
        function [gm, points] = createBase(obj)
            % initialize
            if ~obj.Thruster.Cathode.IsInMiddle
                R0 = 0;
            else
                R0 = obj.Thruster.Cathode.Radius;
            end
            R = [];

            % chambers
            chamberFacePoints = zeros(0,3); % [before-anode, anode , after-anode]  
            for i = 1 : length(obj.Thruster.Chambers)
                chamber = obj.Thruster.Chambers{i};
                % start chamber
                R0 = R0 + chamber.StartRadius;
                R = [R, R0];
                % start anode
                R0 = R0 + chamber.Width * (1 - chamber.AnodeWidth) / 2;
                R = [R, R0];
                chamberFacePoints(i, 1) = (R(end) + R(end-1)) / 2;
                % end anode
                R0 = R0 + chamber.Width * chamber.AnodeWidth;
                R = [R, R0];
                chamberFacePoints(i, 2) = (R(end) + R(end-1)) / 2;
                % end chamber
                R0 = R0 + chamber.Width * (1 - chamber.AnodeWidth) / 2;
                R = [R, R0];
                chamberFacePoints(i, 3) = (R(end) + R(end-1)) / 2;
            end

            % free region
            R0  = R0 + obj.Thruster.FinalWidth;
            R = [R, R0];
            R0  = obj.Radius;
            R = [R, R0];
            freeRegionPoint = (R(end) + R(end-1)) / 2;

            % create geometry
            gm = multicylinder(R,obj.Height);
            points.Chambers = chamberFacePoints;
            points.Free = freeRegionPoint;
        end
        
        % merge all shared faces
        function gm = mergeAll(~, gm)
            while (gm.NumCells > 1)
                % pick two cell randomly
                id1 = randi(gm.NumCells);
                id2 = randi(gm.NumCells);
                if id1 == id2
                    continue
                end

                % check mutual faces
                faceIDs1 = cellFaces(gm,id1);
                faceIDs2 = cellFaces(gm,id2);
                Lia = ismember(faceIDs1,faceIDs2);
                
                % merge faces
                if any(Lia)
                    gm = mergeCells(gm, [id1,id2]);
                end
            end
        end
    
        % get face ides of points
        function points = getFaceIDes(obj, gm, points)
            % get chamber face points
            P = points.Chambers;
            coords = [zeros(numel(P), 1), P(:), obj.Height * ones(numel(P), 1)];
            FaceID = nearestFace(gm, coords);
            points.ChambersFaceIDes = reshape(FaceID, size(P));

            % get free region face point
            P = points.Free;
            coords = [zeros(numel(P), 1), P(:), obj.Height * ones(numel(P), 1)];
            FaceID = nearestFace(gm, coords);
            points.FreeFaceIDes = reshape(FaceID, size(P));
        end
    
        % extrude chambers and free space
        function gm = extrudeChambersAndFree(obj, gm, points)
            % chambers 
            FaceIDes = points.ChambersFaceIDes;
            for i = 1 : length(obj.Thruster.Chambers)
                chamber = obj.Thruster.Chambers{i};

                % extrude heights
                H1 = chamber.Depth * chamber.MagentRegionFrom; % before magnet
                H2 = chamber.Depth * (chamber.MagnetRegionTo - chamber.MagentRegionFrom); % magnet
                H3 = chamber.Depth * (1 - chamber.MagnetRegionTo); % after magnet
                H = [H1, H2, H3];

                % points
                
                % extrude chambers 
                gm = extrude(gm, FaceIDes(i, :), H);
            end

            % free 
            FaceIDes = points.FreeFaceIDes;
            H = obj.Thruster.MaxDepth + 0.05;
            gm = extrude(gm, FaceIDes, H);
        end
    
        % create void of cathode
        function [gm, points] = voidCathode(obj, gm, points)
            % create cathode geometry
            cathode = obj.Thruster.Cathode;
            gmCathode = multicylinder(cathode.Radius, cathode.Height);

            % points
            cathodePoints = [0, 0, 0; ...
                0, cathode.Radius, cathode.Height/2; ...
                0, 0, cathode.Height];
            cathodeTipPoints = [0, 0, 0]; 

            % translate cathode to correct geometry
            delta = 0.01;
            if cathode.IsInMiddle
                T = [0, ...
                    0, ...
                    obj.Height - cathode.Height - delta];
            else
                T = [0, ...
                    obj.Thruster.PanelRadius + cathode.Radius, ...
                    obj.Height - cathode.Height - delta];
            end
            gmCathode = translate(gmCathode, T);

            % make cathode void
            gm = addVoid(gm, gmCathode);

            % add translated points\
            points.CathodePoints = cathodePoints + T;
            points.CathodeTipPoints = cathodeTipPoints + T;
        end
    
        % get all face points
        function points = getAllFaces(obj, gm, points)
            thrusterBodyPoints = [];
            
            % chamber 
            for i = 1 : length(obj.Thruster.Chambers)
                chamber = obj.Thruster.Chambers{i};

                % anode
                anodePoints(i, :) = [0, ...
                    points.Chambers(i,2), ...
                    chamber.Depth + obj.Height];

                % magnet in
                magnetInPoints(i, :) = [0, ...
                    points.Chambers(i,2) - chamber.Width/2, ...
                    chamber.Depth * (chamber.MagentRegionFrom + chamber.MagnetRegionTo)/2 + obj.Height];

                % magnet out
                magnetOutPoints(i, :) = [0, ...
                    points.Chambers(i,2) + chamber.Width/2, ...
                    chamber.Depth * (chamber.MagentRegionFrom + chamber.MagnetRegionTo)/2 + obj.Height];

                % thruster body
                thrusterBodyPoints = [thrusterBodyPoints; [0, ...
                    points.Chambers(i,1) - chamber.Width * (1-chamber.AnodeWidth)/2 - 0.01, ...
                    obj.Height]]; % top center face
                thrusterBodyPoints = [thrusterBodyPoints; [0, ...
                    points.Chambers(i,2) - chamber.Width/2, ...
                    obj.Height + chamber.Depth * (1 - chamber.MagnetRegionTo)/2]]; % after magnet in
                thrusterBodyPoints = [thrusterBodyPoints; [0, ...
                    points.Chambers(i,2) - chamber.Width/2, ...
                    obj.Height + chamber.Depth * (1 - chamber.MagentRegionFrom/2)]]; % before magnet in
                thrusterBodyPoints = [thrusterBodyPoints; [0, ...
                    points.Chambers(i,1), ...
                    obj.Height + chamber.Depth]]; % before anode
                thrusterBodyPoints = [thrusterBodyPoints; [0, ...
                    points.Chambers(i,3), ...
                    obj.Height + chamber.Depth]]; % after anode
                thrusterBodyPoints = [thrusterBodyPoints; [0, ...
                    points.Chambers(i,2) + chamber.Width/2, ...
                    obj.Height + chamber.Depth * (1 - chamber.MagentRegionFrom/2)]]; % before magnet out
                thrusterBodyPoints = [thrusterBodyPoints; [0, ...
                    points.Chambers(i,2) + chamber.Width/2, ...
                    obj.Height + chamber.Depth * (1 - chamber.MagnetRegionTo)/2]]; % after magnet out
            end
            
            % thruster final body
            thrusterBodyPoints = [thrusterBodyPoints; [obj.Thruster.PanelRadius, ...
                0, ...
                obj.Height + obj.Thruster.MaxDepth/2]]; % side outer face
            thrusterBodyPoints = [thrusterBodyPoints; [0, ...
                obj.Thruster.PanelRadius - obj.Thruster.FinalWidth/2, ...
                obj.Height]]; % top outer face


            % set points
            points.AnodePoints = anodePoints;
            points.MagnetInPoints = magnetInPoints;
            points.MagnetOutPoints = magnetOutPoints;
            points.ThrusterBodyPoints = thrusterBodyPoints;

            % find faces
            points.AnodeFaceIDes = nearestFace(gm, anodePoints);
            points.MagnetInFaceIDes = nearestFace(gm, magnetInPoints);
            points.MagnetOutFaceIDes = nearestFace(gm, magnetOutPoints);
            points.ThrusterBodyFaceIDes = nearestFace(gm, thrusterBodyPoints);
            points.CathodeFaceIDes = nearestFace(gm, points.CathodePoints);
            points.CathodeTipFaceIDes = nearestFace(gm, points.CathodeTipPoints);
        end
    end
end

