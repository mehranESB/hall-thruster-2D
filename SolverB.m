classdef SolverB
    properties % mesh properties
        DetailedLength = 0.01
        Hmin = 0.002
        Hmax = 0.03
        Hgrad = 1.1
    end

    properties
        CorePermeability = 10000;
    end

    properties
        Env
        Model
        Points
        SysGeometry
        Results
    end

    methods
      
        function obj = SolverB(env)
            obj.Env = env;
            
            % create electro statis model
            model = createpde("electromagnetic","magnetostatic");
            model.VacuumPermeability = 1.2566370614E-6;
            obj.Model = model;

        end

        function obj = solve(obj)
            % set model geometry
            obj = setGeometry(obj);

            % set boundary condition
            obj = setBC(obj);
            electromagneticProperties(obj.Model,"RelativePermeability",1.0); 
            obj = setCorePermeability(obj);

            % generate mesh
            obj = genMesh(obj);

            % solve problem 
            obj.Results = solve(obj.Model);
        end

        function plotResults(obj)
            figure;
            
            pdeplot(obj.Model, ...
                "XYData",obj.Results.MagneticPotential, ...
                "Contour","on", ...
                "FlowData", [obj.Results.MagneticField.Hx, obj.Results.MagneticField.Hy], ...
                "FaceAlpha", 0.3, ...
                "Levels",20);
            hold on
            pdegplot(obj.SysGeometry);
            xlabel("r")
            ylabel("z")
            title("magnetic potential and magnetic field")

        end
    end
    
    % helper methods
    methods
        function obj = setGeometry(obj)
        
            % create whole geometry
            [gd, ns, points] = createGeometry(obj.Env);
            obj.Points = points;
            obj.SysGeometry = decsg(gd);

            % create remaining space geometry
            mask = false(1, length(ns));
            sf = '(';
            for i = 1 : length(ns)-1
                if contains(ns{i}, 'magnetCore') || contains(ns{i}, 'magnetCoil')
                    sf = [sf, ns{i}, '+'];
                    mask(i) = true;
                end
            end
            sf(end) = ')';
            sf = [sf, '+space'];
            mask(end) = true;
            ns = char(ns(mask));
            [dl,~] = decsg(gd(:,mask),sf,ns');


            % set model geometry
            geometryFromEdges(obj.Model, dl);
            
        end
        
        function obj = setBC(obj)
        
            points = obj.Points;
            model = obj.Model;
            pg = model.Geometry;

            % last point of space 
            EdgeIDs_space = nearestEdge(pg, points.Rect.space);
            electromagneticBC(model,"MagneticPotential",0,...
                 "Edge",EdgeIDs_space);

            % coils
            for i = 1 : length(obj.Env.Thruster.Chambers)
                totalCarrent = obj.Env.Thruster.Chambers{i}.MagnetPower;
                
                electromagneticSource(model,"CurrentDensity", -computeCurrentDensity(totalCarrent, points.Rect.("magnetCoilTopR" + i)), ...
                    "Face",nearestFace(pg, points.Face.("magnetCoilTopR" + i)));
                electromagneticSource(model,"CurrentDensity", computeCurrentDensity(totalCarrent, points.Rect.("magnetCoilDownR" + i)), ...
                    "Face",nearestFace(pg, points.Face.("magnetCoilDownR" + i)));
                electromagneticSource(model,"CurrentDensity", computeCurrentDensity(totalCarrent, points.Rect.("magnetCoilTopL" + i)), ...
                    "Face",nearestFace(pg, points.Face.("magnetCoilTopL" + i)));
                electromagneticSource(model,"CurrentDensity", -computeCurrentDensity(totalCarrent, points.Rect.("magnetCoilDownL" + i)), ...
                    "Face",nearestFace(pg, points.Face.("magnetCoilDownL" + i)));

            end

            function CD = computeCurrentDensity(totalCarrent, rectPoints)
            
                x = rectPoints(:,1); % Extract X-coordinates
                y = rectPoints(:,2); % Extract Y-coordinates
                
                % Ensure the polygon is closed by appending the first point at the end
                x = [x; x(1)];
                y = [y; y(1)];
                
                % Compute area using Shoelace formula
                A = 0.5 * abs(sum(x(1:end-1) .* y(2:end) - x(2:end) .* y(1:end-1)));
                
                % compute current density
                CD = totalCarrent / A;
            end
        
        end
    
        function obj = setCorePermeability(obj)

            points = obj.Points;
            model = obj.Model;
            pg = model.Geometry;

            faceIdes = [];
            for i = 1 : length(obj.Env.Thruster.Chambers)
                FaceID_magnetCore_right = nearestFace(pg, points.Face.("magnetCoreR" + i));
                FaceID_magnetCore_left = nearestFace(pg, points.Face.("magnetCoreL" + i));

                faceIdes = [faceIdes, [FaceID_magnetCore_right, FaceID_magnetCore_left]];
            end

            electromagneticProperties(obj.Model,"RelativePermeability",obj.CorePermeability, "Face", faceIdes); 
        end

        function obj = genMesh(obj)
            
            points = obj.Points;
            model = obj.Model;
            pg = model.Geometry;

            ides = [];

            % cores
            for i = 1 : length(obj.Env.Thruster.Chambers)
                FaceID_magnetCore_right = nearestFace(pg, points.Face.("magnetCoreR" + i));
                FaceID_magnetCore_left = nearestFace(pg, points.Face.("magnetCoreL" + i));
    
                ides = [ides, [FaceID_magnetCore_right, FaceID_magnetCore_left]];
            end
                
            % generate mesh
            generateMesh(model, "Hface", {ides, obj.DetailedLength}, ...
                "Hmin", obj.Hmin, "Hmax", obj.Hmax, "Hgrad", obj.Hgrad);
        end
    end
end

