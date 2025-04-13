classdef SolverE

    properties % mesh properties
        DetailedLength = 0.002
        Hmin = 0.002
        Hmax = 0.04
        Hgrad = 1.3
    end

    properties
        Env
        Model
        Points
        SysGeometry
        Results
    end

    methods
        
        function obj = SolverE(env)
            obj.Env = env;
            
            % create electro statis model
            model = createpde("electromagnetic","electrostatic");
            model.VacuumPermittivity  = 8.8541878128E-12;
            obj.Model = model;
        end

        function obj = solve(obj)
             % set model geometry
            obj = setGeometry(obj);

            % set boundary condition
            obj = setBC(obj);
            electromagneticProperties(obj.Model,"RelativePermittivity",1.0); 

            % generate mesh
            obj = genMesh(obj);

            % solve problem 
            obj.Results = solve(obj.Model);
        end
    
        function plotResults(obj)
            figure;
            
            pdeplot(obj.Model, ...
                "XYData",obj.Results.ElectricPotential, ...
                "Contour","on", ...
                "FlowData", [obj.Results.ElectricField.Ex, obj.Results.ElectricField.Ey], ...
                "FaceAlpha", 0.3, ...
                "Levels",20);
            hold on
            pdegplot(obj.SysGeometry);
            xlabel("r")
            ylabel("z")
            title("electric potential and electric field")

        end

    end
    
    % helper methods
    methods (Access = private)
        
        function obj = setGeometry(obj)
        
            % create whole geometry
            [gd, ns, points] = createGeometry(obj.Env);
            obj.Points = points;
            obj.SysGeometry = decsg(gd);
            
            % create remaining space geometry
            sf = [ns{end}, '-('];
            for i = 1 : length(ns)-1
                sf = [sf, ns{i}, '+'];
            end
            sf(end) = ')';
            ns = char(ns);
            [dl,~] = decsg(gd,sf,ns');

            % set model geometry
            geometryFromEdges(obj.Model, dl);
            
        end
        
        function obj = setBC(obj)
        
            points = obj.Points;
            model = obj.Model;
            pg = model.Geometry;

            % cathode 
            EdgeID_cathode = nearestEdge(pg, points.Edge.cathodeContact);
            electromagneticBC(model,"Voltage",obj.Env.Thruster.Cathode.Voltage,"Edge",EdgeID_cathode);

            % anodes
            for i = 1 : length(obj.Env.Thruster.Chambers)
                chamber = obj.Env.Thruster.Chambers{i};

                EdgeID_anode_right = nearestEdge(pg, points.Edge.("chamberContactR" + i));
                EdgeID_anode_left = nearestEdge(pg, points.Edge.("chamberContactL" + i));

                electromagneticBC(model,"Voltage",chamber.AnodeVoltage,"Edge",[EdgeID_anode_right, EdgeID_anode_left]);
            end
        
        end
    
        function obj = genMesh(obj)
            
            points = obj.Points;
            model = obj.Model;
            pg = model.Geometry;

            ides = [];

            % cathode 
            EdgeID_cathode = nearestEdge(pg, points.Edge.cathodeContact);
            ides = [ides, EdgeID_cathode];

            % anodes
            for i = 1 : length(obj.Env.Thruster.Chambers)
                EdgeID_anode_right = nearestEdge(pg, points.Edge.("chamberContactR" + i));
                EdgeID_anode_left = nearestEdge(pg, points.Edge.("chamberContactL" + i));
                ides = [ides, [EdgeID_anode_right, EdgeID_anode_left]];
            end
                
            % generate mesh
            generateMesh(model, "Hedge", {ides, obj.DetailedLength}, ...
                "Hmin", obj.Hmin, "Hmax", obj.Hmax, "Hgrad", obj.Hgrad);
        end
    end
end

