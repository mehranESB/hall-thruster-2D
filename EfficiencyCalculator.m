classdef EfficiencyCalculator

    properties (Constant = true)
        DiffusionConstant = 0.0625
        GasDensityConstant = 1.396e+23
    end

    properties
        ElectrostaticSolver
        MagnetostaticSolver
    end

    methods
    
        % constractor
        function obj = EfficiencyCalculator(solverE, solverB)

            obj.ElectrostaticSolver = solverE;
            obj.MagnetostaticSolver = solverB;
            
        end

        % get the middle point for chambers
        function [x,y] = getChamberPotins(obj)
        
            % get one the solvers
            solver = obj.ElectrostaticSolver;

            % exctract necessary info
            chambers = solver.Env.Thruster.Chambers;
            points = solver.Points;

            % initialize coordinates 
            x = zeros(1, length(chambers));
            y = zeros(1, length(chambers));

            for i = 1 : length(chambers)

                anodeEdge = points.Edge.("chamberContactL" + i);

                x(i) = anodeEdge(1);
                y(i) = 0;

            end
            
        end

        % interpolate electric field 
        function [Ex, Ey] = interpE(obj, x, y)
            E = obj.ElectrostaticSolver.Results.interpolateElectricField(x,y);
            
            Ex = E.Ex;
            Ey = E.Ey;
        end

        % interpolate magnetic field 
        function [Bx, By] = interpB(obj, x, y)
            H = obj.MagnetostaticSolver.Results.interpolateMagneticField(x,y);
            
            mu = 1.2566370614E-6;
            
            Bx = mu * H.Hx;
            By = mu * H.Hy;
        end

        % get diffusion coef for chambers 
        function D = getDiffusionBohmCoef(obj, temprature)
            
            % get chamber points 
            [x,y] = getChamberPotins(obj);
            
            % compute diffusion coefficients 
            D = zeros(1, length(x));
            for i = 1 : length(x)
        
                % compute magnetic field in radious direction
                [Bx, ~] = interpB(obj, x(i), y(i));

                % compute bohm diffusion coefficient
                D(i) = obj.DiffusionConstant * temprature / abs(Bx);

            end

        end

        % get gas density upper bound 
        function n = getGasDensityUpperBound(obj, efficiency, ionVelocity)
            
            % get chamber points 
            [x,y] = getChamberPotins(obj);
            
            % compute gas density ub
            n = zeros(1, length(x));
            for i = 1 : length(x)
        
                % compute magnetic field in radious direction
                [Bx, ~] = interpB(obj, x(i), y(i));

                % compute electric field in z direction
                [~, Ey] = interpE(obj, x(i), y(i));

                % compute density
                n(i) = obj.GasDensityConstant * (Bx^2) * efficiency * ionVelocity / abs(Ey);

            end

        end

    end

end