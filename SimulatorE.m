classdef SimulatorE
    properties
        q = -1.602e-19;     % Electron charge (Coulombs)
        m = 9.109e-31;      % Electron mass (kg)
        % Optionally, you may store a default time step for plotting or postprocessing
        time_step = 0.001;  % time step in seconds (if needed for other purposes)
    end

    properties
        EResults  % Data for the electric field
        BResults  % (Unused here, but available for extensions)
    end

    methods
        function obj = SimulatorE(eres)
            obj.EResults = eres;
        end

        %% Simulation using ode45
        function [t, positions, velocities] = simulateOde45(obj, x0, y0, vx0, vy0, tspan)
            % Inputs:
            %   x0, y0    - Initial position (meters)
            %   vx0, vy0  - Initial velocity (m/s)
            %   tspan     - [t_start, t_end] (seconds)
            %
            % Outputs:
            %   t         - Time vector from ode45
            %   positions - [x, y] positions at each time step
            %   velocities- [vx, vy] velocities at each time step

            % Initial state vector: [x; y; vx; vy]
            Y0 = [x0; y0; vx0; vy0];
            
            % (Optional) Define event function to stop integration if field becomes NaN
            options = odeset('Events', @(t, Y) obj.fieldExitEvent(t, Y));
            
            % Call ode45 to solve the ODE system
            [t, Y, te, ye, ie] = ode45(@(t, Y) obj.electronODE(t, Y), tspan, Y0, options);
            
            % Extract positions and velocities from the state vector
            positions = Y(:, 1:2);
            velocities = Y(:, 3:4);
            
            % Plot the trajectory
            plot(positions(:, 1), positions(:, 2), 'm');
        end

        %% ODE function defining the electron dynamics
        function dYdt = electronODE(obj, t, Y)
            % Unpack the state vector
            x = Y(1);
            y = Y(2);
            vx = Y(3);
            vy = Y(4);
            
            % Interpolate the electric field at the current position
            Eintrp = interpolateElectricField(obj.EResults, x, y);
            Ex = Eintrp.Ex;
            Ey = Eintrp.Ey;
            
            % Check if the electric field is defined at (x,y)
            if isnan(Ex) || isnan(Ey)
                % If not, you can either set the acceleration to zero or choose to stop integration
                % Here we set the acceleration to zero.
                ax = 0;
                ay = 0;
            else
                % Compute acceleration from F = qE and a = F/m
                ax = -0.17587 * Ex;
                ay = -0.17587 * Ey;
            end
            
            % Return the derivative of the state vector
            dYdt = [vx; vy; ax; ay];
        end

        %% (Optional) Event function to stop integration if electric field is not defined
        function [value, isterminal, direction] = fieldExitEvent(obj, t, Y)
            % This event function stops the integration when the interpolated
            % electric field returns NaN.
            x = Y(1);
            y = Y(2);
            Eintrp = interpolateElectricField(obj.EResults, x, y);
            % When the field is defined, we want the event value to be positive.
            % When the field is NaN, we set value to zero to trigger the event.
            if isnan(Eintrp.Ex) || isnan(Eintrp.Ey)
                value = 0;   % Event triggered: stop integration.
            else
                value = 1;   % Continue integration.
            end
            isterminal = 1;  % Stop the integration when the event occurs.
            direction = 0;   % The zero can be approached from any direction.
        end

    end
end
