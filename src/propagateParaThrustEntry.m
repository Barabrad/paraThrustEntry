function [x, h, Vx, Vy, Ax, Ay, t, rho, m, ranOut] = ...
    propagateParaThrustEntry(h_i, h_f, V_i, gamma, rho0, A, g0, m_i, ...
    CdS_sum, dt, mdot, Isp, m_p, t_prop)
% This function finds the height and velocity for a parachute + thrusting
% entry.
% h_i = inital height (m)
% h_f = final height (m)
% V_i = inital velocity (m/s)
% gamma = flight path angle, relative to horizon (º)
% rho0 = surface density
% A = 1/(scale height), m^-1
% g0 = planet gravitational acceleration, m/s^2
% m_i = EV initial mass, kg
% CdS_sum = sum(Cd.*S) for each body (EV and chutes)
% dt = propagator time step, s
% mdot = propellant usage rate, kg/s
% Isp = specific impulse, s
% m_p = propellant mass, kg
% - Note: m_p is assumed to be 90% of the initial mass by default.
% t_prop = time interval to stop propagating after, s
% - Note: t_prop is Inf by default, which will result in propagation ending
%   after the final height is reached.

    % Account for missing arguments
    if (nargin < 14)
        t_prop = Inf;
    end
    if (nargin < 13)
        m_p = 0.9*m_i;
    end
    % Make gamma positive for convenience
    gamma = abs(gamma);
    % Account for t_prop being smaller than dt
    if (dt > t_prop)
        dt = t_prop/2; % Run at least 2 steps
    end

    % Constants
    m_cutoff = m_i - m_p;
    g_E = 9.81; % m/s^2
    % Pre-allocate arrays, assuming the entry will take at most 60 seconds
    len_preAlloc = ceil(60/dt);
    [h, x, Vx, Vy, Ax, Ay, t, rho, m] = deal(nan(len_preAlloc, 1));
    % Make iterators, knowing we will fill in the first slot
    h_curr = h_i;
    x_curr = 0;
    Vx_curr = V_i*cosd(gamma);
    Vy_curr = V_i*sind(gamma);
    Ax_curr = g0*cosd(gamma);
    Ay_curr = g0*sind(gamma);
    m_curr = m_i;
    t_curr = 0; % No propagation yet
    index = 1;
    % Fill in the first slots
    h(index) = h_curr;
    x(index) = x_curr;
    Vx(index) = Vx_curr;
    Vy(index) = Vy_curr;
    Ax(index) = Ax_curr;
    Ay(index) = Ay_curr;
    t(index) = t_curr;
    m(index) = m_curr;
    % Propagate with a while loop that accounts for time and height
    ranOut = false;
    while and((t_prop - t_curr > dt/2), (h_curr > h_f))
        % Calculate next values
        rho_calc = rhoModelExp(rho0, A, h_curr);
        % Simplify calculations
        RCDSdiv2M = rho_calc*CdS_sum/(2*m_curr);
        if or(isinf(RCDSdiv2M), isnan(RCDSdiv2M))
            warning("Division error occurred when mass hit zero.");
            ranOut = true; % In case the cutoff mass is negative
            break; % Stop before things get weird with Infs and NaNs
        end
        a_Dx = sign(Vx_curr)*RCDSdiv2M*(Vx_curr^2);
        a_Dy = sign(Vy_curr)*RCDSdiv2M*(Vy_curr^2);
        if (m_curr > m_cutoff)
            a_T = mdot*g_E*Isp/m_curr;
            a_Tx = a_T*cosd(gamma);
            a_Ty = a_T*sind(gamma);
            m_curr = m_curr - mdot*dt;
        else
            [a_Tx, a_Ty] = deal(0); % No thrust if no propellant
            ranOut = true;
            m_curr = m_cutoff;
        end
        % In the reference frame, positive is down and right
        % Thrust and drag act opposite of falling velocity
        Ax_curr = -a_Dx - a_Tx;
        Ay_curr = g0 - a_Dy - a_Ty;
        Vx_next = Vx_curr + Ax_curr*dt;
        Vy_next = Vy_curr + Ay_curr*dt;
        x_next = x_curr + Vx_curr*dt;
        h_next = h_curr - Vy_curr*dt;
        if (h_next >= h_i)
            % There is way too much thrust, and there's going to be another
            % free-fall from a higher altitude
            %warning("Too much thrust. Exiting loop...");
            ranOut = true;
            break;
        end
        % Update arrays
        rho(index) = rho_calc;
        t_curr = t_curr + dt;
        index = index + 1;
        x(index) = x_next;
        h(index) = h_next;
        Vx(index) = Vx_next;
        Vy(index) = Vy_next;
        Ax(index) = Ax_curr;
        Ay(index) = Ay_curr;
        t(index) = t_curr;
        m(index) = m_curr;
        % Update current values
        x_curr = x_next;
        h_curr = h_next;
        Vx_curr = Vx_next;
        Vy_curr = Vy_next;
        gamma = atan2d(Vy_curr, Vx_curr); % The angle should approach 90º
        if (gamma < 0)
            % Note that the rocket is now moving upwards
            gamma = abs(gamma); % Other parts of the code will catch this
        end
        if (isnan(gamma))
            % This case should not have been reached
            warning("gamma became NaN. Investigate this case.");
        end
    end

    % Remove any leftover NaN values from the arrays
    if (index < len_preAlloc)
        x = x(1:index);
        h = h(1:index);
        Vx = Vx(1:index);
        Vy = Vy(1:index);
        Ax = Ax(1:index);
        Ay = Ay(1:index);
        t = t(1:index);
        m = m(1:index);
        rho = rho(1:index);
    end

    % Find the most recent rho value (not in array yet)
    rho(index) = rhoModelExp(rho0, A, h_curr);

end