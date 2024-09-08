function mdot = solveParaThrustForMDot(h_i, h_f, V_i, gamma, rho0, A, ...
    g0, m_i, CdS_sum, dt, Isp, V_f, m_p, tols)
% This function propagates a parachute and thrusting descent to find out
% what mdot would satisfy the final conditions.

    % Account for missing arguments
    if (nargin < 14)
        % Since the function is solved discretely, a larger y-tol is needed
        tols = [1e-10, 0.1];
    end
    if (nargin < 13)
        % Assume propellant mass is the all the mass for a worst case
        m_p = m_i;
    end

    % Initial guess of bounds, kg/s
    lowerMDot = 0;
    upperMDot = 100;
    mdot = NaN;
    params = [h_i; h_f; V_i; gamma; rho0; A; g0; m_i; CdS_sum; dt; Isp; ...
        m_p; Inf];
    while (isnan(mdot))
        try
            bounds = [lowerMDot, upperMDot];
            mdot = bisectSolver(@ParaThrustWrapper, params, bounds, ...
                tols, V_f);
            % If the code gets here, the solver worked
            if (isnan(mdot))
                % If the code gets here, the tolerances were not met even
                % though the solver converged. Thus, the bounds will not
                % change, and an infinite loop will occur
                error("An infinite loop will occur. Check tolerances.");
            end
        catch
            % The bounds guess is wrong
            lowerMDot = upperMDot;
            upperMDot = upperMDot + 100;
        end
    end

end