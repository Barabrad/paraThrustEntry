function Vy_end = ParaThrustWrapper(mdot, param_arr)
% This function acts as a wrapper that the bisection solver can use, with
% the parachute with thrust propagation.
    % Get parameters from array
    h_i = param_arr(1);
    h_f = param_arr(2);
    V_i = param_arr(3);
    gamma = param_arr(4);
    rho0 = param_arr(5);
    A = param_arr(6);
    g0 = param_arr(7);
    m_i = param_arr(8);
    CdS_sum = param_arr(9);
    dt = param_arr(10);
    Isp = param_arr(11);
    m_p = param_arr(12);
    t_prop = param_arr(13);
    % Run the propagation
    % Output: [x, h, Vx, Vy, Ax, Ay, t, rho, m, ranOut]
    [~, ~, ~, Vy_arr, ~, ~, ~, ~, ~, didRunOut] = ...
        propagateParaThrustEntry(h_i, h_f, V_i, gamma, rho0, A, g0, ...
            m_i, CdS_sum, dt, mdot, Isp, m_p, t_prop);
    % Get the final velocity
    if (didRunOut)
        % mdot was too high, so treat as if velocity was too low
        Vy_end = deal(-Inf);
    else
        Vy_end = Vy_arr(end);
    end
end