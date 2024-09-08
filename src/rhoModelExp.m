function rho_calc = rhoModelExp(rho_0, A, h_calc)
% This function finds density at a given altitude, and was provided in HW.
% rho0 = atmosphere density at surface, kg/m^3
% A = scale height, m^-1
% h_calc = current height, m

    rho_calc = rho_0*exp(-A*h_calc);

end