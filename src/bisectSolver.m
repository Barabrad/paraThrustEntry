function x_cross = bisectSolver(fxn, param_arr, x_bounds, tols, y_cross)
% This function finds where a function (given the constant parameters in an
% array) crosses y_cross where x_cross is in the finite range x_bounds.
% tols = [tolerance in x; tolerance in y], assumed to be 1e-6 for both by
% default.
% The input function fxn shall be structured like fxn(x, paramArr).
% Note that y_cross is assumed to be zero by default.
% If there is no zero found within the bounds, x_cross will be NaN.
% NOTE: The asymptote checks were added to this file after the semester.

    % Account for missing arguments
    tol_def = 1e-6;
    if (nargin < 5)
        y_cross = 0;
    end
    if (nargin < 4)
        tols = [tol_def; tol_def];
    else
        if (length(tols) == 1)
            tols = [tols; tols];
        end
        if (prod(tols) == 0)
            warning("Tolerances cannot be zero. Replacing with 1e-6...")
            zero_tols = (tols == 0);
            tols(zero_tols) = tol_def;
        end
    end

    % Establish bounds
    x_min = x_bounds(1);
    x_max = x_bounds(2);
    x_guess = (x_max + x_min)/2;
    % Establish tolerances (accounting for user shenanigans)
    x_tol = abs(tols(1));
    y_tol = abs(tols(2));
    tol_vec = char("<" + x_tol + ", " + y_tol + ">");
    if (abs(x_max - x_guess) <= x_tol)
        warning("Oh no! Cringe!");
        error("The x-tolerance should be far less than half of the range");
    end
    % Establish variables used in the loop
    x_prev = x_min - 2*x_tol; % An impossible value for the beginning

    % Determine the positive and negative ends
    dy_min = fxn(x_min, param_arr) - y_cross;
    dy_max = fxn(x_max, param_arr) - y_cross;
    % Account for asymptotes
    if or(isinf(dy_min), isnan(dy_min))
        x_min = x_min + (x_tol/2);
        dy_min = fxn(x_min, param_arr) - y_cross;
    end
    if or(isinf(dy_max), isnan(dy_max))
        x_max = x_max - (x_tol/2);
        dy_max = fxn(x_max, param_arr) - y_cross;
    end
    if (sign(dy_max*dy_min) == 1)
        % Both ends have the same sign
        error("f(x_min) and f(x_max) have the same sign.")
    end
    if (dy_min < dy_max)
        x_neg = x_min;
        x_pos = x_max;
    else
        x_neg = x_max;
        x_pos = x_min;
    end
    % Check the ends, and allow bypass of the loop if necessary
    if (abs(dy_min) <= y_tol)
        x_guess = x_min;
        dy_guess = dy_min;
        x_prev = x_guess;
    elseif (abs(dy_max) <= y_tol)
        x_guess = x_max;
        dy_guess = dy_max;
        x_prev = x_guess;
    end

    % https://i.imgflip.com/1ka4ul.jpg?a470520
    % Account for the x- and y-tolerances
    % (Note that with the way the while loop is structured, it will also be
    % broken when abs(dy_guess) <= y_tol since x_guess will equal x_prev)
    while (abs(x_guess - x_prev) > x_tol)
        % Find the difference in y
        dy_guess = fxn(x_guess, param_arr) - y_cross;
        % Prepare for next loop
        x_prev = x_guess;
        if (abs(dy_guess) > y_tol)
            if (dy_guess > 0)
                % x_guess was an overestimate
                x_pos = x_guess;
            else
                % x_guess was an underestimate
                x_neg = x_guess;
            end
            x_guess = (x_pos + x_neg)/2;
        end
    end

    if (abs(dy_guess) > y_tol)
        x_cross = NaN;
        currErrY = char(string(dy_guess));
        warningMsg = ['With a tolerance vector ', tol_vec, ...
            ', the x-tolerance was reached before the ', ...
            'y-tolerance.\n', ...
            'abs(f(x_guess) - y_tol) currently is ', currErrY, '\n', ...
            'The returned value is NaN.\n'];
        warning(warningMsg);
    else
        x_cross = x_guess;
    end

end