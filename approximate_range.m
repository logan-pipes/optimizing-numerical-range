%%  approximate_range
%   Displays inner and outer bounds on the numerical range of A,
%   and generates random points on the interior.
%
%   This function has one required argument:
%     A: A square matrix.
%
%   This function has two optional arguments:
%     MAX_BOUND_ITERATIONS: An integer specifying how many vertices the
%                           boundary polygons should contain.
%     MAX_POINT_ITERATIONS: An integer specifying how many interior points
%                           to generate.
%

%   author: Logan Pipes (ldpipes@mta.ca)
%   last updated: August 5, 2022

function approximate_range(A, MAX_BOUND_ITERATIONS, MAX_POINT_ITERATIONS)
    % OPTIONAL VARIABLE DEFAULTS
    if nargin < 3
        MAX_POINT_ITERATIONS = 20000;
    end
    if nargin < 2
        MAX_BOUND_ITERATIONS = 500;
    end

    % CONSTANTS
    e = exp(1);
    n = size(A,1);

    % COMPUTE INNER AND OUTER BOUND ON NUMERICAL RANGE
    if MAX_BOUND_ITERATIONS ~= 0
        inner_boundary = complex(zeros(MAX_BOUND_ITERATIONS+1,1));
        outer_boundary = complex(zeros(MAX_BOUND_ITERATIONS+1,1));
        cdif = cos(-2*pi/MAX_BOUND_ITERATIONS);
        sdif = sin(-2*pi/MAX_BOUND_ITERATIONS);
        
        % FIRST ITERATION PREP
        H = (A + A')/2;
        [x,lambda_phi] = eigs(H,1,'largestreal');
        x = x/norm(x);
        p = x' * A * x;
        inner_boundary(1) = p;
        
        % PRIMARY COMPUTATION
        for iteration=1:MAX_BOUND_ITERATIONS
            theta = 2*pi*iteration/MAX_BOUND_ITERATIONS; % Calculate theta
            R = e^(1i*theta) * A; % Rotate the matrx
            H = (R + R')/2; % Get the Hermitian part
            [x, lambda_theta] = eigs(H,1,'largestreal'); % Get the largest eigenvalue and corresponding eigenvector
            x = x/norm(x); % Normalize it
            p = x' * A * x; % Get the corresponding point in the numerical range (it's a boundary point)
            % The next line calculates the intersection of the tangent line to the numerical range
            % at the current angle with the tangent line at the last angle
            pt = e^(-1i*theta)*(lambda_theta + 1i*(cdif*lambda_theta - lambda_phi)/sdif);
            inner_boundary(iteration+1) = p; % Store the data for plotting
            outer_boundary(iteration) = pt;
            lambda_phi = lambda_theta; % Update previous angle's largest eigenvector
        end
        outer_boundary(MAX_BOUND_ITERATIONS+1)=outer_boundary(1); % Close the polygon
    end

    % GENERATE INTERIOR POINTS
    if MAX_POINT_ITERATIONS ~= 0
        interior = zeros(1, MAX_POINT_ITERATIONS);
        for iteration=1:MAX_POINT_ITERATIONS
            v = 2*rand(n,1)-1 + 2i*rand(n,1)-1i; % Generate a random vector
            v = v/norm(v); % Normalize it
            z = v' * A * v; % Get the corresponding point in the numerical range
            interior(iteration) = z; % Store the data for plotting
        end
    end


    % OUTPUT
    holdstate = ishold;
    if MAX_BOUND_ITERATIONS ~= 0
        plot(real(inner_boundary), imag(inner_boundary));
        hold on
        plot(real(outer_boundary), imag(outer_boundary));
    end
    if MAX_POINT_ITERATIONS ~= 0
        plot(interior,'.');
    end
    % Restore the hold state to what it was previously
    if holdstate == 1
        hold on
    else
        hold off
    end
end
