%%  bound_mu
%   Computes the minimum or maximum real scalar c such that c(1+i)
%   is in the numerical range of A.
%
%   Note that the paper classifies W^{1+i}_{max}(A)
%   as bound_mu(A + 1i*PartialTranspose(A)), not as bound_mu(A).
%   
%   This function has one required argument:
%     A: A square matrix.
%
%   c = bound_mu(A) outputs the largest real scalar c such that c(1+i)
%   is in the numerical range of A.
%
%   This function has one optional argument:
%     MINMAX: A string (from the set {'min', 'minimum', 'max', 'maximum'})
%             or a boolean (0 for min, 1 for max) denoting whether the bound
%             should be an upper or lower bound. (default 1)
%
%   c = bound_mu(A, MINMAX) outputs the largest real scalar c
%   if MINMAX is 1, 'max', or 'maximum', or the smallest
%   if MINMAX is 0, 'min', or 'minimum', such that c(1+i) is in the numerical range of A.

%   author: Logan Pipes (ldpipes@mta.ca)
%   last updated: August 5, 2022

function c = bound_mu(A, MINMAX)
    % OPTIONAL VARIABLE DEFAULTS
    if nargin < 2 % Default second parameter
        MINMAX = 1;
    elseif strcmp(MINMAX, 'min') || strcmp(MINMAX, 'minimum') || (isnumeric(MINMAX) && MINMAX == 0) % Sanitize
        MINMAX = 0;
    else
        MINMAX = 1;
    end

    % CONSTANTS
    e = exp(1);
    THETA_TOLERANCE = 1e-15; % Approximately 52 iterations

    lower = pi/4; % Set appropriate bounds for computing W_max or W_min
    if MINMAX == 1
        upper = -3*pi/4;
    else
        upper = 5*pi/4;
    end

    % CHECK STARTING BOUNDS
    R = e^(1i*pi/4) * A; % Rotate
    H = (R + R')/2; % Get the Hermitian part
    [vecs,~] = eigs(H,2,'bothendsreal'); % Get the largest and smallest eigenvalue and corresponding eigenvector pairs
    l_vec = vecs(:,2);
    u_vec = vecs(:,1);
    l_vec = l_vec/norm(l_vec); % Normalize the eigenvectors
    u_vec = u_vec/norm(u_vec);
    l_pt = complex(l_vec' * A * l_vec); % Get the corresponding point in the unrotated numerical range (they're boundary points)
    u_pt = complex(u_vec' * A * u_vec);

    if imag(l_pt) == real(l_pt)
        % Tangent on bottom
        c = real(l_pt);
    elseif imag(u_pt) == real(u_pt)
        % Tangent on top
        c = real(u_pt);
    elseif imag(l_pt) > real(l_pt) || imag(u_pt) < real(u_pt)
        % Disjoint
        if MINMAX == 1
            c = -Inf;
        else
            c = Inf;
        end
    else

    % INTERSECTION
        while (abs(upper-lower) > THETA_TOLERANCE)
            mid = (lower+upper)/2;
            R = e^(1i*mid) * A; % Rotate
            H = (R + R')/2; % Get the Hermitian part
            [x,~] = eigs(H,1,'largestreal'); % Get the eigenvector corresponding to the largest eigenvalue
            x = x/norm(x); % Normalize it
            p = x' * A * x; % Get the corresponding point in the unrotated numerical range (a boundary point)
            if imag(p) <= real(p) % Binary search condition
                lower = mid;
            else
                upper = mid;
            end
        end

        R = e^(1i*lower) * A; % Get vertices determined by binary search for intersection approximation
        H = (R + R')/2;
        [x,~] = eigs(H,1,'largestreal');
        x = x/norm(x);
        l = complex(x' * A * x);

        R = e^(1i*upper) * A;
        H = (R + R')/2;
        [x,~] = eigs(H,1,'largestreal');
        x = x/norm(x);
        u = complex(x' * A * x);

        outer_upper = (real(u) + imag(u) + (imag(u)-real(u))*(cot(upper - pi/4)))/2; % Get outer approximations
        outer_lower = (real(l) + imag(l) + (real(l)-imag(l))*(cot(pi/4 - lower)))/2;
        if MINMAX == 1 % Pick the best outer bound (i.e. the innermost one)
            outer_approx = min(outer_upper, outer_lower);
        else
            outer_approx = max(outer_upper, outer_lower);
        end
        c = outer_approx; % Assign it as the value to output (outer bound will always be valid, inner bound will not ever)

        % This determines the point on the line from u to l with imag(p) == real(p)
        % but beacause u and l can be very close, it's innacurate due to division by an infinitesimal
        %inner_approx = (real(l)*imag(u) - imag(l)*real(u))/(real(l)-imag(l)-real(u)+imag(u));
        %new_pt = [real(l)-real(u), real(u)-real(l); imag(l)-imag(u), real(u)-real(l)]\[0; imag(l)-imag(u)-imag(u)*(real(l)-real(u))]
        %inner_approx = new_pt(1);
        %disp("Error within " + abs(outer_approx - inner_approx));
    end
end
