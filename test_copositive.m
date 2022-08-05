%%  test_copositive
%   Uses semidefinite program to test for copositivity of the input matrix.
%   
%   This function has one required argument:
%     W: A real symmetric matrix.
%
%   c = test_copositive(W) outputs a positive number if W is not copositive.
%   If W is copositive, it may output any real number.

%   requires: cvx (http://cvxr.com/cvx/), QETLAB (qetlab.com)
%   author: Logan Pipes (ldpipes@mta.ca)
%   last updated: August 5, 2022

function c = test_copositive(W)
    n = size(W,1);
    E = zeros(n^2, n);
    I = eye(n);
    for j=1:n
        E((j-1)*n+1:j*n, j) = I(:,j);
    end
    W = E*W*E';
    
    % All 8 size preserving
    A = {W; PartialTranspose(W); PartialTranspose(Realignment(W)); Realignment(PartialTranspose(Realignment(W))); PartialTranspose(Realignment(PartialTranspose(W,1))); PartialTranspose((Realignment(W'))); Swap(W); PartialTranspose(Swap(W))}; 

    B = A;

    cvx_begin sdp quiet
    variable p(size(A,1))
    variable c
    maximize c
    subject to
    for k=1:size(A,1)
        B{k} = p(k)*A{k};
    end
    sum(cat(3, B{:}), 3) - c*eye(size(W)) >= 0;
    sum(p) == 1;
    cvx_end
end
