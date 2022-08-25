%%  ternary_p
%   Computes the optimal value of p such that the maximum eigenvalue of
%   p*B + (1-p)*PartialTranspose(B) is minimized.
%
%   This function has one required argument:
%     B: A square matrix.
%
%   p = ternary_p(B) outputs the optimal value of p such that the largest
%   eigenvalue of p*B + (1-p)*PartialTranspose(B) is as small as possible.
%
%   [p, ev] = ternary_p(B) outputs the optimal value of p, as well as the largest eigenvalue.
%
%   This function has one optional argument:
%     dims: A vector of the sizes of the systems B is comprised of.
%           (by default assumes dimensions are equal)
%
%   [p, ev] = ternary_p(B, dims) outputs the same as above,
%   but with dims specified to the partial transposition if the dimensions are not equal.

%   requires: QETLAB (qetlab.com)
%   author: Logan Pipes (ldpipes@mta.ca)
%   last updated: August 25, 2022

function [p,ev] = ternary_p(B, dims)
    % OPTIONAL VARIABLE DEFAULTS
    if nargin < 2 % Default second parameter
        dims = sqrt(size(B));
    end

    TOL = 0.001;

    ptB = PartialTranspose(B,2,dims);
    lo = -1000;
    hi = 1000;

    while hi-lo > TOL
        left = (2*lo+hi)/3;
        right = (lo+2*hi)/3;
        eig_left = eigs(left*B + (1-left)*ptB,1,'largestreal');
        eig_right = eigs(right*B + (1-right)*ptB,1,'largestreal');
        if eig_left > eig_right
            lo = left;
        else
            hi = right;
        end
    end
    p = (hi+lo)/2;
    ev = eigs(p*B + (1-p)*ptB,1,'largestreal');
end
