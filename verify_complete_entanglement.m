%%  verify_complete_entanglement
%   Returns 1 if the completely entangled subspace of maximal dimension
%   as per https://arxiv.org/pdf/quant-ph/0409032 is unable to be determined
%   to be completely entangled via our method, and a value less than 1 if it is.
%   
%   This function has one required argument:
%     dims: A vector of p integers representing the dimensions
%           of the tensor factor systems.
%
%   c = verify_complete_entanglement(dims) outputs W^{1+i}_{max} of
%   the projection matrix onto the largest completely entangled subspace
%   of the multipartite system of dimensions specified by dims.

%   requires: cvx (http://cvxr.com/cvx/), QETLAB (qetlab.com)
%   author: Logan Pipes (ldpipes@mta.ca)
%   last updated: August 12, 2022

function c = verify_complete_entanglement(dims)
    if min(dims) <= 1
        c = 0;
        return
    end
    P = max_complete_entangled_proj(dims);
    parties = numel(dims);
    A = cell(2^parties,1);
    for subset=1:2^parties - 1 % All partial transposes
        A{subset} = PartialTranspose(P, find(fliplr(dec2bin(subset,parties))=='1'), dims);
    end
    A{2^parties} = P; % And identity
    B = A; % copy for modification, will be A but with the convex weighting applied

    I = iden(size(P,1),1);
    cvx_begin sdp quiet
    variable p(size(A,1))
    variable c
    minimize c
    subject to
    for k=1:size(A,1)
        B{k} = p(k)*A{k};
    end
    c*I - sum(cat(3, B{:}), 3) >= 0;
    sum(p) == 1;
    cvx_end
end
