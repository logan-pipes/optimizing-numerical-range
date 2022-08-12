%%  max_complete_entangled_proj
%   Returns the projection matrix to the largest completely entangled subspace
%   of M_dims(0) \otimes M_dims(1) \otimes \cdots \otimes M_dims(p)
%   as per https://arxiv.org/pdf/quant-ph/0409032.
%   
%   This function has one required argument:
%     dims: A vector of p integers representing the dimensions
%           of the tensor factor systems.
%
%   P = max_complete_entangled_proj(dims) outputs a sparse matrix P
%   that is the projection matrix onto the largest completely entangled subspace
%   of the tensor space of the given dimensions.

%   requires: QETLAB (qetlab.com)
%   author: Logan Pipes (ldpipes@mta.ca)
%   last updated: August 11, 2022

function P = max_complete_entangled_proj(dims)
    p = numel(dims); % length of dims
    prods = ones(1,p);
    for jj=p-1:-1:1 % backwards from second last index
        prods(jj) = prods(jj+1)*dims(jj+1);
    end % products(j) is the product of dims(j+1:)
    anti_diags = zeros(1,prod(dims));
    index_vec = zeros(p,1); % tensor indices
    for iter=1:prod(dims)
        anti_diags(prods * index_vec + 1) = sum(index_vec)+p;
        index_vec = update_odometer(index_vec, dims);
    end % anti_diags holds the "anti-diagonal index"---the sum of the indices of the SBV tensored together to get that index

    % fill V with spanning vectors
    sz = round(norm(sum(anti_diags' == p:sum(dims)))^2); % sum of squares of number of entries on a given diagonal
    row_indices = zeros(1,2*sz);
    iter = 1;
    for ssum = p:sum(dims)
        indices = find(anti_diags == ssum);
        for ii=indices
            for jj=indices
                row_indices(iter:iter+1) = [ii,jj];
                iter = iter+2;
            end
        end
    end
    V = sparse(row_indices, reshape([1:sz; 1:sz], 1,sz*2), reshape([ones(1,sz); -ones(1,sz)], 1,sz*2), prod(dims),sz);

    S = sporth(V); % Trim to linearly independent set
    P = S*S'; % generate projection to desired subspace
    % verify rank one avoidance with our method (SDP with all partial transposes)
end
