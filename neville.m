%-------------------------------------------------------------------------------
% File:         neville.m
% Authors:      Igor Janjic
% Description:  Implements Neville's Method.
%               Approximates a function with the given data using Neville's
%               algorithm by computing the table of Lagrange interpolation
%               approximations Q(i, j). Note that indexing starts at 1.
% Inputs:
%   x           The selected point being approximated.
%   xi          The given nodes.
%   fxi         The function values at the given nodes.
%
% Outputs:
%   Q           The Neville matrix.
%-------------------------------------------------------------------------------

function [Q] = neville(x, xi, fxi)

% Determine n and reserve space for Q.
n = length(xi) - 1;
Q = zeros(n + 1, n + 1);

% Copy nodal function values to first column of Q.
for i = 0:n
    Q(i + 1, 1) = fxi(i + 1);
end

% Determine the other columns of Q.
for i = 1:n
    for j = 1:i
        dQ =(x - xi(i - j + 1)) * Q(i + 1, j) - (x - xi(i + 1)) * Q(i, j);
        h = xi(i + 1) - xi(i - j + 1);
        Q(i + 1, j + 1) = dQ / h;
    end
end
end
