function p  = pp2msspoly(t,s)
% Piecewise polynomial (output of spline) to msspoly
% Input:
%   t     -- 1-by-1 msspoly (variable of spline)
%   s     -- pp object of n-by-1 with N pieces.
% Output:
%   p     -- n-by-N msspoly in t. Each column a piece
% Copied from 'spot' sim_pp2p by Yi-Hsuan Chen 08/28/2023

N   = s.pieces;
o   = s.order;
n   = s.dim;

if nargin < 3, pieces = 1:N; end
if size(pieces,1) ~= 1, error('pieces must be a row'); end
if size(n) ~= [1 1], error('can only transform n-by-1 pieces'); end

M   = length(pieces);
I   = repmat((1:n)',1,M)+ n*repmat(pieces-1,n,1);
m   = monomials(t,0:o-1);
p   = reshape(s.coefs(I,:)*m(o:-1:1),n,M);
end