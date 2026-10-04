%% Bogner-Fox-Schmit shape functions
%-------------------------------------------------------------------------%
%
% (c) 02/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the function vales of the second order derivative of the Bogner-
% Fox-Schmit shape functions with respect to xi (d^2N/dxi^2).
% (C1-continuous plate element)
%
%-------------------------------------------------------------------------%
%% Input
%
% (xi,eta):     Local coordinates (have to be provided as vectors)
%
%-------------------------------------------------------------------------%
%% Output
%
% N:            Matrix of shape functions evaluated at (xi,eta)
%
%-------------------------------------------------------------------------%
function [N] = SF_BFS_dxi2(xi,eta)

% Check whether the coordinates are given as row or column vectors
if ~isrow(xi)
    xi  = reshape(xi,1,[]);
    eta = reshape(eta,1,[]);
end

% Matrix of shape functions
% The first 4 rows correspond to the degrees of freedom associated with
% node 1 --> (1) Transverse displacement, (2) First order derivative w.r.t.
% xi, (3) First order derivative w.r.t. eta, (4) Mixed second order
% derivative of the deflection.
% Row 5  to  8 correspond to the degrees of freedom associated with node 2.
% Row 9  to 12 correspond to the degrees of freedom associated with node 3.
% Row 13 to 16 correspond to the degrees of freedom associated with node 4.
N = [       (3*xi.*(eta - 1).^2.*(eta + 2))/8; ...
      ((3*xi - 1).*(eta - 1).^2.*(eta + 2))/8; ...
            (3*xi.*(eta - 1).^2.*(eta + 1))/8; ...
      ((3*xi - 1).*(eta - 1).^2.*(eta + 1))/8; ...
           -(3*xi.*(eta - 1).^2.*(eta + 2))/8; ...
      ((3*xi + 1).*(eta - 1).^2.*(eta + 2))/8; ...
           -(3*xi.*(eta - 1).^2.*(eta + 1))/8; ...
      ((3*xi + 1).*(eta - 1).^2.*(eta + 1))/8; ...
            (3*xi.*(eta + 1).^2.*(eta - 2))/8; ...
     -((3*xi + 1).*(eta + 1).^2.*(eta - 2))/8; ...
           -(3*xi.*(eta - 1).*(eta + 1).^2)/8; ...
      ((3*xi + 1).*(eta - 1).*(eta + 1).^2)/8; ...
           -(3*xi.*(eta + 1).^2.*(eta - 2))/8; ...
     -((3*xi - 1).*(eta + 1).^2.*(eta - 2))/8; ...
            (3*xi.*(eta - 1).*(eta + 1).^2)/8; ...
      ((3*xi - 1).*(eta - 1).*(eta + 1).^2)/8];

end