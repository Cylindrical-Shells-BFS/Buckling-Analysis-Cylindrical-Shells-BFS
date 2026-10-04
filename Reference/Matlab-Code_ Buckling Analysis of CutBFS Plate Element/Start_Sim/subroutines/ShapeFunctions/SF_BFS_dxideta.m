%% Bogner-Fox-Schmit shape functions
%-------------------------------------------------------------------------%
%
% (c) 02/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the function vales of the mixed second order derivative 
% (d^2N/dxi deta) of the Bogner-Fox-Schmit shape functions.
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
function [N] = SF_BFS_dxideta(xi,eta)

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
N = [                     (9*(eta.^2 - 1).*(xi.^2 - 1))/16; ...
              -(3*(eta.^2 - 1).*(- 3*xi.^2 + 2*xi + 1))/16; ...
             -(3*(xi.^2 - 1).*(- 3*eta.^2 + 2*eta + 1))/16; ...
     ((- 3*eta.^2 + 2*eta + 1).*(- 3*xi.^2 + 2*xi + 1))/16; ...
                         -(9*(eta.^2 - 1).*(xi.^2 - 1))/16; ...
                 (3*(eta.^2 - 1).*(3*xi.^2 + 2*xi - 1))/16; ...
              (3*(xi.^2 - 1).*(- 3*eta.^2 + 2*eta + 1))/16; ...
      -((- 3*eta.^2 + 2*eta + 1).*(3*xi.^2 + 2*xi - 1))/16; ...
                          (9*(eta.^2 - 1).*(xi.^2 - 1))/16; ...
                -(3*(eta.^2 - 1).*(3*xi.^2 + 2*xi - 1))/16; ...
               -(3*(xi.^2 - 1).*(3*eta.^2 + 2*eta - 1))/16; ...
         ((3*eta.^2 + 2*eta - 1).*(3*xi.^2 + 2*xi - 1))/16; ...
                         -(9*(eta.^2 - 1).*(xi.^2 - 1))/16; ...
               (3*(eta.^2 - 1).*(- 3*xi.^2 + 2*xi + 1))/16; ...
                (3*(xi.^2 - 1).*(3*eta.^2 + 2*eta - 1))/16; ...
      -((3*eta.^2 + 2*eta - 1).*(- 3*xi.^2 + 2*xi + 1))/16];

end