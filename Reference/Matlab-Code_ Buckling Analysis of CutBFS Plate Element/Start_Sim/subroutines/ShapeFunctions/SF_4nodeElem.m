%% 4-node element shape functions
%-------------------------------------------------------------------------%
%
% (c) 03/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the function vales of the of the 4-node element shape functions.
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
function [N] = SF_4nodeElem(xi,eta)


% Check whether the coordinates are given as row or column vectors
if ~isrow(xi)
    xi  = reshape(xi,1,[]);
    eta = reshape(eta,1,[]);
end

% Matrix of shape functions
N = [1/4*(1-xi).*(1-eta);...
     1/4*(1+xi).*(1-eta);...
     1/4*(1+xi).*(1+eta);...
     1/4*(1-xi).*(1+eta)];

end