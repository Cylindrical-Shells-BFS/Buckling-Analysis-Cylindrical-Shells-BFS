%% Mass matrix: Bogner-Fox-Schmit plate element
%-------------------------------------------------------------------------%
%
% (c) 02/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the (consistent) mass matrix of a Bogner-Fox-Schmit plate
% element.
%
%-------------------------------------------------------------------------%
%% Input
%
% (xElem,yElem):    Global coordinates of the nodes
%  h:               Plate thickness
%  rho:             Mass density
%
%-------------------------------------------------------------------------%
%% Output
%
% Me:            	Elemental mass matrix
%
%-------------------------------------------------------------------------%
function [Me] = MassMat(xElem,yElem,h,rho)

% Numerical integration - Gauss-Legendre quadrature rule
% Gauss points and weights - 1d
[xG,wG] = GaussQuad(4,-1,1);
% Gauss points - 2d
[xiMat,etaMat] = meshgrid(xG,xG);
xiG  = reshape(xiMat,[],1);
etaG = reshape(etaMat,[],1);
% Number of integration points
nIP = numel(xiG);
% Weights - 2d
[wxiMat,wetaMat] = meshgrid(wG,wG);
wxiG  = reshape(wxiMat,[],1);
wetaG = reshape(wetaMat,[],1);
wG    = wxiG.*wetaG;

% Mass matrix
Me = zeros(16,16);
% Loop over all integration points
for i = 1:nIP
    % Compute shape functions
    [N] = SF_BFS(xiG(i),etaG(i));
    % Jacobian matrix -- J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
    J = JacobianMatrix(xiG(i),etaG(i),xElem,yElem);
    % Determinant of the Jacobi matrix
    detJ = det(J);
    % Mass matrix
    Me = Me + N*N'*rho*h*detJ*wG(i);
end

end