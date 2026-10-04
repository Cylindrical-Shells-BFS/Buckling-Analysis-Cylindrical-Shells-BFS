%% Gravity load
%-------------------------------------------------------------------------%
%
% (c) 03/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the gravity load (self-weight).
%
%-------------------------------------------------------------------------%
%% Input
%
% (xElem,yElem):    Global coordinates of the nodes
%  h:               Plate thickness
%  rho:             Mass density
%  g:               Acceleration due to gravity
%
%-------------------------------------------------------------------------%
%% Output
%
% Fe:            	Elemental load vector due to gravity
%
%-------------------------------------------------------------------------%
function [Fe] = Gravity(xElem,yElem,h,rho,g)

% Compute Fe
Fe = zeros(16,1);

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

% Loop over all integration points
for i = 1:nIP
    % Compute shape functions
    [N] = SF_BFS(xiG(i),etaG(i));
    % Jacobian matrix -- J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
    J = JacobianMatrix(xiG(i),etaG(i),xElem,yElem);
    % Determinant of the Jacobi matrix
    detJ = det(J);
    % Gravity load
    Fe = Fe + N*rho*g*h*detJ*wG(i);
end
fmax = max(abs(Fe));
Fe(abs(Fe)<fmax*1e-10) = 0;

end