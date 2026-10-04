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
%  Geometry:        Structure that contains geometric- and meta-
%                   information about the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% Fe:            	Elemental load vector due to gravity
%
%-------------------------------------------------------------------------%
function [Fe] = GravityInclusion(xElem,yElem,h,rho,Geometry)

% Acceleration
g = Geometry.g;

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

% Check which material has to be used
xc = 1/4*sum(xElem);
yc = 1/4*sum(yElem);
[FCMval,IdxDomain,Hole] = PointInMembershipTest(xc,yc,Geometry);

% Material
if FCMval == 1
    % The centroid of the cell is located in one of the geometric regions
    % defined by the variable Geometry.
    if Hole == 0
        % Inclusion
        dens = rho{Geometry.material{IdxDomain}};
    elseif Hole == 2
        % Physical domain
     	dens = rho{Geometry.base_mat};
    end
elseif FCMval == 0
    % The centroid of the element is not located in one of the geometric
    % regions defined by the variable Geometry. Therefore it is located in
    % the base mesh.
    dens = rho{Geometry.base_mat};
end

% Loop over all integration points
for i = 1:nIP
    % Compute shape functions
    [N] = SF_BFS(xiG(i),etaG(i));
    % Jacobian matrix -- J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
    J = JacobianMatrix(xiG(i),etaG(i),xElem,yElem);
    % Determinant of the Jacobi matrix
    detJ = det(J);
    % Gravity load
    Fe = Fe + N*dens*g*h*detJ*wG(i);
end
fmax = max(abs(Fe));
Fe(abs(Fe)<fmax*1e-10) = 0;

end