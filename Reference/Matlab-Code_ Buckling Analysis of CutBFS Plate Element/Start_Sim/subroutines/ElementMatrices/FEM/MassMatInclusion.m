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
%  Geometry:        Structure that contains geometric- and meta-
%                   information about the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% Me:            	Elemental mass matrix
%
%-------------------------------------------------------------------------%
function [Me] = MassMatInclusion(xElem,yElem,h,rho,Geometry)

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
    % Mass matrix
    Me = Me + N*N'*dens*h*detJ*wG(i);
end

end