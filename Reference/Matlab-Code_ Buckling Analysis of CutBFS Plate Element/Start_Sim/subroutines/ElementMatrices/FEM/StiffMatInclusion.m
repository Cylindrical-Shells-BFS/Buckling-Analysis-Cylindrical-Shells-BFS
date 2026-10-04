%% Stiffness matrix: Bogner-Fox-Schmit plate element
%-------------------------------------------------------------------------%
%
% (c) 02/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the stiffness matrix of a Bogner-Fox-Schmit plate element.
%
%-------------------------------------------------------------------------%
%% Input
%
% (xElem,yElem):    Global coordinates of the nodes
%  h:               Plate thickness
%  D:               Plate elasticity matrix
%  Geometry:        Structure that contains geometric- and meta-
%                   information about the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% Ke:            	Elemental stiffness matrix
%
%-------------------------------------------------------------------------%
function [Ke] = StiffMatInclusion(xElem,yElem,h,D,Geometry)

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

% Stiffness matrix
Ke = zeros(16,16);

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
        C = D{Geometry.material{IdxDomain}};
    elseif Hole == 2
        % Physical domain
     	C = D{Geometry.base_mat};
    end
elseif FCMval == 0
    % The centroid of the element is not located in one of the geometric
    % regions defined by the variable Geometry. Therefore it is located in
    % the base mesh.
    C = D{Geometry.base_mat};
end

% Loop over all integration points
for i = 1:nIP
    % Jacobian matrix -- J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
    J = JacobianMatrix(xiG(i),etaG(i),xElem,yElem);
    % Determinant of the Jacobi matrix
    detJ = det(J);
    % Inverse of the Jacobian matrix
    Jinv = [J(2,2), -J(1,2); -J(2,1), J(1,1)]/det(J);
    % Matrix of the second order derivatives of the shape functions w.r.t.
    % the local coordinates
    Bl = [ SF_BFS_dxi2(xiG(i),etaG(i))'   ; ...
           SF_BFS_deta2(xiG(i),etaG(i))'  ; ...
          2*SF_BFS_dxideta(xiG(i),etaG(i))'];
    % Matrix of the second order derivatives of the shape functions w.r.t.
    % the global coordinates
    Bg = [  Bl(1,:)*Jinv(1,1)^2         + ...
          2*Bl(3,:)*Jinv(1,1)*Jinv(1,2) + ...
            Bl(2,:)*Jinv(1,2)^2         ; ...
            Bl(1,:)*Jinv(2,1)^2         + ...
          2*Bl(3,:)*Jinv(2,1)*Jinv(2,2) + ...
            Bl(2,:)*Jinv(2,2)^2         ; ...
            Bl(1,:)*Jinv(1,1)*Jinv(2,1) + ...
            Bl(3,:)*Jinv(2,1)*Jinv(1,2) + ...
            Bl(3,:)*Jinv(1,1)*Jinv(2,2) + ...
            Bl(2,:)*Jinv(1,2)*Jinv(2,2)];
    % Stiffness matrix
    Ke = Ke + Bg'*C*Bg*h*detJ*wG(i);
end

end