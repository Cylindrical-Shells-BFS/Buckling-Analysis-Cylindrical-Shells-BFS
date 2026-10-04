%% Stiffness matrix (FCM): Bogner-Fox-Schmit plate element
%-------------------------------------------------------------------------%
%
% (c) 04/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the stiffness matrix of a Bogner-Fox-Schmit plate element.
%
%-------------------------------------------------------------------------%
%% Input
%
% (xElem,yElem):    Global coordinates of the nodes
%  h:               Plate thickness
%  D:               Plate elasticity matrix
%  QT:              Quadtree data
%  alpha:           Value of the indicator function (ill-conditioning)
%  Geometry:        Structure that contains geometric- and meta-
%                   information about the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% Ke:            	Elemental stiffness matrix
%
%-------------------------------------------------------------------------%
function [Ke] = StiffMatFCM(xElem,yElem,h,D,QT,alpha,Geometry)

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

% Number of subcells for the numerical integration
nSC = numel(QT.domains);

% Local coordinates of the (corner) vertices of the subcells
xiQT  = QT.XI;
etaQT = QT.ETA;

% Loop over subcells
for j = 1:nSC
    % Mapping from subcell to reference domain
    % Centroid of the rectangular subcell
    xiC  = 0.5*(xiQT(2,j)+xiQT(1,j));
	etaC = 0.5*(etaQT(4,j)+etaQT(1,j));
    % Height and width of the subcell
	hxi = abs(xiQT(2,j)-xiQT(1,j));
	heta = abs(etaQT(4,j)-etaQT(1,j));
	% Jacobian matrix and determinant for the subcell
	Jsc    = 1/2*[hxi 0; 0 heta];
	detJsc = det(Jsc);
    
    % Loop over all integration points
    for i = 1:nIP
        % Compute the local coordinate of the integration point
        xi  = xiC  + 0.5*xiG(i)*hxi;
        eta = etaC + 0.5*etaG(i)*heta;
        % Jacobian matrix -- J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
        J = JacobianMatrix(xi,eta,xElem,yElem);
        % Determinant of the Jacobi matrix
        detJ = det(J);
        % Inverse of the Jacobian matrix
        Jinv = [J(2,2), -J(1,2); -J(2,1), J(1,1)]/det(J);
        % Matrix of the second order derivatives of the shape functions
        % w.r.t. the local coordinates
        Bl = [ SF_BFS_dxi2(xi,eta)'   ; ...
               SF_BFS_deta2(xi,eta)'  ; ...
              2*SF_BFS_dxideta(xi,eta)'];
        % Matrix of the second order derivatives of the shape functions
        % w.r.t. the global coordinates
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
        
        % Global coordinates of the integration points
        [xIP,yIP] = Local2GlobalCoordinatesQuad(xi,eta,xElem,yElem);
%         plot(xIP,yIP,'bx')
        % Is the integration point in the fictitious or physical domain
        [FCM,Idx,Hole] = PointInMembershipTest(xIP,yIP,Geometry);
        % Penalty-factor
        if FCM == 1 && Hole == 1
            % The integration point is located in one of the geometric
            % domains defined in the variable Geometry (FCM == 1) and the
            % domain is a void region (Hole == 1)
            val = alpha;
            C = D{Geometry.base_mat};
        elseif FCM == 1 && Hole == 0
            % The integration point is located in one of the geometric
            % domains defined in the variable Geometry (FCM == 1) and the
            % domain is an inclusion (Hole == 0)
            val = 1;
            C = D{Geometry.material{Idx}};
        elseif FCM == 1 && Hole == 2
            % The integration point is located in one of the geometric
            % domains defined in the variable Geometry (FCM == 1) and the
            % domain is the physical domain (Hole == 2)
            val = 1;
            C = D{Geometry.material{Idx}};
        elseif FCM == 0 && Hole == 0
            % The integration point is not located in one of the geometric
            % domains defined in the variable Geometry (FCM == 0) and
            % therefore, it is located in the base material
            val = 1;
            C = D{Geometry.base_mat};
        elseif FCM == 0 && Hole == 1
            % The integration point is not located in one of the geometric
            % domains defined in the variable Geometry (FCM == 0) and
            % therefore, it is located in the base material
            val = 1;
            C = D{Geometry.base_mat};
        elseif FCM == 0 && Hole == 2
            % The integration point is not located in one of the geometric
            % domains defined in the variable Geometry (FCM == 0) and
            % the geometry is the physical domain. Therefore, the point is
            % located outside of the computational domain.
            val = alpha;
            C = D{Geometry.base_mat};
        end
        
        % Stiffness matrix
        Ke = Ke + val*Bg'*C*Bg*h*detJ*detJsc*wG(i);
    end
end
    
end