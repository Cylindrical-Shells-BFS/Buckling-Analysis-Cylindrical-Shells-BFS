%% Stiffness matrix (BFCM): Bogner-Fox-Schmit plate element
%-------------------------------------------------------------------------%
%
% (c) 05/05/2019, Sascha Duczek, UNSW Sydney, Australia
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
%  flagBFCM:        Quadtree decomposition method
%  Geometry:        Structure that contains geometric- and meta-
%                   information about the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% Ke:            	Elemental stiffness matrix
%
%-------------------------------------------------------------------------%
function [Ke] = StiffMatBFCM(xElem,yElem,h,D,QT,alpha,flagBFCM,Geometry)

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
nSC = numel(QT.alphaB);

% Local coordinates of the (corner) vertices of the subcells
xiQT  = QT.XI;
etaQT = QT.ETA;

% Alpha values for the BFCM
if flagBFCM == 1
    alphaB = QT.alphaB;
end

% Does the variable Geometry define the physical domain?
nGeo = numel(Geometry.type);
flag_physical = 0;
for i = 1:nGeo
    if strcmp(Geometry.domain{i},'Physical')
        flag_physical = 1;
        break
    end
end

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
     	% Is the integration point in the fictitious or physical domain
        [FCM,Idx,Hole] = PointInMembershipTest(xIP,yIP,Geometry);
            
        if flagBFCM == 1
            % Compute the stiffness matrix for the whole subcell (physical
            % domain)
            if alphaB(j) == 11
                val = 1;
                % Material properties
                C = D{Geometry.material{Idx}};
            % Compute the stiffness matrix for the whole subcell
            %(fictitious domain)
            elseif alphaB(j) == -11
                val = -1 + alpha;
                % Material properties
                C = D{Geometry.material{Idx}};
            % Compute the stiffness matrix for the physical part of the
            % subcell
            elseif alphaB(j) == 1
                % Penalty-factor
                if FCM == 1
                    % Inside a geometric domain defined in Geometry
                    if Hole == 0
                        % Inclusion
                        val = 1;
                        % Material properties
                        C = D{Geometry.material{Idx}};
                    elseif Hole == 1
                        % Void
                        continue;
                    elseif Hole == 2
                        % Physical domain
                        val = 1;
                        % Material properties
                        C = D{Geometry.material{Idx}};
                    end
                elseif FCM == 0
                    % Outside of all geometric domains defined in Geometry
                    if flag_physical == 1
						% If a physical domain is part of the geometry
						% description, FCM = 0 means that we are in the
						% fictitious domain.
                        continue;
                    else
						% If no physical domain is part of the geometry
						% description, FCM = 0 means that we are in the
						% domain with the base material.
                        val = 1;
                        % Material properties
                        C = D{Geometry.material{Geometry.base_mat}};
                    end 
                end
            % Compute the stiffness matrix for the fictitious part of the
            % subcell
            elseif alphaB(j) == -1
                % Penalty-factor
                if FCM == 1
                    % Inside a geometric domain defined in Geometry
                    if Hole == 0
                        % Inclusion
                        continue;
                    elseif Hole == 1
                        % Void
                        val = -1 + alpha;
                        % Material properties
                        C = D{Geometry.material{Geometry.base_mat}};
                    elseif Hole == 2
                        % Physical domain
                        continue;
                    end
                elseif FCM == 0
                    % Outside a geometric domain defined in Geometry
                    if flag_physical == 1
                        val = -1 + alpha;
                        % Material properties
                        C = D{Geometry.material{Geometry.base_mat}};
                    else
                        continue;
                    end 
                end
            end
        else
            % Penalty-factor
            if FCM == 1
                % Inside a geometric domain defined in Geometry
                if Hole == 0
                    % Inclusion
                    val = 1;
                    % Material properties
                    C = D{Geometry.material{Idx}};
                elseif Hole == 1
                    % Void
                    val = alpha;
                    % Material properties
                    C = D{Geometry.material{Geometry.base_mat}};
                elseif Hole == 2
                    % Physical domain
                    val = 1;
                    % Material properties
                    C = D{Geometry.material{Idx}};
                end
            elseif FCM == 0
                % Outside a geometric domain defined in Geometry
                if flag_physical == 1
                    val = alpha;
                    % Material properties
                    C = D{Geometry.material{Geometry.base_mat}};
                else
                    val = 1;
                    % Material properties
                    C = D{Geometry.material{Geometry.base_mat}};
                end 
            end
        end
        
        % Stiffness matrix
        Ke = Ke + val*Bg'*C*Bg*h*detJ*detJsc*wG(i);
    end
end
    
end