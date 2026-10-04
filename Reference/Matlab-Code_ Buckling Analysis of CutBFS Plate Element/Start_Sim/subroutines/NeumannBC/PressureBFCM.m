%% Pressure load (BFCM)
%-------------------------------------------------------------------------%
%
% (c) 06/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the pressure load.
%
%-------------------------------------------------------------------------%
%% Input
%
% (xE,yE):      Global coordinates of the nodes
%  QT:       	Quadtree data
%  mag:      	Magnitude of the pressure load
%  alpha:     	Value of the indicator function (ill-conditioning)
%  flag_func:  	Constant pressure load or functional dependency
%  Geometry:  	Structure that contains geometric- and meta-information 
% 				about the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% Fe:            	Elemental load vector due to pressure
%
%-------------------------------------------------------------------------%
function [Fe] = PressureBFCM(xE,yE,QT,mag,alpha,flag_func,Geometry)

% Set up the numerical integration
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

% Compute Fe
Fe = zeros(16,1);

% Number of subcells for the numerical integration
nSC = numel(QT.alphaB);

% Local coordinates of the (corner) vertices of the subcells
xiQT  = QT.XI;
etaQT = QT.ETA;
% Alpha values for the BF
alphaB = QT.alphaB;

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
        % Compute shape functions
        [N] = SF_BFS(xi,eta);
        % Jacobian matrix -- J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
        J = JacobianMatrix(xi,eta,xE,yE);
        % Determinant of the Jacobi matrix
        detJ = det(J);
        if flag_func == 1
            % Global coordinates of the integration point
            [xIP,yIP] = Local2GlobalCoordinatesQuad(xi,eta,xE,yE);
            % Amplitude of the pressure function
            amp = mag*AmplitudeFunc(xIP,yIP);
        else
            % Constant pressure
            amp = mag;
        end
        
        % Global coordinates of the integration points
        [xIP,yIP] = Local2GlobalCoordinatesQuad(xi,eta,xE,yE);
        % Is the integration point in the fictitious or physical domain
        [FCM,~,Hole] = PointInMembershipTest(xIP,yIP,Geometry);
        
        % Compute the gravity load vector for the whole subcell
        % (physical domain)
        if alphaB(j) == 11
            val = 1;
        % Compute the gravity load vector for the whole subcell
        % (fictitious domain)
        elseif alphaB(j) == -11
            val = -1 + alpha;
        % Compute the gravity load vector for the physical part of the
        % subcell
        elseif alphaB(j) == 1
            % Penalty-factor
            if FCM == 1
                % Inside a geometric domain defined in Geometry
                if Hole == 0
                	% Inclusion
                 	val = 1;
             	elseif Hole == 1
                  	% Void
                	continue;
               	elseif Hole == 2
                 	% Physical domain
                    val = 1;
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
                end 
            end
        % Compute the gravity load vector for the fictitious part of
        % the subcell
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
             	elseif Hole == 2
                  	% Physical domain
                   	continue;
                end
          	elseif FCM == 0
              	% Outside a geometric domain defined in Geometry
              	if flag_physical == 1
                        val = -1 + alpha;
                else
                 	continue;
             	end 
            end
        end

        % Gravity load
        Fe = Fe + val*N*amp*detJ*detJsc*wG(i);
    end
end
fmax = max(abs(Fe));
Fe(abs(Fe)<fmax*1e-10) = 0;
    
end