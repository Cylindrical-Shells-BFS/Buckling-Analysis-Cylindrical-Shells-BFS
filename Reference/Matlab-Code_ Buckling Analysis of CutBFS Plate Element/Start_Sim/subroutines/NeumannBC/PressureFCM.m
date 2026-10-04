%% Pressure load (FCM)
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
function [Fe] = PressureFCM(xE,yE,QT,mag,alpha,flag_func,Geometry)

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
        % Compute shape functions
        [N] = SF_BFS(xi,eta);
        % Jacobian matrix -- J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
        J = JacobianMatrix(xi,eta,xE,yE);
        % Determinant of the Jacobi matrix
        detJ = det(J);
        % Global coordinates of the integration point
    	[xIP,yIP] = Local2GlobalCoordinatesQuad(xi,eta,xE,yE);
        if flag_func == 1
            % Amplitude of the pressure function
            amp = mag*AmplitudeFunc(xIP,yIP);
        else
            % Constant pressure
            amp = mag;
        end
        
        % Is the integration point in the fictitious or physical domain
        [FCM,~,Hole] = PointInMembershipTest(xIP,yIP,Geometry);
        % Penalty-factor
        if FCM == 1 && Hole == 1
            % The integration point is located in one of the geometric
            % domains defined in the variable Geometry (FCM == 1) and the
            % domain is a void region (Hole == 1)
            val = alpha;
        elseif FCM == 1 && Hole == 0
            % The integration point is located in one of the geometric
            % domains defined in the variable Geometry (FCM == 1) and the
            % domain is an inclusion (Hole == 0)
            val = 1;
        elseif FCM == 1 && Hole == 2
            % The integration point is located in one of the geometric
            % domains defined in the variable Geometry (FCM == 1) and the
            % domain is the physical domain (Hole == 2)
            val = 1;
        elseif FCM == 0 && Hole == 0
            % The integration point is not located in one of the geometric
            % domains defined in the variable Geometry (FCM == 0) and
            % therefore, it is located in the base material
            val = 1;
        elseif FCM == 0 && Hole == 1
            % The integration point is not located in one of the geometric
            % domains defined in the variable Geometry (FCM == 0) and
            % therefore, it is located in the base material
            val = 1;
        elseif FCM == 0 && Hole == 2
            % The integration point is not located in one of the geometric
            % domains defined in the variable Geometry (FCM == 0) and
            % the geometry is the physical domain. Therefore, the point is
            % located outside of the computational domain.
            val = alpha;
        end

        % Gravity load
        Fe = Fe + val*N*amp*detJ*detJsc*wG(i);
    end
end
fmax = max(abs(Fe));
Fe(abs(Fe)<fmax*1e-10) = 0;
    
end