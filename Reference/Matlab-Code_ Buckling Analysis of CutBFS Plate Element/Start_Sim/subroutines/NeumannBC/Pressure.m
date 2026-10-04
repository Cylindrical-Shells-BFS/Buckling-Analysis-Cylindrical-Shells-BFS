%% Pressure load
%-------------------------------------------------------------------------%
%
% (c) 06/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the pressure load.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:        Nodal topology (connectivity/coincidence) matrix
% TopoD:        DOF topology (connectivity/coincidence) matrix
% (x,y):        Global coordinates of the nodes
% Elset:        Element set for the pressure load
%  mag:      	Magnitude of the pressure load
%  flag_func:  	Constant pressure load or functional dependency
%
%-------------------------------------------------------------------------%
%% Output
%
% Fp:            	Elemental load vector due to pressure
%
%-------------------------------------------------------------------------%
function [Fp] = Pressure(TopoN,TopoD,x,y,Elset,mag,flag_func)

% Number of elements with applied pressure load
nE = numel(Elset);

% Number of degrees of freedom
nDOF = max(max(TopoD));

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

% Compute Fp
Fp = zeros(nDOF,1);

% Loop over all elements with an applied pressure load
for j = 1:nE
    % Nodal coordinates of the vertices
    xElem = x(TopoN(Elset(j),:));
    yElem = y(TopoN(Elset(j),:));

    % Elemental force vector
    Fe = zeros(16,1);

    % Loop over all integration points
    for i = 1:nIP
        % Compute shape functions
        [N] = SF_BFS(xiG(i),etaG(i));
        % Jacobian matrix -- J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
        J = JacobianMatrix(xiG(i),etaG(i),xElem,yElem);
        % Determinant of the Jacobi matrix
        detJ = det(J);
        if flag_func == 1
            % Global coordinates of the integration point
            [xIP,yIP] = Local2GlobalCoordinatesQuad(xiG(i),etaG(i),...
                        xElem,yElem);
            % Amplitude of the pressure function
            amp = mag*AmplitudeFunc(xIP,yIP);
        else
            % Constant pressure
            amp = mag;
        end
        
        % Gravity load
        Fe = Fe + N*amp*detJ*wG(i);
    end
    fmax = max(abs(Fe));
    Fe(abs(Fe)<fmax*1e-10) = 0;
    
    % Assemble global force vector
    Fp(TopoD(Elset(j),:)) = Fp(TopoD(Elset(j),:)) + Fe;

end

end