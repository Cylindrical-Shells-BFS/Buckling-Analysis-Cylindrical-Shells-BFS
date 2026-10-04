%% Mapping from the local to the global domain - Quadrilateral
%-------------------------------------------------------------------------%
%
% (c) 25.08.2016, Sascha Duczek, OvGU Magdeburg, Germany
%
% Compute the global coordinates in a quadrilateral finite element.
%
%-------------------------------------------------------------------------%
%% Input
%
% (xi,eta):	Local coordinates of the points
% (xE,yE):	Global coordinates of the nodes of the quadrilateral element
%
%-------------------------------------------------------------------------%
%% Output
%
% (xg,yg):  Global coordinates of the points
%
%-------------------------------------------------------------------------%
function [xg,yg] = Local2GlobalCoordinatesQuad(xi,eta,x,y)

% Number of points
n = numel(xi);

% Shape functions of a 4-node element
[N] = SF_4nodeElem(xi,eta);

% Multiply each column of a matrix by the corresponding element of a vector
X = repmat(x,1,n);
Y = repmat(y,1,n);

% Matrix to save the coordinates of the background nodes
xg = sum(N.*X)'; 
yg = sum(N.*Y)';

end