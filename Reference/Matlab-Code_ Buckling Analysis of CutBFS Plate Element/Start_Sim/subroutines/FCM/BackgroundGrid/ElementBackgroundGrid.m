%% Element background grid
%-------------------------------------------------------------------------%
%
% (c) 01.09.2016, Sascha Duczek, OvGU Magdeburg, Germany
%
% Compute the global coordinates of a fine element background grid.
%
%-------------------------------------------------------------------------%
%% Input
%
% [x,y]:    Global coordinates of the nodes of a finite element
% n:        Number of division for the elemental background grid
%
%-------------------------------------------------------------------------%
%% Output
%
% Coor:  	Global coordinates of the elemental background grid
%
%-------------------------------------------------------------------------%
function [Coor] = ElementBackgroundGrid(x,y,n)

% Grid points
[xi,eta] = BackgroundGridQuadrilateral(n);

% Shape functions of a 4-node element
% Interior points
[NI] = SF_4nodeElem(xi.Interior,eta.Interior);
% Boundary points
[NB] = SF_4nodeElem(xi.Boundary,eta.Boundary);

% Multiply each column of a matrix by the corresponding element of a vector
% Global coordinates of the interior points
XI = repmat(x,1,(n-2)*(n-2));
YI = repmat(y,1,(n-2)*(n-2));
% Global coordinates of the boundary points
XB = repmat(x,1,4*(n-1));
YB = repmat(y,1,4*(n-1));

% Matrix to save the coordinates of the background nodes
Coor.Interior = [sum(NI.*XI), ; sum(NI.*YI)]';
Coor.Boundary = [sum(NB.*XB), ; sum(NB.*YB)]';

end