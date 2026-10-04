%% Mesh generator
%-------------------------------------------------------------------------%
%
% (c) 03/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Generate a simple structured mesh.
%
%-------------------------------------------------------------------------%
%% Input
%
% ndof:     Number of degrees of freedom per nodes
% nx:       Number of element in x-direction
% ny:       Number of element in y-direction
% lx:       Dimension x-direction
% ly:       Dimension y-direction
%
%-------------------------------------------------------------------------%
%% Output
%
% TopoN:    Nodal topology (connectivity/coincidence) matrix
% TopoD:    DOF topology (connectivity/coincidence) matrix
% nNodes:   Number of nodes
% nDOF:     Number of degrees of freedom (DOF)
% xCoor:    Global x-coordinates of the nodes
% yCoor:    Global y-coordinates of the nodes
%
%-------------------------------------------------------------------------%
function [TopoN,TopoD,nNodes,nDOF,xCoor,yCoor] = MeshGenerator(ndof,nx,...
                                                 ny,lx,ly)

% Number of elements
nElem = nx*ny;
% Nodal topology data (connectivity matrix)
TopoN = zeros(nElem,4);
counter1 = 0;
counter2 = 0;
for i = 1:nElem
    TopoN(i,:) = [1, 2, nx+3, nx+2] + counter1 + counter2;
    counter2 = counter2 + 1;
    if mod(i,nx) == 0
        counter1 = counter1 + nx + 1;
        counter2 = 0;
    end
end
% Number of nodes
nNodes = max(max(TopoN));
% Number of degrees of freedom
nDOF = nNodes*ndof;
% Nodal coordinates
xCoor = zeros(nNodes,1);
yCoor = zeros(nNodes,1);
for i = 1:ny+1
    xCoor(1+(nx+1)*(i-1):(nx+1)*i) = linspace(0,lx,nx+1);
    yCoor(1+(nx+1)*(i-1):(nx+1)*i) = ly/ny*(i-1);
end
% DOF topology data (connectivity matrix)
Topo_aux1 = kron(TopoN,ones(1,ndof)*ndof);
Topo_aux2 = kron(ones(size(TopoN)),(ndof-1:-1:0));
TopoD = Topo_aux1 - Topo_aux2;

end