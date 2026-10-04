%% Topology data
%-------------------------------------------------------------------------%
%
% (c) 04/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Derive the degree of freedom topology data from the nodal topology
% matrix.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:    Nodal topology (connectivity/coincidence) matrix
%
%-------------------------------------------------------------------------%
%% Output
%
% TopoD:    DOF topology (connectivity/coincidence) matrix
%
%-------------------------------------------------------------------------%
function [TopoD] = TopologyDOF(TopoN)

% Number of degrees of freedom per node
ndof = 4;

% DOF topology data (connectivity matrix)
Topo_aux1 = kron(TopoN,ones(1,ndof)*ndof);
Topo_aux2 = kron(ones(size(TopoN)),(ndof-1:-1:0));
TopoD = Topo_aux1 - Topo_aux2;

end