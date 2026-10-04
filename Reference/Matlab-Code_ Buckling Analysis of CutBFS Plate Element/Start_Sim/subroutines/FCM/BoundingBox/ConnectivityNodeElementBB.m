%% For each node all elements are listed
%-------------------------------------------------------------------------%
%
% (c) 12.09.2017, Sascha Duczek, OvGU Magdeburg, Germany
%
% Determines all elements that share a node. The results is matrix where
% the row denotes the node number and the entries in that row correspond to
% the finite elements that share this node.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:	Nodal topology (connectivity/coincidence) matrix
% Nodes:    Nodes that are checked
%
%-------------------------------------------------------------------------%
%% Output
%
% Patch: 	Matrix of element numbers
%
%-------------------------------------------------------------------------%
function [Patch] = ConnectivityNodeElementBB(TopoN,Nodes)

% Determine the number of nodes that are checked
if isempty(Nodes)
    % All nodes in the whole FE-model are checked
    NumNodes = max(max(TopoN));
else
    % Only specific nodes are checked
    NumNodes = numel(Nodes);
end

% Array containing the element numbers sharing a specific node
Patch = zeros(NumNodes,8);

% Loop over all nodes
for i = 1:NumNodes
	[row,~] = find(TopoN == Nodes(i));
	Patch(i,1:size(row,1)) = row;
end

end