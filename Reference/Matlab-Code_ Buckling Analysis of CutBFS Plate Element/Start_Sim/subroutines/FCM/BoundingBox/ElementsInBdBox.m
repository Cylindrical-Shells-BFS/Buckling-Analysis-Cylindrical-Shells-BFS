%% Determine all elements in the bounding box of the fictitious domain
%-------------------------------------------------------------------------%
%
% (c) 10.09.2017, Sascha Duczek, OvGU Magdeburg, Germany
%
% The subroutine determines all elements whose nodes (at least one) are 
% located withing the bounding box of a void or inclusion. This reduces the
% numerical costs to determine cut and fictitious elements.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:		Nodal topology (connectivity/coincidence) matrix
% xCoor: 		Global x-coordinates of all nodes in the FE model
% yCoor:  		Global y-coordinates of all nodes in the FE model
% BdBox: 		Bounding box of the fictitious domain
%
%-------------------------------------------------------------------------%
%% Output
%
% ElemBdBox:    Elements within the bounding box of the fictitious
%               domain/inclusion
%
%-------------------------------------------------------------------------%
function [ElemBdBox] = ElementsInBdBox(TopoN,xCoor,yCoor,BdBox)

% Number of bounding boxes
nBB = numel(BdBox.coor);

% Initialize the result variable
ElemBdBox = cell(nBB,1);

for i = 1:nBB
	% Determine all nodes whose coordinates are less or equal to xmax, ymax
    % and greater or equal to xmin, ymin. All elements that are within the
    % bounding box are tested. This is a significant reduction of the
    % computational effort for voids and inclusions. If the geometric
    % domain describes a physical boundary (outer boundary of the
    % computational domain) the costs are not much reduced. For the
    % physical boundary one would need to derive an interior bounding box
    % which is not straight forward for arbitrary shapes.
 	idxXmax = find(xCoor <= BdBox.coor{i}(2,1));
   	idxXmin = find(xCoor >= BdBox.coor{i}(1,1));
   	idxYmax = find(yCoor <= BdBox.coor{i}(2,2));
  	idxYmin = find(yCoor >= BdBox.coor{i}(1,2));
        
	% Intersect all nodal sets
 	idxX  = intersect(idxXmax,idxXmin);
  	idxY  = intersect(idxYmax,idxYmin);
 	idxXY = intersect(idxX,idxY);
    
    % Determine the elements that are associated a node
    [Patch] = ConnectivityNodeElementBB(TopoN,idxXY);

    % Determine all elements that share the identified nodes
    ElemBdBox{i} = unique(Patch);
    ElemBdBox{i} = ElemBdBox{i}(ElemBdBox{i} ~= 0);
end

end