%% Determine cut and fictitious elements
%-------------------------------------------------------------------------%
%
% (c) 01.09.2016, Sascha Duczek, OvGU Magdeburg, Germany
%
% Determines all elements that are cut by the boundary of the void or
% inclusion. Elements that are entirely located in a void region are
% fictitious elements and can be discarded in the analysis.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:        Nodal topology (connectivity/coincidence) matrix
% xCoor:        x-coordinates of the nodes
% yCoor:        y-coordinates of the nodes
% ElemBdBox: 	Elements within the bounding box of the fictitious
%               domain/inclusion
% Geometry:		Structure that contains geometric- and meta-information
% 				about the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% Cut_Elem: 	Cell array containing a vector for each void region/
%               inclusion listing the numbers of all finite elements that
%               are intersected by the physical boundary of the domain.
% Fict_Elem: 	Cell array containing a vector for each void region/
%               inclusion listing the numbers of all finite elements that
%               are located in the fictitious domain.
%
%-------------------------------------------------------------------------%
function [CutElem,FictElem] = Cut_Elements(TopoN,xCoor,yCoor,ElemBdBox,...
                              Geometry)

% Number of bounding boxes
nBB = numel(ElemBdBox);

% Initialize result variables
CutElem  = cell(nBB,1);
FictElem = cell(nBB,1);
HoleElem = cell(nBB,1);

% We need to check all elements that are within the bounding box of the
% individual void regions or inclusions
for k = 1:nBB                            
    % Number of finite elements in the ith bounding box
    num_Elem = numel(ElemBdBox{k});
    % Vector containing flags whether the ith element is cut or physical
    elem_vec_cut = zeros(num_Elem,1);
    % Vector containing flags whether the ith element is fictitious
    elem_vec_fict = zeros(num_Elem,1);
    % Geometry data for current domain
    geo.type{1}   = Geometry.type{k};
    geo.domain{1} = Geometry.domain{k};
    geo.shape{1}  = Geometry.shape{k};
    geo.param{1}  = Geometry.param{k};
    geo.points{1} = Geometry.points{k};

    counter_cut  = 0;
    counter_fict = 0;
    for i = 1:num_Elem
        % Nodal coordinates
        x = xCoor(TopoN(ElemBdBox{k}(i),:));
        y = yCoor(TopoN(ElemBdBox{k}(i),:));

        % Number of nodes for the ith finite element
        nNodes = numel(x);
        
        % Which points are located in the void region/inclusion?
		[FCM_val_vec,~,Hole] = PointInMembershipTest(x,y,geo);
        FCM_val  = sum(FCM_val_vec);
        Hole_val = sum(Hole)/numel(Hole);

        % Cut elements are all elements that belong to two different
		% domains, i.e. either physical domain and void region or
		% physical domain and inclusion
        if FCM_val >= 1 && FCM_val <= nNodes - 1
            counter_cut = counter_cut + 1;
            elem_vec_cut(counter_cut) = ElemBdBox{k}(i);
        % Fictitious element (all nodes are within the geometric domain and
        % the domain is a void region); the nearest domain indicated by
        % IdxDomain should also be the domain the is investigated in this
        % step
        elseif FCM_val == nNodes && Hole_val == 1
            counter_fict = counter_fict + 1;
            elem_vec_fict(counter_fict) = ElemBdBox{k}(i);
     	% Fictitious element (all nodes are outside of the geometric domain
        % and the domain is a physical region); the nearest domain 
        % indicated by IdxDomain should also be the domain the is 
        % investigated in this step
        elseif FCM_val == 0 && Hole_val == 2
            counter_fict = counter_fict + 1;
            elem_vec_fict(counter_fict) = ElemBdBox{k}(i);
        end
    end

    % Delete all zero components in the CutElem vector
    idx_cut  = find(elem_vec_cut == 0,1,'first');
    idx_fict = find(elem_vec_fict == 0,1,'first');
    if isempty(idx_cut)
        CutElem{k}  = elem_vec_cut;
    else
        CutElem{k}  = elem_vec_cut(1:idx_cut-1);
    end
    % Delete all zero components in the FictElem vector
    if isempty(idx_fict)
        FictElem{k} = elem_vec_fict;
    else
        FictElem{k} = elem_vec_fict(1:idx_fict-1);
    end
end

end