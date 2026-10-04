function [XI,ETA,domains,tree,t,varargout] = modQuadtree(xCoord, ...
          yCoord, kMax, Geometry, checkLeafCells, globalCoord)
		  
% Created by: Marton Petoe
% Date: 12.02.2019

%% Description
% The following code creates a kMax-deep quadtree decomposition of a given 
% cell. The code originally was created by Christoph Neuhaus (07.05.2013),
% Sacha Duczek (23.05.2014) and then modified and further extended by 
% Marton Petoe (12.02.2019).

% The checking, which domains do the sub-cells belong to, is no longer done
% within this code, but by the external function cellFlag.m.

%% Visualization
% For visualizing the results, use the function "displayQuadtree.m".

%% Input and output arguments
% INPUT -------------------------------------------------------------------
%   * xCoord, yCoord:  Global coordinates of the quadrilateral element
%   * kMax:            Refinement level of the quadtree decomposition
%   * Geometry:        Structure that contains geometric- and meta-
%                      information about the given boundary.
%   * checkLeafCells:  Optional boolean parameter. If its not given, or it 
%                      is given with the value "false", then flags of the
%                      leaf-cells remain unchecked. If it is given and its 
%                      value is "true", the leaf-cells will be also 
%                      checked. For further information please see the 
%                      section "Checking the leaf-cells" at the end of this
%                      function.
%
% OUTPUT ------------------------------------------------------------------
%   * XI, ETA:  [4 x nSC] matrices containing the local coordinates of the 
%               nSC sub-cells in column-wise order
%   * domains:  Integer vector containing the domain labels (flags) of the
%               nSC sub-cells. Following encoding is used:
%               0 = physical, 1 = fictitious, 2 = cut, 3 = not checked
%   * tree:     [kMax+1 x nSC] matrix, containing the locations of the nSC
%               sub-cells on the quadtree.
%   * t:        Time cost of the algorithm.

%% Initializing the output and internal variables 
tic

% At the very beginning, there is only one sub-cell, which in fact
% corresponds to the entire cell. Its local coordinates are:
XI(1:4,1) = [-1 1 1 -1];
ETA(1:4,1) = [-1 -1 1 1];

% Although it is likely, that the algorithm is only called, once it is
% certain, that a given cell is cut, for safety reasons we declare the cell
% at the beginning of the quadtree decomposition to be "not checked".
domains = 3; % 0 = physical, 1 = fictitious, 2 = cut, 3 = not checked

% The refinement level kMax defines the depth of the quadtree. We need the
% kMax+1 digits for successfully labeling every branch and leaf of the 
% tree. The cell at the refinement level k=0 gets the ID 1.
tree = zeros(kMax+1,1);
tree(1,1) = 1;

% The matrices XI and ETA containing the local coordinates of the sub-cells
% are dynamically changing arrays. At every refinement level k, the found
% new sub-cells will be added to the arrays, while the parent-cells will be 
% deleted. The values idxStart and idxEnd define a vector of integers, 
% referring to the IDs of the newly found sub-cells at the refinement 
% level k-1 that need to be examined at the refinement level k. Since at
% the starting level k=0 we only have a single sub-cell, both values are 1.
idxStart = 1;
idxEnd = 1;

%% Quadtree Decomposition

% Proceed the quadtree decomposition until the maximal three depth kMax is
% reached.
for k = 1:kMax

    scFound = 0; % Found sub-cells at the current refinement level k
    scDelete = []; % Sub-cells to delete
    
    for i = idxStart:idxEnd
        % Iterate through the sub-cells. For the first step there is no sub-
        % cell created, hence the sub-cell is the entire cell.
        
        flag = cellFlag(XI(:,i), ETA(:,i), xCoord, yCoord, Geometry, ...
		       kMax, k-1);
        
        if flag == 2
            % The flag value 2 stands for a cut sub-cell.
            
            [scX,scY] = DecompositionQuadrilateral(XI(:,i),ETA(:,i));
            
            % Create IDs for the 4 sub-cells, that mark their locations on
            % the quadtree.
            branch = [tree(:,i), tree(:,i), tree(:,i), tree(:,i)];
            branch(k+1,:) = [1 2 3 4];
            
            % Add the coordinates of the 4 sub-cells
            XI(:,end+1:end+4) = scX;
            ETA(:,end+1:end+4) = scY;
            % The new sub-cells are unchecked    
            domains(end+1:end+4) = [3,3,3,3];        
            tree(:,end+1:end+4) = branch;
            
            % The deleting of the indicated sub-cells will be carried out
            % later
            scDelete(end+1) = i;
                       
            % Current number of found sub-cells at the refinement level k
            scFound = scFound + 4;                
        else
            % If a given sub-cell is not cut, we can clearly assign a
            % domain to it according to its flag.
            domains(i) = flag;     
        end
    end
    
    % Delete the sub-cells and their meta-data indicated by the vector
    % scDelete.
    XI(:,scDelete) = [];
    ETA(:,scDelete) = [];
    domains(scDelete) = [];
    tree(:,scDelete) = [];
    
    % Compute iteration interval for the next loop
    idxEnd = size(XI,2);
    idxStart = idxEnd-scFound+1;
end

%% Checking the leaf-cells
% At this stage, every sub-cell has a clear domain flag, except the
% sub-cells on the last level of refinement, also known as leaf-cells. This
% is due to the mechanics of the above code: The flags of the sub-cells
% created on the refinement level k will be defined on the refinement level
% k+1. 

% This is not necessary a problem, if the resulting sub-cells will be 
% directly used for the numerical integration in the context of the FCM.
% The integration points will be distributed in the leaf-cells and tested,
% whether they belong to the physical or to the fictitious domain.

% However, in some cases the numerical integration is not carried out 
% directly on the sub-cells resulting from the quadtree decomposition, but
% on a set of rectangular domains with a lower cardinality. These 
% rectangular domains are obtained by covering all the physical sub-cells 
% resulting from the quadtree decomposition. Therefore it is of our 
% interest to find ALL the physical sub-cells and so an examination of the
% the leaf-cells is necessary.   

% The variables idxStart and idxEnd were already defined at the end of the
% last refinement level in the previous section.
if exist('checkLeafCells','var') && checkLeafCells == true
    for i = idxStart:idxEnd
        domains(i) = cellFlag(XI(:,i), ETA(:,i), xCoord, yCoord, ...
		             Geometry, kMax, k);
    end
end

t = toc;

if exist('globalCoord','var')
	if globalCoord == true
	   XX = zeros(4,size(XI,2));
	   YY = XX;
	   for i = 1:size(XI,2)
		   [xg,yg] = Local2GlobalCoordinatesQuad(XI(:,i),ETA(:,i),...
                     xCoord,yCoord);
		   XX(:,i) = xg;
		   YY(:,i) = yg;
	   end
	   
	   varargout{1} = XX;
	   varargout{2} = YY;
	end
end

end