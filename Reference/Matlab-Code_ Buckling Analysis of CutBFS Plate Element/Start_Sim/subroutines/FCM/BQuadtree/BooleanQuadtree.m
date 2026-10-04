function [XI,ETA,alphaB,tree,t,varargout] = BooleanQuadtree(xCoord, ...
          yCoord, kMax, Geometry, globalCoord)
% Created by: Marton Petoe
% Date: 13.02.2019

%% Description
% The following algorithm performs a boolean-operations-based quadtree 
% decomposition of a given cell containing a discontinuity. The resulting
% boolean sub-cells will be then used in the context of the Boolean Finite
% Cell Method (B-FCM) for computing the integral of a discontinuous 
% function. The goal of this method is to use a significantly lower number
% of sub-cells, which consequently results in a lower number of integration
% points compared to the traditional FCM deploying a traditional quadtree-
% decomposition-based adaptive integration scheme.

% The algorithm is entirely based on the work "An extension of the finite 
% cell method using boolean operations" by A. Abedian and A. Düster.
% According to the authors, the algorithm has 3 main features:
%    1. If a sub-cell has the value alphaB = 1, it means that the 
%       integration of the PHYSICAL domain within the sub-cell is of
%       interest. If the sub-cell is cut, 4 children are created and
%       tested, and it is determined which domains (physical, fictitious,
% 		cut) they belong to. If (i) the ratio of the area of the PHYSICAL
% 		domain within the parent-cell to the area of the parent-cell is
% 		higher than a predefined tolerance, and (ii) at least 1 of the 4
% 		sub-cells is located entirely in the PHYSICAL domain, the parent
% 		cell is saved for integration, while all the the sub-cells get the
% 		values alphaB = -1.
%    2. If a sub-cell has the value alphaB = -1, it means that the 
%       integration of the FICTITIOUS domain within the sub-cell is of
%       interest. If the sub-cell is cut, 4 children are created and
%       tested, and it is determined which domains (physical, fictitious,
% 		cut) they belong to. If (i) the ratio of the area of the FICTITIOUS
%       domain within the parent-cell to the area of the parent-cell is
% 		higher than a predefined tolerance, and (ii) at least 1 of the 4
% 		sub-cells is located entirely in the FICTITIOUS domain, the parent
% 		cell is saved for integration, while all the the sub-cells get the
% 		values alphaB = 1.
%    3. The third feature should prevent the ill-conditioning during the 
%       numerical integration by distributing integration points in the 
%       entire cell and considering only the ones in the fictitious domain 
%       with a small penalty value alphaB << 1.

% ------------ Comments to the above mentioned three steps ----------------

% Comments to features #1 and #2:
%   If neither feature #1 nor #2 can be applied due to the violation of the
%   conditions (i) and/or (ii), the cut parent-cell splits into 4 sub-cells
%   having the same alphaB value as the parent-cell, then the parent cell
%   is deleted as in the traditional quadtree decomposition.
%
%   If the above mentioned conditions are fulfilled, the parent-cell will
%   not be deleted, but it will be saved with an alphaB value of 11 or -11.
%   These values not only indicate for the algorithm computing the integral
%   which domain is of interest, but they also indicate, that we wish to 
%   integrate the entire sub-cell while neglecting the discontinuity within
%   the sub-cell. For the integration of sub-cells with such alphaB
%   values no inside-outside tests have to be ran during the numerical
%   integration.

% Comments to feature #3:
%   The method applied in feature #3 further reduces the number of used
%   integration points. Additional description of this technique can be
% 	found in the following two papers by Abedian et al.:  
%   - "Finite cell-method compared to the h-version finite element method
%     for elasto-plastic problems"
%   - "The finite cell method for the J2 flow theory of plasticity"
%   Although the integration itself is not carried out by this algorithm,
%   we can already delete some sub-cells, whose contribution will be 
%   considered during the integration process by applying the method 
%   described in the feature #3 and in the above mentioned papers. In every
%   case, physical sub-cells with an alphaB value of -1 and fictitious sub-
%   cells with an alphaB value of 1 can be directly deleted. 
%
%   Note that this method is not equivalent to deleting every sub-cell
%   lying in the fictitious domain. Matter of fact, in certain cases the
%   algorithm results in sub-cells lying completely in the fictitious
%   domain, however, they are very much needed for an accurate integration.
%  
%   If we do not wish to apply feature #3, we may include the lines in the
%   code, that are currently commented out and marked with a "!" symbol.

%% Limitations
% According to the condition (i) in feature #1 and #2 (see previous
% section), the ratio between the area of the physical/fictitious domain to
% the area of the entire sub-cell should be evaluated. It is not 
% implemented yet at the current stage of this function. Therefore, 
% condition (i) is ignored and condition (ii) is modified so, that it 
% requires 2 sub-cells to be completely physical/fictitious in order for 
% features #1 and #2 to apply. 

%% Visualization
% For visualizing the results, use the function "displayBooleanQuadtree.m".

%% Input and output arguments
% INPUT -------------------------------------------------------------------
%   * xCoord, yCoord:  Global coordinates of the quadrilateral element
%   * kMax:            Refinement level of the quadtree decomposition
%   * points:          Global coordinates of the polygon defining the
%                      physical boundary
%
% OUTPUT ------------------------------------------------------------------
%   * XI, ETA:  [4 x nSC] matrices containing the local coordinates of the 
%               nSC sub-cells in column-wise order
%   * alphaB:   Integer vector containing the boolean parameters of the 
%               individual sub-cells. Following values are possible: 
%               1:  The integration of the physical domain is of interest
%               11: The integration of the physical domain is of interest
%                   while ignoring the discontinuity within the sub-cell
%              -1:  The integration of the fictitious domain is of interest
%              -11: The integration of the fictitious domain is of interest
%                   while ignoring the discontinuity within the sub-cell
%   * tree:     [kMax+1 x nSC] matrix, containing the locations of the nSC
%               sub-cells on the quadtree.
%   * t:        Time cost of the algorithm.


%% Initializing the output and internal variables
tic;

% At the very beginning, there is only one sub-cell, which in fact
% corresponds to the entire cell. Its local coordinates are:
XI(1:4,1) = [-1 1 1 -1];
ETA(1:4,1) = [-1 -1 1 1];

% According to the theory of the B-FCM (see description above), the
% parent-cell at the level k=0 always have the alphaB value 1, which means,
% that the integration of the physical domain is of interest. 
alphaB = 1;

% The refinement level kMax defines the depth of the quadtree. We need the
% kMax+1 digits for successfully labeling every branch and leaf of the tree.
% The cell at the refinement level k=0 gets the ID 1.
tree = zeros(kMax+1,1);
tree(1,1) = 1;

% The vector scFlags contains the flags of the sub-cells, indicating
% whether they are physical, fictitious or cut. One can assume, that the
% function "booleanQuadtree.m" will be called for a given cell, only if it
% is cut. Therefore, one can be certain, that the root-cell at k=0 has a 
% flag marking a cut cell. However, for the sake of safety, before we would
% perform the quadtree decomposition, we check the flag of the root-cell, 
% which corresponds to the first value in the vector scFlags. 
scFlags = cellFlag(XI(:,1), ETA(:,1), xCoord, yCoord, Geometry, kMax, 0);

% Since the function "cellFlag.m" uses a different encoding for the domains
% that we wish to use in this code, we have to transform those values.
% The function cellFlag.m uses the following encoding in regards of the
% domains: 0 = physical, 1 = fictitious, 2 = cut.
% For the sake of simplicity of this function, it is  beneficial to write
% the flags in the following form: 1 = physical, -1 = fictitious, 0 = cut.
if scFlags == 1
    scFlags = -1;
elseif scFlags == 0 
    scFlags = 1;
elseif scFlags == 2
    scFlags = 0;
end

% !  booleanBranch = 0;

% The constant "threshold" indicates the number of sub-cells that have to be 
% either entirely physical or fictitious among the 4 children of a given 
% parent-cell in order for features #1 and #2 to apply. For further info
% please see the sections 'Description' and 'Limitations'.
treshold = 2;

% The matrices XI and ETA containing the local coordinates of the sub-cells
% are dynamically changing arrays. At every refinement level k, the found
% new sub-cells will be added to the arrays, while the parent-cells will be 
% deleted. The values idxStart and idxEnd define a vector of integers, 
% referring to the IDs of the newly found sub-cells at the refinement 
% level k-1 that need to be examined at the refinement level k. Since at
% the starting level k=0 we only have a single "sub-cell", both values are 1.
idxStart = 1;
idxEnd = 1;

%% Boolean quadtree decomposition

% Proceed the quadtree decomposition until the maximal three depth kMax is
% reached.
for k = 1:kMax
    
    scFound = 0; % Found sub-cells at the current refinement level k
    scDelete = []; % Sub-cells to delete
    
    % Iterate through the sub-cells, created on the previous refinement
    % level k-1.
    for i = idxStart:idxEnd
        
        % A given sub-cells can only undergo the following procedures, when
        % it is cut, i.e. its flag is equal to 0.
        if scFlags(i) == 0
            
            % Create 4 sub-cells of the current cut parent-cell
            [scX,scY] = DecompositionQuadrilateral(XI(:,i),ETA(:,i));
            
            % Create IDs for the 4 sub-cells, that mark their locations on
            % the quadtree.
            branch = [tree(:,i), tree(:,i), tree(:,i), tree(:,i)];
            branch(k+1,:) = [1 2 3 4];
            
            % Iterate through the 4 sub-cells and check which domains do 
            % they belong to.
            flags = [0, 0, 0, 0]; % Initial array
            for j = 1:4
                flags(j) = cellFlag(scX(:,j), scY(:,j), xCoord, yCoord,...
                           Geometry, kMax, k);
 
                % The function cellFlag.m uses the following encoding in 
                % regards of the domains:
                % 0 = physical, 1 = fictitious, 2 = cut.
                % For the sake of simplicity of this function, it is
                % beneficial to write the flags in the following form:
                % 1 = physical, -1 = fictitious, 0 = cut.
                if flags(j) == 1
                    flags(j) = -1;
                elseif flags(j) == 0
                    flags(j) = 1;
                elseif flags(j) == 2
                    flags(j) = 0;
                end
            end
            
            % If one of the following two criteria is fulfilled, the given
            % parent-cell will be saved for integration and the alphaB 
            % value of the sub-cells changes its sign.
            % The alphaB values 11 and -11 mean, that the whole sub-cell is
            % integrated as either entirely physical or fictitious. Unlike
            % in the case of the alphaB values 1 and -1, the discontinuity
            % in the sub-cell will be ignored during the integration.
            
            if alphaB(i) == 1 && length(flags(flags == 1)) >= treshold
                alphaB(i) = 11;
                scAlphaB = -1;
            elseif alphaB(i) == -1 && length(flags(flags == -1)) >= treshold
                alphaB(i) = -11;
                scAlphaB = 1;
            else
                % If none of the above conditions are fulfilled, the main
                % advantages of the B-FCM cannot be applied in this 
                % iteration. In this case, the sub-cells inherit the alphaB 
                % value of the parent-cell and then parent-cell gets
                % deleted. 
                scAlphaB = alphaB(i);
                % The deleting of the indicated sub-cells will be carried
                % out later
                scDelete(end+1) = i; 
            end
            
% !          if booleanBranch(i) == 1 || alphaB(i) == 11 || alphaB(i) == -11
                scX(:,flags == -scAlphaB) = [];
                scY(:,flags == -scAlphaB) = [];
                branch(:,flags == -scAlphaB) = [];
                flags(flags == -scAlphaB) = [];
                
                % Number of sub-cells we keep in the current iteration
                nsc = size(scX,2); % = size(scY,2)
% !               b = 1; % Boolean value for the vector booleanBranch
% !           else
% !               nsc = 4;
% !               b = 0; % Boolean value for the vector booleanBranch
% !           end
                  
            % Add sub-cells and their meta-data to the corresponding arrays
            XI(:,end+1:end+nsc) = scX;
            ETA(:,end+1:end+nsc) = scY;
            alphaB(end+1:end+nsc) = scAlphaB * ones(1,nsc);
            tree(:,end+1:end+nsc) = branch;
% !           booleanBranch(end+1:end+nsc) = b * ones(1,nsc);        
            scFlags(end+1:end+nsc) = flags;
            
            % Current number of found sub-cells at the refinement level k
            scFound = scFound + nsc;
		else
           % If a sub-cell i has a scFlag(i) value 1/-1, it means  that the
           % integration of the physical/fictitious domain is of interest 
           % with a scaling factor of 1/-1. However, if the sub-cells is
           % not cut, we can assign alphaB = 11/-11 to the sub-cell and so
           % we can spare some inside-outside tests for the given sub-cell
           % during the numerical integration.
           alphaB(i) =  sign(scFlags(i)) * 11;
        end
    end
    
    % Delete the sub-cells and their meta-data indicated by the vector
    % scDelete.
    XI(:,scDelete) = [];
    ETA(:,scDelete) = [];
    alphaB(scDelete) = [];
    tree(:,scDelete) = [];
    scFlags(scDelete) = [];
% !   booleanBranch(scDelete) = [];
    
    % Compute iteration interval for the next loop
    idxEnd = size(XI,2);
    idxStart = idxEnd-scFound+1;
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