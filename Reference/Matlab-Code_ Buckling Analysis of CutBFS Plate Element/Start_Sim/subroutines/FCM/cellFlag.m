function [ flag ] = cellFlag( xiCoord, etaCoord, xCoord, yCoord, ...
					Geometry, kMax, k)
% Created by: Marton Petoe
% Date: 11.02.2019

% Original authors: 
% -Sascha Duczek (23.05.2014) 
% -Christoph Neuhaus (07.05.2013)

%% Description 
% The following code takes a sub-cell at a given refinement level (depth) k
% of a quadtree decomposition and checks, whether it is cut by the boundary
% or not. 

% The code originally was created by Christoph Neuhaus (section "Test 1") 
% and extended by Sascha Duczek (section "Test 2"), as mentioned above, and
% then it was slightly modified and placed into the current separate code 
% by Marton Petoe.

%% Input and output arguments
% INPUT: ------------------------------------------------------------------
% xiCoord, etaCoord: Local coordinates of the corners of the given sub-cell
% xCoord, yCoord:    Global coordinates of the corners of the given cell
% Geometry:          Structure that contains geometric- and meta-
%                    information about the given boundary.
% kMax:              Maximal refinement level of the quadtree decomposition
% k:                 Current refinement level (0 <= k <= kMax)

% OUTPUT: -----------------------------------------------------------------
% flag: Integer value indicating which domain does the given sub-cell 
%       belong to. Following encoding is used:
%       0 = physical, 1 = fictitious, 2 = cut

%% Initialization
% Number of different geometry sets
nGeo = numel(Geometry.type);

%% Test 1
% The algorithm first tests the four corners of the given sub-cell. Once 
% there exists two corners with different cornerFlags, it is certain that
% the sub-cell is cut and the entire function can be terminated. The figure
% below depicts such a case.
%
%   x-----------x
%   |           |
%   |   ____    /
%   |  |    \__/|
%   |  |        |
%   x--|--------o

for t = 1:4
    
    % Shape functions of the 4-node element 
	[N] = SF_4nodeElem(xiCoord(t),etaCoord(t));
    
    % Global coordinates (visualization of the decomposition)
    x = dot(N,xCoord(1:4));
    y = dot(N,yCoord(1:4));

    % Check for the given corner which domain does it belong to
    [cornerFlag,IdxDomain,Hole] = PointInMembershipTest(x,y,Geometry);
    if cornerFlag == 0
        IdxDomain = -1;
    end
    
    if t == 1
        % Set reference value equal to the flag of the first corner
        refFlag = cornerFlag;
        refDomain = IdxDomain;
        refHole = Hole;
        if Hole == 0
            % Inclusion --> physical domain
            flag = 0;
        elseif Hole == 1
            % Void/cutout/hole --> fictitious domain
            if cornerFlag == 1
                flag = 1;
            else
                flag = 0;
            end
        elseif Hole == 2
            % Physical domain
            if cornerFlag == 1
                flag = 0;
            else
                flag = 1;
            end
        end
        
    end
    
    if refFlag ~= cornerFlag || refDomain ~= IdxDomain
        % If any of the corner flags of the next nodes does not match
        % the reference value, the sub cell is cut.
        flag = 2;
        break;
    end
end

% If the variable 'flag' haven't been created, no two corner were found
% with different corner flags. In this case, the given cell undergo further
% procedures in section "Test 2" in order to test, whether it is cut or not.

%% Test 2
% If all corners have the same flag, it does not necessarily mean that the
% sub-cell is not cut. The entire boundary can be located within the sub-
% cell, or the boundary may cuts only a single edge of the sub-cell. The
% figures below depict examples for both cases.

%       x-----------x       x-----------x
%       |   _____   |       |  ______   |
%       |  /     |  |       |  \     \  |
%       | /_     |  |       |   \    _| |
%       |   \____|  |       |   |   |   |
%       x-----------x       x---|---|---x
%
% Consequently, these further tests that can resolve such cases result in a
% higher computational time!

for i = 1:nGeo
    % Vertices describing the geometry
    points = Geometry.points{i};
    
    if ~exist('flag','var')

        % Since the cells in the FCM undergo a geometry mapping with a
        % constant Jacobian, the mapping from the local to global space can
        % be carried out easily without the need of any shape functions and
        % matrix-vector operations.
        x1 = linearScalarMap(xiCoord(1),-1,1,min(xCoord), max(xCoord));
        y1 = linearScalarMap(etaCoord(1),-1,1,min(yCoord), max(yCoord));
        x2 = linearScalarMap(xiCoord(2),-1,1,min(xCoord), max(xCoord));
        y2 = linearScalarMap(etaCoord(2),-1,1,min(yCoord), max(yCoord));
        x3 = linearScalarMap(xiCoord(3),-1,1,min(xCoord), max(xCoord));
        y3 = linearScalarMap(etaCoord(3),-1,1,min(yCoord), max(yCoord));
        x4 = linearScalarMap(xiCoord(4),-1,1,min(xCoord), max(xCoord));
        y4 = linearScalarMap(etaCoord(4),-1,1,min(yCoord), max(yCoord));

        % Bounding box of the cell
        max_x = max([x1 x2 x3 x4]);
        min_x = min([x1 x2 x3 x4]);
        max_y = max([y1 y2 y3 y4]);
        min_y = min([y1 y2 y3 y4]);

        if ~isempty(points)
            % Test, whether there are any vertices of the boundary defined
            % by a polygonal chain in "points", that are located within the
            % given sub-cell.
            [in_val,~] = inpolygon(points(:,1),points(:,2),...
                [min_x; max_x; max_x; min_x],...
                [min_y; min_y; max_y; max_y]);

            if any(in_val)
                % At least one vertex is in the sub-cell and therefore it
                % must be cut and we assign the corresponding flag.
                flag = 2;
            else
                % No vertex is in the sub-cell and therefore it is either
                % physical or fictitious, depending on the reference corner
                % flag.
                if refHole == 0
                    flag = 0;
                elseif refHole == 1
                    flag = 1;
                elseif refHole == 2
                    flag = 0;
                end
            end

        elseif isempty(points)
            % Create nDiv^2 grid points in given sub-cell
            p = ceil(5-(3*k)/kMax);
            nDiv = (2^p)+1;

            [Coor] = ElementBackgroundGrid([min_x; max_x; max_x; min_x],...
            [min_y; min_y; max_y; max_y], nDiv);

            % Boundary points
            % Test, which domains do the grid points belong to
            [cornerFlag] = PointInMembershipTest(Coor.Boundary(:,1),...
                           Coor.Boundary(:,2),Geometry);

            if any(cornerFlag) && sum(cornerFlag) < 4*(nDiv-1)
                flag = 2;
                return;
            end

            % Interior points
            [cornerFlag] = PointInMembershipTest(Coor.Interior(:,1),...
                           Coor.Interior(:,2),Geometry);

            if any(cornerFlag) && sum(cornerFlag) < (nDiv-2)*(nDiv-2)
                flag = 2;
            else
                % If the vector cornerFlag contains only ones or only
                % zeros, the sub-cell is not cut, but either physical or
                % fictitious, depending on the reference corner flag.
                if refHole == 0
                    flag = 0;
                elseif refHole == 1
                    flag = 1;
                elseif refHole == 2
                    flag = 0;
                end
            end
        end
    end
end

end