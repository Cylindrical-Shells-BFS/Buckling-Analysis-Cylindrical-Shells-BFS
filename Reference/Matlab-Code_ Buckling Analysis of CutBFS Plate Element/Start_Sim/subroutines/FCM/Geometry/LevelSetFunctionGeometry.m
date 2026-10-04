%% Level set function
%-------------------------------------------------------------------------%
%
% (c) 05.09.2016, Sascha Duczek, OvGU Magdeburg, Germany
%
% Level set function for the geometry description of the computational 
% domain.
%
%-------------------------------------------------------------------------%
%% Input
%
% [x,y]:            Coordinates of evaluation points (compute value of the
% 					level set function)
% Geometry:         Structure that contains geometric- and meta-
%                   information about the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% LevelSetValue:	Values of the level set functions at the evaluation
% 					points (< 0: inside; > 0: outside)
% IdxDomain:        Domain with the shortest distance to the evaluation
%					point (if several cutout/inclusions are involved)
% Hole:        		Inclusion (0), void/cutout(1), physical domain (2)
% BdBox:            Bounding box(es) of the defined geometries
%   * BdBox.coor    Coordinates of the diagonal describing the rectangular
%                   bounding box
%   * BdBox.type    Type of the bounding box; defined for a 'Physical' or
%                   'Fictitious' region (domain); if the domain is a
%                   physical domain the bounding box is actually the
%                   maximum rectangle that can be inscribed in the geometry
%
%-------------------------------------------------------------------------%
function [LevelSetValue,IdxDomain,Hole,BdBox] = ...
          LevelSetFunctionGeometry(x,y,Geometry)

% Number of different geometry sets
nGeo = numel(Geometry.type);

% Initialize variables
LSV = cell(1,nGeo);
Hole = zeros(nGeo,1);

for i = 1:nGeo
    % Type of the domain: 'Implicit' or 'Polygon'
    if strcmp(Geometry.type{i},'Polygon')
        xPoints = Geometry.points{i}(:,1);
        yPoints = Geometry.points{i}(:,2);
        
        % Bounding box
        BdBox.coor{i} = [min(xPoints),min(yPoints);...
                         max(xPoints),max(yPoints)];
        BdBox.type{i} = Geometry.domain{i};
        
        % Test whether the evaluation points are inside the polygon
        [LevelSetValue1,~] = DistanceFunctionPolygon(x,y,xPoints,...
                             yPoints,true);
                         
        % Inclusion (=0), hole (=1), or physical domain (=2)
        if  strcmp(Geometry.domain{i},'Inclusion')
            Hole1 = 0;
        elseif strcmp(Geometry.domain{i},'Fictitious')
            Hole1 = 1;
        elseif  strcmp(Geometry.domain{i},'Physical')
            Hole1 = 2;
        end
    elseif strcmp(Geometry.type{i},'Implicit')
        if strcmp(Geometry.shape{i},'Circle')
            % Centroid of the circle
            xc  = Geometry.param{i}(1);
            yc  = Geometry.param{i}(2);
            % Radius of the circle
            rad = Geometry.param{i}(3);

            % Bounding box
            BdBox.coor{i} = [xc-rad,yc-rad;xc+rad,yc+rad];
            BdBox.type{i} = Geometry.domain{i};
%             % Inner square
%             BdBox.coor{i} = [xc-rad/sqrt(2),yc-rad/sqrt(2);...
%                           	 xc+rad/sqrt(2),yc+rad/sqrt(2)];
%             BdBox.type{i} = Geometry.domain{i};

            [LevelSetValue1,~] = DistanceFunctionCircle([x,y],xc,yc,rad);
            
            % Inclusion (=0), hole (=1), or physical domain (=2)
            if  strcmp(Geometry.domain{i},'Inclusion')
                Hole1 = 0;
            elseif strcmp(Geometry.domain{i},'Fictitious')
                Hole1 = 1;
            elseif  strcmp(Geometry.domain{i},'Physical')
                Hole1 = 2;
            end
        elseif strcmp(Geometry.shape{i},'Ellipse')
            % Centroid of the circle
            xc    = Geometry.param{i}(1);
            yc    = Geometry.param{i}(2);
            % Dimension of the major axis
            a     = Geometry.param{i}(3);
            % Dimension of the minor axis
            b     = Geometry.param{i}(4);
            % Angle (counted counterclockwise from x-axis)
            gamma = Geometry.param{i}(5);

            % Bounding box
            BdBox.coor{i} = [xc-a,yc-a;xc+a,yc+a];
            BdBox.type{i} = Geometry.domain{i};

            [LevelSetValue1,~] = DistanceFunctionEllipse([x,y],xc,yc,a,...
                                 b,gamma);
            
            % Inclusion (=0), hole (=1), or physical domain (=2)
            if  strcmp(Geometry.domain{i},'Inclusion')
                Hole1 = 0;
            elseif strcmp(Geometry.domain{i},'Fictitious')
                Hole1 = 1;
            elseif  strcmp(Geometry.domain{i},'Physical')
                Hole1 = 2;
            end
        elseif strcmp(Geometry.shape{i},'QuadrangleRot')
            % Coordinates of the vertices of the quadrangle
            xp = [Geometry.param{i}(1); Geometry.param{i}(2); ...
                  Geometry.param{i}(3); Geometry.param{i}(4)];
            yp = [Geometry.param{i}(5); Geometry.param{i}(6); ...
                  Geometry.param{i}(7); Geometry.param{i}(8)];
        	% Angle (counted counterclockwise from x-axis)
            gamma = Geometry.param{i}(9);
            
            % Translate the centroid of the quadrangle to the origin of the
            % coordinate system
            PPtrans = [xp - sum(xp)/4, yp - sum(yp)/4];

            % Rotate the coordinate system
            PProt = RotateCoos(PPtrans,gamma);

            % Translate the points back
            PP = [PProt(:,1) + sum(xp)/4, PProt(:,2) + sum(yp)/4];
            
            % Bounding box
            BdBox.coor{i} = [min(PP(:,1)),min(PP(:,2));...
                             max(PP(:,1)),max(PP(:,2))];
            BdBox.type{i} = Geometry.domain{i};

            [LevelSetValue1,~] = DistanceFunctionQuadrangle([x,y],PP);
            
            % Inclusion (=0), hole (=1), or physical domain (=2)
            if  strcmp(Geometry.domain{i},'Inclusion')
                Hole1 = 0;
            elseif strcmp(Geometry.domain{i},'Fictitious')
                Hole1 = 1;
            elseif  strcmp(Geometry.domain{i},'Physical')
                Hole1 = 2;
            end
        elseif strcmp(Geometry.shape{i},'Quadrangle')
            % Coordinates of the vertices of the quadrangle
            xp = Geometry.param{i}(:,1);
            yp = Geometry.param{i}(:,2);

            % Bounding box
            BdBox.coor{i} = [min(xp),min(yp);max(xp),max(yp)];
            BdBox.type{i} = Geometry.domain{i};

            [LevelSetValue1,~] = DistanceFunctionQuadrangle([x,y],...
                                 [xp,yp]);
            
            % Inclusion (=0), hole (=1), or physical domain (=2)
            if  strcmp(Geometry.domain{i},'Inclusion')
                Hole1 = 0;
            elseif strcmp(Geometry.domain{i},'Fictitious')
                Hole1 = 1;
            elseif  strcmp(Geometry.domain{i},'Physical')
                Hole1 = 2;
            end
        elseif strcmp(Geometry.shape{i},'Rectangle')
            % Coordinates of the diagonal vertices
            x1  = Geometry.param{i}(1,1);
            y1  = Geometry.param{i}(1,2);
            x2  = Geometry.param{i}(2,1);
            y2  = Geometry.param{i}(2,2);

            % Bounding box
            BdBox.coor{i} = [x1,y1;x2,y2];
            BdBox.type{i} = Geometry.domain{i};

            [LevelSetValue1,~] = DistanceFunctionRectangle([x,y],x1,x2,...
                                 y1,y2);
            
            % Inclusion (=0), hole (=1), or physical domain (=2)
            if  strcmp(Geometry.domain{i},'Inclusion')
                Hole1 = 0;
            elseif strcmp(Geometry.domain{i},'Fictitious')
                Hole1 = 1;
            elseif  strcmp(Geometry.domain{i},'Physical')
                Hole1 = 2;
            end
        elseif strcmp(Geometry.shape{i},'Line')
            % Coordinates of the line vertices
            x1  = Geometry.param{i}(1);
            y1  = Geometry.param{i}(2);
            x2  = Geometry.param{i}(3);
            y2  = Geometry.param{i}(4);

            % Bounding box
            BdBox.coor{i} = [x1,y1;x2,y2];
            BdBox.type{i} = Geometry.domain{i};

            [LevelSetValue1,~] = DistanceFunctionLine([x,y],xi,x2,y1,y2);
                        
            % Inclusion (=0), hole (=1), or physical domain (=2)
            if  strcmp(Geometry.domain,'Inclusion')
                Hole1 = 0;
            elseif strcmp(Geometry.domain,'Fictitious')
                Hole1 = 1;
            elseif  strcmp(Geometry.domain,'Physical')
                Hole1 = 2;
            end
        end
        
    end
    % LevelSetValue specifies the distance to the closest boundary
    LSV{i} = LevelSetValue1;
    % Hole specifies if the ith domain a hole or an inclusion
    Hole(i) = Hole1;
end

% Combine the results for multiple domains
if nGeo == 1
    % LevelSetValue specifies the distance to the closest boundary
    LevelSetValue = LevelSetValue1;
    % Idx identifies the domain with the smallest distance to the
    % considered point --> can be used to determine the material properties
    % for this point
    IdxDomain = ones(size(x));
else
	% Boolean operations
    for i = 1:nGeo-1
        if i == 1
            % Level set values for domain i
            LevelSetValue1 = LSV{i};
        else
            LevelSetValue1 = LevelSetValue12;
        end
        % Level set value for domain i+1
        LevelSetValue2 = LSV{i+1};
        
        % Determine the LevelSetValues that are less than zero for both
        % domains
        LSV1_logical = LevelSetValue1 < 0;
        LSV2_logical = LevelSetValue2 < 0;
        
        % If both values are zero a different Boolean operation is used; if
        % a point is in more than one domain the smallest distance to the
        % boundary is of interest
        idxLSV1 = find(LSV1_logical);
        idxLSV2 = find(LSV2_logical);
        idx1 = intersect(idxLSV1,idxLSV2);
        idx2 = setdiff(idxLSV1,idxLSV2);
        idx3 = setdiff(1:numel(LSV1_logical),[idx1;idx2]);
        
        % Initialize result variable
        % Level set values
        LevelSetValue12 = zeros(size(LevelSetValue1));
        % IdxDomain identifies the domain with the smallest distance to the
        % considered point --> can be used to determine the material
        % properties for this point
        if i == 1
            IdxDomain = zeros(size(LevelSetValue1));
        end
        % Negative signs of the Level Set Values, i.e. the points are
        % located in two domains; smallest distance is now of interest to
        % determine in which domain the points are located.
        if~isempty(idx1)
            % Union of the two domains
            [LevelSetValue3,Idx] = ...
             DistanceFunctionBooleanOperationIntersection(...
             LevelSetValue1(idx1),LevelSetValue2(idx1));
             LevelSetValue12(idx1) = LevelSetValue3;
             if i == 1
                 IdxDomain(idx1) = Idx;
             else
                Idd = Idx == 2;
                IdxDomain(idx1(Idd)) = i+1;
            end
        end
        % Different signs of the Level Set Values; simple union operation
        % suffices; smaller value indicates the closest domain.
        if ~isempty(idx2)
            % Union of the two domains
            [LevelSetValue3,Idx] = ...
             DistanceFunctionBooleanOperationUnion(LevelSetValue1(idx2),...
             LevelSetValue2(idx2));
             LevelSetValue12(idx2) = LevelSetValue3;
             if i == 1
                 IdxDomain(idx2) = Idx;
             else
                Idd = Idx == 2;
                IdxDomain(idx2(Idd)) = i+1;
            end
        end
        % Positive signs of the Level Set Values; simple union operation
        % suffices; smaller value indicates the closest domain.
        if ~isempty(idx3)
            % Union of the two domains
            [LevelSetValue3,Idx] = ...
             DistanceFunctionBooleanOperationUnion(LevelSetValue1(idx3),...
             LevelSetValue2(idx3));
             LevelSetValue12(idx3) = LevelSetValue3;
             if i == 1
                 IdxDomain(idx3) = Idx;
             else
                Idd = Idx == 2;
                IdxDomain(idx3(Idd)) = i+1;
            end
        end
    end
    
    % Final level set value
    LevelSetValue = LevelSetValue12;
   
    % Inclusion, void, or physical domain
    Hole = Hole(IdxDomain);
end

end