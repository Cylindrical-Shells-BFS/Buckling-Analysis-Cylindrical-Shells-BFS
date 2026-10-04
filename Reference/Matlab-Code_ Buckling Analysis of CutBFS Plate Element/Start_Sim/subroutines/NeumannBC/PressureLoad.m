%%  Pressure load - FCM
%-------------------------------------------------------------------------%
%
% (c) 04/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the pressure load.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:        Nodal topology (connectivity/coincidence) matrix
% TopoD:        DOF topology (connectivity/coincidence) matrix
% Elset:        Element set for the pressure load
% CutElem:      Elements that are intersected by the boundary of the
%               computational domain
% mag:          Magnitude of the pressure load
% xCoor:        Global x-coordinates of the nodes
% yCoor:        Global y-coordinates of the nodes
% k:            Number of hierarchic subdivisions of cut elements
% alpha:        Value of the indicator function (avoid ill-conditioning)
% flagBFCM:     Quadtree decomposition method
% flag_func: 	Constant pressure load or functional dependency
% Geometry:     Structure that contains geometric- and meta-information 
% 				about the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% F:            Global load vector (due to pressure load)
%
%-------------------------------------------------------------------------%
function [F] = PressureLoad(TopoN,TopoD,Elset,CutElem,mag,xCoor,yCoor,k,...
               alpha,flagBFCM,flag_func,Geometry)

% Number of elements with pressure load
nElem = numel(Elset);

% Number of degrees of freedom
nDOF = max(max(TopoD));

% Save results in cell arrays
Fcell     = cell(nElem,1);

% Number of void regions/inclusions
nD = numel(CutElem);

% Compute the element matrices
for i = 1:nElem
    % Nodal coordinates
    xE = xCoor(TopoN(Elset(i),:));
    yE = yCoor(TopoN(Elset(i),:));
    
    % Determine if the current element is a cut or standard element
    for j = 1:nD
        ID = find(CutElem{j} == Elset(i), 1);
        if ~isempty(ID)
            break;
        end
    end
    
    if isempty(ID)  % FEM
        % Elemental force vector due to pressure
        [Fe] = Pressure((1:4),(1:16),xE,yE,1,mag,flag_func);
    else            % FCM
         % Check whether the cut cell is belonging to an inclusion or a
        % void/physical region. The difference is that in the first case
        % two different material properties are used within the cell.
        % Therefore the standard quadtree decomposition is used (no
        % distinct advantage of using the Boolean quadtree decomposition
        % twice).
        [FCMval,IdxDomain,~] = PointInMembershipTest(xE,yE,Geometry);
        idx = find(FCMval == 0,1,'first');
        
        % The centroid of the cell is located in one of the geometric
        % domains that have been defined in the main file.
        if ~isempty(idx)
            domain = Geometry.domain{IdxDomain(idx)};
            if strcmp(domain,'Inclusion')
                flagBFCM = 0;
            end
        end
        
        % Standard or Boolean quadtree decomposition
        if flagBFCM == 1
            % Quadtree decomposition of the cut element (Boolean FCM)
            [XI,ETA,alphaB,tree,~] = BooleanQuadtree(xE,yE,k,Geometry);
            % Struct with all quadtree data
            QT.XI     = XI;
            QT.ETA    = ETA;
            QT.alphaB = alphaB;
            QT.tree   = tree;
            % displayBooleanQuadtree(QT)

            % Elemental force vector due to pressure loading
            [Fe] = PressureBFCM(xE,yE,QT,mag,alpha,flag_func,Geometry);
        else
            % Quadtree decomposition of the cut element (Standard quadtree)
            [XI,ETA,domains,tree,~] = modQuadtree(xE,yE,k,Geometry,true);
            % Struct with all quadtree data
            QT.XI      = XI;
            QT.ETA     = ETA;
            QT.domains = domains;
            QT.tree    = tree;
            % displayQuadtree(QT)

            % Elemental force vector due to pressure loading
            [Fe] = PressureFCM(xE,yE,QT,mag,alpha,flag_func,Geometry);
        end
    end
    
	% Save elemental matrices
	Fcell{i}    = Fe;
end

% Global force vector
F = zeros(nDOF,1);
for i = 1:nElem
    F(TopoD(i,:),1) = F(TopoD(i,:),1) + Fcell{i};
end
fmax = max(abs(F));
F(abs(F)<fmax*1e-10) = 0;

end