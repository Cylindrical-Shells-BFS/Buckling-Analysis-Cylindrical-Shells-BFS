%% Assembly of the system matrices and load vector - FCM
%-------------------------------------------------------------------------%
%
% (c) 04/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% The overall (global) stiffness, mass matrices and load vector are
% assembled.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:    Nodal topology (connectivity/coincidence) matrix
% TopoD:    DOF topology (connectivity/coincidence) matrix
% CutElem:  Elements that are intersected by the boundary of the
%           computational domain
% xCoor:    Global x-coordinates of the nodes
% yCoor:    Global y-coordinates of the nodes
% h:       	Plate thickness
% D:      	Plate elasticity matrix
% rho:    	Mass density
% k:        Number of hierarchic subdivisions of cut elements
% alpha:    Value of the indicator function (avoid ill-conditioning)
% flagBFCM: Quadtree decomposition method
% Geometry: Structure that contains geometric- and meta-information about
% 			the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% K:        Global stiffness matrix
% M:        Global mass matrix
% F:        Global load vector (due to gravity load)
%
%-------------------------------------------------------------------------%
function [K,M,F] = AssemblyProcedureFCM(TopoN,TopoD,CutElem,xCoor,yCoor,...
                   h,D,rho,k,alpha,flagBFCM,Geometry)

% Start time measurement
tic;

% Number of elements
nElem = size(TopoD,1);

% Number of degrees of freedom
nDOF = max(max(TopoD));

% Save results in cell arrays
Kcell     = cell(nElem,1);
Mcell     = cell(nElem,1);
Fcell     = cell(nElem,1);
topo_row  = cell(nElem,1);
topo_col  = cell(nElem,1);

% Number of void regions/inclusions
nD = numel(CutElem);

% Initialize the waitbar
anz = waitbar(0,'Computing FE element matrices');

% If inclusions are present in the model the correct material properties
% have to be determined
nGeo = numel(Geometry.type);
flag_inclusion = 0;
for i = 1:nGeo
    if strcmp(Geometry.domain{i},'Inclusion')
        flag_inclusion = 1;
        break;
    end
end

% Material properties of the base material
Dbase = D{Geometry.base_mat};
rhobase = rho{Geometry.base_mat};

% Acceleration due to gravity
g = Geometry.g;

% Compute the element matrices
for i = 1:nElem
    % Nodal coordinates
    xElem = xCoor(TopoN(i,:));
    yElem = yCoor(TopoN(i,:));
    
    % Determine if the current element is a cut or standard element
    for j = 1:nD
        ID = find(CutElem{j} == i, 1);
        if ~isempty(ID)
            break;
        end
    end
    
    if isempty(ID)  % FEM
        if flag_inclusion == 0
            % Elemental stiffness matrix
            [Ke] = StiffMat(xElem,yElem,h,Dbase);
            % Elemental mass matrix
            [Me] = MassMat(xElem,yElem,h,rhobase);
            % Elemental force vector due to gravity
            [Fe] = Gravity(xElem,yElem,h,rhobase,g);
        else
            % Elemental stiffness matrix
            [Ke] = StiffMatInclusion(xElem,yElem,h,D,Geometry);
            % Elemental mass matrix
            [Me] = MassMatInclusion(xElem,yElem,h,rho,Geometry);
            % Elemental force vector due to gravity
            [Fe] = GravityInclusion(xElem,yElem,h,rho,Geometry);
        end
    else            % FCM
        % Check whether the cut cell is belonging to an inclusion or a
        % void/physical region. The difference is that in the first case
        % two different material properties are used within the cell.
        % Therefore the standard quadtree decomposition is used (no
        % distinct advantage of using the Boolean quadtree decomposition
        % twice).
        [FCMval,IdxDomain,~] = PointInMembershipTest(xElem,yElem,Geometry);
        idx = FCMval == 1;
        ID = unique(IdxDomain(idx));
        % Number of domains the element is intersected by
        nIntersect = numel(ID);
        
        % The centroid of the cell is located in one of the geometric
        % domains that have been defined in the main file.
        flagBFCMi = flagBFCM;
        for j = 1:nIntersect
            domain = Geometry.domain{ID(j)};
            if strcmp(domain,'Inclusion')
                flagBFCMi = 0;
                break;
            end
        end
        
        % Standard or Boolean quadtree decomposition
        if flagBFCMi == 1
            % Quadtree decomposition of the cut element (Boolean FCM)
            [XI,ETA,alphaB,tree,~] = BooleanQuadtree(xElem,yElem,k,...
                                     Geometry);
            % Struct with all quadtree data
            QT.XI     = XI;
            QT.ETA    = ETA;
            QT.alphaB = alphaB;
            QT.tree   = tree;
            % displayBooleanQuadtree(QT)
            
            % Elemental stiffness matrix
            [Ke] = StiffMatBFCM(xElem,yElem,h,D,QT,alpha,flagBFCMi,...
                   Geometry);
            % Elemental mass matrix
            [Me] = MassMatBFCM(xElem,yElem,h,rho,QT,alpha,flagBFCMi,...
                   Geometry);
            % Elemental force vector due to gravity
            [Fe] = GravityBFCM(xElem,yElem,h,rho,QT,alpha,flagBFCMi,....
                   Geometry);
        else
            % Quadtree decomposition of the cut element (Standard quadtree)
            [XI,ETA,domains,tree,~] = modQuadtree(xElem,yElem,k,...
                                      Geometry,true);
            % Struct with all quadtree data
            QT.XI      = XI;
            QT.ETA     = ETA;
            QT.domains = domains;
            QT.tree    = tree;
            % displayQuadtree(QT)
            
            % Elemental stiffness matrix
            [Ke] = StiffMatFCM(xElem,yElem,h,D,QT,alpha,Geometry);
            % Elemental mass matrix
            [Me] = MassMatFCM(xElem,yElem,h,rho,QT,alpha,Geometry);
            % Elemental force vector due to gravity
            [Fe] = GravityFCM(xElem,yElem,h,rho,QT,alpha,Geometry);
        end
    end
    
	% Save elemental matrices
	Kcell{i}    = reshape(Ke,[],1);
	Mcell{i}    = reshape(Me,[],1);
	Fcell{i}    = Fe;
	topo_row{i} = reshape(repmat(TopoD(i,:),16,1),[],1);
  	topo_col{i} = reshape(repmat(TopoD(i,:)',16,1),[],1);
	
	% Waitbar
	waitbar(i/nElem,anz,['Computing FE element matrices: ' ...
    num2str(i/nElem*100,'%5.1f') '%'])
end

% Close the waitbar
close(anz)

% Concatenate all the cell arrays into a matrix
k   = cell2mat(Kcell);
m   = cell2mat(Mcell);
row = cell2mat(topo_row);
col = cell2mat(topo_col);

% Save the result vectors in sparse matrices
% Global stiffness matrix
K = sparse(row,col,k,nDOF,nDOF);
% % % % Delete small components
% % % kmax = max(max(abs(K)));
% % % K(abs(K)<kmax*1e-10) = 0;

% Global mass matrix
M = sparse(row,col,m,nDOF,nDOF); 
% % % % Delete small components
% % % mmax = max(max(abs(M)));
% % % M(abs(M)<mmax*1e-10) = 0;

% Global force vector
F = zeros(nDOF,1);
for i = 1:nElem
    F(TopoD(i,:),1) = F(TopoD(i,:),1) + Fcell{i};
end
fmax = max(abs(F));
F(abs(F)<fmax*1e-10) = 0;

% Stop time measurement
toc

end