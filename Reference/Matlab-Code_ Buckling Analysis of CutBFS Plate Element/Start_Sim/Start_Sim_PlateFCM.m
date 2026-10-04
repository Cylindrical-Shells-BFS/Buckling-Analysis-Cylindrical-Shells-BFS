%% Analysis of plates - Finite Cell Method
%-------------------------------------------------------------------------%
% FE analysis of plate bending using Bogner-Fox-Schmit shape functions, 
% which guarantee C^1 inter-element continuity.
%
% Theoretical background:
% Yu. Vetyukov, Finite element modeling of Kirchhoff-Love shells as smooth
% material surfaces, ZAMM-Journal of Applied Mathematics and Mechanics,
% Vol. 94, Issue 1-2, 2014, pp. 150-163
%
% Yu. Vetyukov, Nonlinear mechanics of thin-walled structures: asymptotics,
% direct approach and numerical analysis, Springer, 2014
%
% Plate and shell elements (general):
% Zienkiewicz, O. C. & Taylor, R. L.: The Finite Element Method
% Volume 2: Solid Mechanics, Butterworth Heinemann, 2000
%
% Mathematica implementation:
% (c) 2019 Yury Vetyukov, TU Wien, Austria
%
% Matlab implementation:
% (c) 2019 Sascha Duczek, UNSW Sydney, Australia
%-------------------------------------------------------------------------%
%% Settings
clc
clearvars
close all
format long

% For easier debugging
dbstop if error

% Collect subfolders as string
pa = genpath('.\subroutines');
addpath(pa); 

% Output on the command window
Line= '-';
LineDisp = repmat(Line,1,75);
%--------------------------------------------------------------------------
disp([repmat(Line,1,28) 'Start of Simulation!' repmat(Line,1,27)])
disp('Bogner-Fox-Schmit plate element -  Finite cell method')
disp(LineDisp)

%% Material properties and parameters
disp('Define material properties and parameters')
% Define the material properties for the base material and possible
% inclusions
% Number of degrees of freedom per nodes
ndof = 4;
% Young's modulus
E{1} = 2.1e11;      % Base material
% E{2} = 7.0e10;      % Inclusion material #1
% Poisson's ratio
nu{1} = 0.3;        % Base material
% nu{2} = 0.3;        % Inclusion material #1
% Mass density
rho{1} = 7800;      % Base material
% rho{2} = 2700;      % Inclusion material #1
% Number of different materials
nMat = numel(E);
% Dimensions
lx = 1;
ly = 1;
% Number of elements in x- and y-direction (even number --> point force)
m = 9;
nx = 2^m;
ny = 2^m;
nElem = nx*ny;
% Plate thickness
h = 1e-3;
% Plate stiffness parameters
D = cell(nMat,1);
for i = 1:nMat
    D{i} = (E{i}*h^2)/(12*(1-nu{i}^2))*[1,nu{i},0;nu{i},1,0;0,0,(1-nu{i})/2];
end
% Surface load (pressure)
flag_func = 0;  % Pressure amplitude is a function of space: yes (1)
press = -1;

% Preconditioning (diagonal)
flag_precond = true;

%% FCM: Settings
% Parameters of the fictitious domain (see LevelSetFunctionGeometry.m)
%---------------------------Geometry Definition---------------------------%
% The geometry can be a void/hole, a physical domain (outer boundary of the
% computational domain), or an inclusion. The following values have to be
% assigned:'Fictitious', 'Physical', or 'Inclusion'
Geometry.domain{1} = 'Fictitious';
% Geometry.domain{2} = 'Fictitious';
% Geometry.domain{3} = 'Physical';

% The geometry can be either described by an implicit function or by a list
% of points (polygonal line): 'Implicit' or 'Polygon'
Geometry.type{1}   = 'Implicit';
% Geometry.type{2}   = 'Implicit';
% Geometry.type{3}   = 'Implicit';

% Shape of the domain; possible are a circle, an ellipse, a quadrangle, a
% rectangle or a line: 'Circle', 'Ellipse', 'Rectangle', 'Quadrangle', or
% 'QuadrangleRot' (please LevelSetFunctionGeometry.m for the definition of
% the parameters!!!)
Geometry.shape{1}  = 'Circle';
% Geometry.shape{2}  = 'Circle';
% Geometry.shape{3}  = 'Circle';

% Geometry parameters describing the defined shape
Geometry.param{1}  = [0.5, 0.5, 0.15];
% Geometry.param{2}  = [0.65, 0.65, 0.15];
% Geometry.param{3}  = [0.5,  0.5,  0.4];
% hold on
% th = 0:pi/500:2*pi;
% xunit = Geometry.param{3}*cos(th) + Geometry.param{1};
% yunit = Geometry.param{3}*sin(th) + Geometry.param{2};
% plot(xunit,yunit,'r','LineWidth',2);

% Index to the material properties that should be used
Geometry.material{1} = 1;
% Geometry.material{2} = 1;
% Geometry.material{3} = 1;

% Index to the base material (physical domain is defined by the boundary of
% the finite element mesh)
Geometry.base_mat    = 1;

% If a polygonal domain is chosen as type, the vertices are defined in the
% field points
Geometry.points{1} = [];
% Geometry.points{2} = [];
% Geometry.points{3} = [];

% Define the standard acceleration due to gravity (standard gravity)
Geometry.g = -9.81;
%-------------------------------------------------------------------------%
% Indicator function for FCM
alpha_FCM = 1e-10;

% Refinement level for the numerical integration of the cut cells (finite
% cells that are intersected/cut by the physical boundary of the domain)
level_k = 4;

% Standard (0) or Boolean (1) quadtree decomposition
flagBFCM = 0;

%% Mesh
disp('Mesh generation')
[TopoN,~,nNodes,~,xCoor,yCoor] = MeshGenerator(ndof,nx,ny,lx,ly);

%% FCM: Determine cut and fictitious elements
disp('FCM: Determine cut and fictitious elements')
% Set the bounding box of the fictitious domain (can be multiple)
[~,~,~,BdBox] = LevelSetFunctionGeometry(0,0,Geometry);
% Check which nodes are within the bounding box
[ElemBdBox] = ElementsInBdBox(TopoN,xCoor,yCoor,BdBox);

% Fictitious and cut elements
[CutElem,FictElem] = Cut_Elements(TopoN,xCoor,yCoor,ElemBdBox,Geometry);

% Physical elements (cut elements + standard elements)
PhysElem = setdiff((1:nElem),cell2mat(FictElem));

% Delete fictitious elements and set up a new nodal topology matrix
TopoN = TopoN(PhysElem,:);

% Delete nodes
NodesDel = setdiff((1:nNodes),unique(TopoN));
nNodesDel = numel(NodesDel);

% Nodes that are kept
NodesKeep = setdiff((1:nNodes),NodesDel);
% New nodal coordinates
xCoor = xCoor(NodesKeep);
yCoor = yCoor(NodesKeep);

% New number of nodes in the model
nNodes = nNodes - nNodesDel;

% New number of degrees of freedom
nDOF = nNodes*ndof;

% Adjust nodal topology data
for i = 1:nNodesDel
    % Nodes with nodal numbers that are greater than the nodal number of
    % the deleted node have to be adjusted.
    ID = TopoN > NodesDel(i);
    % The nodal numbers are reduced by 1
    TopoN(ID) = TopoN(ID) - 1;
    NodesDel = NodesDel - 1;
end

% Number of void regions/inclusions
nGeo = numel(CutElem);
for i = 1:nGeo
    % Adjust the numbers of the cut elements
    VecFictElem = unique(cell2mat(FictElem));
    nFictElem   = numel(VecFictElem);
    for j = 1:nFictElem
        % Elements with element numbers that are greater than those of the
        % fictitious elements have to be adjusted.
        ID = CutElem{i} > VecFictElem(j);
        % The element numbers are reduced by 1
        CutElem{i}(ID) = CutElem{i}(ID) - 1;
        VecFictElem = VecFictElem - 1;
    end
end

% New DOF topology data
[TopoD] = TopologyDOF(TopoN);

% Plot the mesh
FE_Grid(xCoor,yCoor,TopoN);

%% Stiffness, mass matrices and load vector (self-weight)
disp('Compute stiffness, mass matrices and load vector')
[K,M,Fg] = AssemblyProcedureFCM(TopoN,TopoD,CutElem,xCoor,yCoor,h,D,rho,...
           level_k,alpha_FCM,flagBFCM,Geometry);

%% Point force
% Find node at (lx/2,ly/2)
idx = xCoor == lx/2;
idy = yCoor == ly/2;
id  = and(idx,idy);
% Corresponding displacement degree of freedom
id = find(id)*ndof-3;

% Global force vector
Fp = zeros(nDOF,1);
% Same magnitude of the force (self-weight)
Fp(id) = rho{1}*Geometry.g*lx*ly*h;

%% Distributed (surface) load
% Elements with pressure loading
Elset_p = (1:size(TopoN,1));
% Compute the pressure loading
[Fpres] = PressureLoad(TopoN,TopoD,Elset_p,CutElem,press,xCoor,yCoor,...
          level_k,alpha_FCM,flagBFCM,flag_func,Geometry);

%% Boundary conditions
% Clamped boundary conditions
%--------------------------------------------------------------------------
% % Find bottom boundary edge
% idx_be = yCoor == 0;
% % Find right boundary edge
% idx_re = xCoor == lx;
% % Find top boundary edge
% idx_te = yCoor == ly;
% % Find left boundary edge
% idx_le = xCoor == 0;
% 
% % Inactive nodes
% InactiveNodes = or(idx_be,or(idx_re,or(idx_te,idx_le)));
% % Active nodes
% ActiveNodes = (1:nNodes)';
% ActiveNodes = ActiveNodes(~InactiveNodes);
% 
% % Active degrees of freedom
% DOF_aux1 = kron(ActiveNodes,ones(1,ndof)*ndof);
% DOF_aux2 = kron(ones(size(ActiveNodes)),(ndof-1:-1:0));
% ActiveDOFs = unique(DOF_aux1 - DOF_aux2);

% Simply-supported boundary conditions
%--------------------------------------------------------------------------
% Find bottom boundary edge
idx_be = yCoor == 0;
% Find right boundary edge
idx_re = xCoor == lx;
% Find top boundary edge
idx_te = yCoor == ly;
% Find left boundary edge
idx_le = xCoor == 0;

% Index to the first degree of freedom of the corresponding node,
%  i.e., displacement degree of freedom (U)
id_be = (find(idx_be)-1)*ndof+1;
id_re = (find(idx_re)-1)*ndof+1;
id_te = (find(idx_te)-1)*ndof+1;
id_le = (find(idx_le)-1)*ndof+1;

% Combine the displacement degrees of freedom of all edges into one vector
% to constrain them
id1 = union(id_be, id_re);
id1 = union(id1, id_te);
id1 = union(id1, id_le);

% Constrain the first derivatives parallel to the edges
id_be2 = (find(idx_be)-1)*ndof+2;
id_re2 = (find(idx_re)-1)*ndof+3;
id_te2 = (find(idx_te)-1)*ndof+2;
id_le2 = (find(idx_le)-1)*ndof+3;

% Combine the rotational degrees of freedom of all edges into one vector
% to constrain them
id2 = union(id_be2, id_re2);
id2 = union(id2, id_te2);
id2 = union(id2, id_le2);

% Collect all DOFs to constrain
BC = [id1; id2];
BC = unique(BC);

% Define active DOFs
ActiveDOFs = (1:nDOF)';
% remove selected DOFs from ActiveDOFs
ActiveDOFs(BC) = [];

%% Solution - Statics
disp('Static solution')
%% Gravity load
% uC = K(ActiveDOFs,ActiveDOFs)\Fg(ActiveDOFs);
%% Point load
% uC = K(ActiveDOFs,ActiveDOFs)\Fp(ActiveDOFs);
%% Pressure (surface) load
uC = K(ActiveDOFs,ActiveDOFs)\Fpres(ActiveDOFs);

% Displacement field
U = zeros(nDOF,1);
U(ActiveDOFs) = uC;

% Strain energy
SE = 1/2*U'*K*U;

if flag_precond
% Preconditioning
[Kcond,Fcond,D] = PreconditioningLSE(K(ActiveDOFs,ActiveDOFs),Fpres(ActiveDOFs));
y = Kcond\Fcond;
uC2 = D'*y;

% Displacement field
U2 = zeros(nDOF,1);
U2(ActiveDOFs) = uC2;

% Strain energy
SE2 = 1/2*U2'*K*U2;
end

%% Visualization
disp('Visualization of the displacement field')
% Number of subdivisions of each element for the visualization
ndiv = 1;
% Plot the displacement field
PlotResults(TopoN,TopoD,xCoor,yCoor,U,ndiv);

%% Solution - Modal
% disp('Eigenvalue problem')
% % Number of eigenvalues requested
% nEV = 15;
% if nEV > numel(ActiveDOFs)
%     nEV = numel(ActiveDOFs);
% end
% % Solve the eigenvalue problem and compute the nEV smallest eigenvalues
% [MS,EV] = eigs(K(ActiveDOFs,ActiveDOFs),M(ActiveDOFs,ActiveDOFs),nEV,'sm');
% 
% % Determine if the eigenvalues are real
% if ~isreal(EV)
%     EV = real(EV);
%     MS = real(MS);
%     % Display warning message
%     warning('Eigenvalues and mode shapes are complex!')
% end
% 
% % Add the fixed degrees of freedom
% ModeShapes = zeros(nDOF,nEV);
% ModeShapes(ActiveDOFs,:) = MS;
% % Eigenfrequencies
% EF = sqrt(diag(EV))/(2*pi);
% 
% % Plot mode shape
% disp('Visualization of a mode shape')
% PlotResults(TopoN,TopoD,xCoor,yCoor,ModeShapes(:,nEV),ndiv);

disp(LineDisp)
disp([repmat(Line,1,28) 'End of Simulation!' repmat(Line,1,28)])
disp(LineDisp)