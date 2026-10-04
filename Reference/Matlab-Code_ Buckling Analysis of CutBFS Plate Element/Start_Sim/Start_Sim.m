%% Analysis of plates
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
% Number of degrees of freedom per nodes
ndof = 4;
% Young's modulus
E = 2.1e11;
% Poisson's ratio
nu = 0.3;
% Mass density
rho = 7800;
% Dimensions
lx = 1;
ly = 1;
% Number of elements in x- and y-direction (even number --> point force)
nx = 10;
ny = 10;
% Plate thickness
h = 1e-3;
% Gravity
g = -9.8;
% Plate stiffness parameters
D = (E*h^2)/(12*(1-nu^2))*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
% Surface load (pressure)
flag_func = 1;  % Pressure amplitude is a function of space: yes (1)
press = -100;

%% Mesh
disp('Mesh generation')
[TopoN,TopoD,nNodes,nDOF,xCoor,yCoor] = MeshGenerator(ndof,nx,ny,lx,ly);

%% Stiffness, mass matrices and load vector (self-weight)
disp('Compute stiffness, mass matrices and load vector')
[K,M,Fg] = AssemblyProcedure(TopoN,TopoD,xCoor,yCoor,h,D,rho,g,1);

%% Point force
% Find node at (lx/2,ly/2)
idx = xCoor == lx/2;
idy = yCoor == ly/2;
id  = and(idx,idy);
% Corresponding displacement degree of freedom
id = find(id)*ndof-3;

% Global force vector
Fp = zeros(nDOF,1);
Fp(id) = rho*g*lx*ly*h;        % Same magnitude of the force (self-weight)

%% Distributed (surface) load
% Elements with pressure loading
Elset_p = (1:nx*ny);
% Compute the pressure loading
[Fpres] = Pressure(TopoN,TopoD,xCoor,yCoor,Elset_p,press,flag_func);

%% Boundary conditions
% All degrees of freedom at the boundaries are fixed
% Find bottom boundary edge
idx_be = yCoor == 0;
% Find right boundary edge
idx_re = xCoor == lx;
% Find top boundary edge
idx_te = yCoor == ly;
% Find left boundary edge
idx_le = xCoor == 0;
% Inactive nodes
InactiveNodes = or(idx_be,or(idx_re,or(idx_te,idx_le)));
% Active nodes
ActiveNodes = (1:nNodes)';
ActiveNodes = ActiveNodes(~InactiveNodes);
% Active degrees of freedom
DOF_aux1 = kron(ActiveNodes,ones(1,ndof)*ndof);
DOF_aux2 = kron(ones(size(ActiveNodes)),(ndof-1:-1:0));
ActiveDOFs = unique(DOF_aux1 - DOF_aux2);

%% Solution - Statics
disp('Static solution')
% Gravity load
% uC = K(ActiveDOFs,ActiveDOFs)\Fg(ActiveDOFs);
% Point load
% uC = K(ActiveDOFs,ActiveDOFs)\Fp(ActiveDOFs);
% Pressure (surface) load
uC = K(ActiveDOFs,ActiveDOFs)\Fpres(ActiveDOFs);

% Displacement field
U = zeros(nDOF,1);
U(ActiveDOFs) = uC;

%% Visualization
disp('Visualization of the displacement field')
% Number of subdivisions of each element for the visualization
ndiv = 5;
% Plot the displacement field
PlotResults(TopoN,TopoD,xCoor,yCoor,U,ndiv);

%% Solution - Modal
disp('Eigenvalue problem')
% Number of eigenvalues requested
nEV = 15;
if nEV > numel(ActiveDOFs)
    nEV = numel(ActiveDOFs);
end
% Solve the eigenvalue problem and compute the nEV smallest eigenvalues
[MS,EV] = eigs(K(ActiveDOFs,ActiveDOFs),M(ActiveDOFs,ActiveDOFs),nEV,'sm');

% Determine if the eigenvalues are real
if ~isreal(EV)
    EV = real(EV);
    MS = real(MS);
    % Display warning message
    warning('Eigenvalues and mode shapes are complex!')
end

% Add the fixed degrees of freedom
ModeShapes = zeros(nDOF,nEV);
ModeShapes(ActiveDOFs,:) = MS;
% Eigenfrequencies
EF = sqrt(diag(EV))/(2*pi);

% Plot mode shape
disp('Visualization of a mode shape')
PlotResults(TopoN,TopoD,xCoor,yCoor,ModeShapes(:,nEV),ndiv);

disp(LineDisp)
disp([repmat(Line,1,28) 'End of Simulation!' repmat(Line,1,28)])
disp(LineDisp)