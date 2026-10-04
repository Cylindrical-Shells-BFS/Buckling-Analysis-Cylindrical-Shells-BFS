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
pa = genpath('..\subroutines');
addpath(pa); 

%% Shape functions
% Local coordinates
xi  = -1:0.05:1;
eta = -1:0.05:1;

% Number of points
nP_1D = numel(xi);

% Grid
[XI,ETA] = meshgrid(xi,eta);
xi_vis  = reshape(XI,1,[]);
eta_vis = reshape(ETA,1,[]);

% Compute shape functions
[N] = SF_BFS(xi_vis,eta_vis);

% Visualization
% Number of shape function to visualize
figure('Name','Shape functions: Bogner-Fox-Schmit plate element',...
       'units','normalized','outerposition',[0 0 1 1])
for i = 1:16
    subplot(4,4,i)
    nSF = i;
    Nvis = reshape(N(nSF,:),nP_1D,nP_1D);
    surf(XI,ETA,Nvis,'FaceLighting','gouraud')
    rotate3d on
end

%% Topology Data - example
% 1-------2-------3
% |       |       |
% |       |       |
% 4-------5-------6
% |       |       |
% |       |       |
% 7-------8-------9
% Sample displacement field
a = 1;
% Element 1
uElem1 = [0,0,0,0,1,a,0,0,0,0,0,0,0,0,0,0]';
% Element 2
uElem2 = [1,a,0,0,1,0,1,0,0,0,0,0,0,0,0,0]';
% Element 3
uElem3 = [0,0,0,0,0,0,0,0,1,a,0,0,0,0,0,0]';
% Element 4
uElem4 = [0,0,0,0,0,0,0,0,1,0,1,0,1,a,0,0]';
% Compute the displacement value for the visualization on a fine grid
Uvis1 = N'*uElem1;
Uvis2 = N'*uElem2;
Uvis3 = N'*uElem3;
Uvis4 = N'*uElem4;
% Visualize the displacement field
figure('Name','4-Elements: Displacement field','units','normalized',...
       'outerposition',[0 0 1 1])
hold on
rotate3d on
% Element 1 (Vetyukov --> green element)
Uvis = reshape(Uvis1,nP_1D,nP_1D);
surf(XI-1,ETA+1,Uvis,'FaceLighting','gouraud')
% Element 2 (Vetyukov --> red element)
Uvis = reshape(Uvis2,nP_1D,nP_1D);
surf(XI+1,ETA+1,Uvis,'FaceLighting','gouraud')
% Element 3 (Vetyukov --> blue element)
Uvis = reshape(Uvis3,nP_1D,nP_1D);
surf(XI-1,ETA-1,Uvis,'FaceLighting','gouraud')
% Element 4 (Vetyukov --> yellow element)
Uvis = reshape(Uvis4,nP_1D,nP_1D);
surf(XI+1,ETA-1,Uvis,'FaceLighting','gouraud')
hold off
axis equal
view(3)
xlim([-2,2])
ylim([-2,2])

%% Material properties and parameters
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
nx = 4;
ny = 4;
% Plate thickness
h = 1e-3;
% Gravity
g = -9.8;
% Plate stiffness parameters
D = (E*h^2)/(12*(1-nu^2))*[1,nu,0;nu,1,0;0,0,(1-nu)/2];

%% Element matrices
% Numerical integration - Gauss-Legendre quadrature rule
% Gauss points and weights - 1d
[xG,wG] = GaussQuad(4,-1,1);
% Gauss points - 2d
[xiMat,etaMat] = meshgrid(xG,xG);
xiG  = reshape(xiMat,[],1);
etaG = reshape(etaMat,[],1);
% Number of integration points
nIP = numel(xiG);
% Weights - 2d
[wxiMat,wetaMat] = meshgrid(wG,wG);
wxiG  = reshape(wxiMat,[],1);
wetaG = reshape(wetaMat,[],1);
wG    = wxiG.*wetaG;

% Mass matrix
% Compute M for one element (identical for all elements - same geometry)
Me = zeros(nSF,nSF);
% Jacobian matrix -- J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
J = [lx/(2*nx), 0; 0, ly/(2*ny)];
% Determinant of the Jacobi matrix
detJ = det(J);
% Loop over all integration points
for i = 1:nIP
    % Compute shape functions
    [N] = SF_BFS(xiG(i),etaG(i));
    % Mass matrix
    Me = Me + N*N'*rho*h*detJ*wG(i);
end

% Stiffness matrix
% Compute K for one element (identical for all elements - same geometry)
Ke = zeros(nSF,nSF);
% Inverse Jacobian matrix (mapping local to global domain)  -- 
% Jinv = [dxi/dx, deta/dx; dxi/dy, deta/dy]
Jinv = [ly/(2*ny), 0; 0, lx/(2*nx)]/det(J);
% Loop over all integration points
for i = 1:nIP
    % Matrix of the second order derivatives of the shape functions w.r.t.
    % the local coordinates
    Bl = [ SF_BFS_dxi2(xiG(i),etaG(i))'   ; ...
           SF_BFS_deta2(xiG(i),etaG(i))'  ; ...
          2*SF_BFS_dxideta(xiG(i),etaG(i))'];
    % Matrix of the second order derivatives of the shape functions w.r.t.
    % the global coordinates
    Bg = [  Bl(1,:)*Jinv(1,1)^2         + ...
          2*Bl(3,:)*Jinv(1,1)*Jinv(1,2) + ...
            Bl(2,:)*Jinv(1,2)^2         ; ...
            Bl(1,:)*Jinv(2,1)^2         + ...
          2*Bl(3,:)*Jinv(2,1)*Jinv(2,2) + ...
            Bl(2,:)*Jinv(2,2)^2         ; ...
            Bl(1,:)*Jinv(1,1)*Jinv(2,1) + ...
            Bl(3,:)*Jinv(2,1)*Jinv(1,2) + ...
            Bl(3,:)*Jinv(1,1)*Jinv(2,2) + ...
            Bl(2,:)*Jinv(1,2)*Jinv(2,2)];
    % Stiffness matrix
    Ke = Ke + Bg'*D*Bg*h*detJ*wG(i);
end

%% Gravity loading
% Compute Fg for one element (identical for all elements - same geometry)
Fe = zeros(nSF,1);
% Loop over all integration points
for i = 1:nIP
    % Compute shape functions
    [N] = SF_BFS(xiG(i),etaG(i));
    % Determinant of the Jacobi matrix
    detJ = lx/(2*nx)*ly/(2*ny);
    % Gravity load
    Fe = Fe + N*rho*g*h*detJ*wG(i);
end
fmax = max(abs(Fe));
Fe(abs(Fe)<fmax*1e-10) = 0;

%% Assembly process
% Number of degrees of freedom per nodes
ndof = 4;
% Number of elements
nElem = nx*ny;
% Nodal topology data (connectivity matrix)
TopoN = zeros(nElem,4);
counter1 = 0;
counter2 = 0;
for i = 1:nElem
    TopoN(i,:) = [1, 2, nx+3, nx+2] + counter1 + counter2;
    counter2 = counter2 + 1;
    if mod(i,nx) == 0
        counter1 = counter1 + nx + 1;
        counter2 = 0;
    end
end
% Number of nodes
nNodes = max(max(TopoN));
% Number of degrees of freedom
nDOF = nNodes*ndof;
% Nodal coordinates
xCoor = zeros(nNodes,1);
yCoor = zeros(nNodes,1);
for i = 1:ny+1
    xCoor(1+(nx+1)*(i-1):(nx+1)*i) = linspace(0,lx,nx+1);
    yCoor(1+(nx+1)*(i-1):(nx+1)*i) = ly/ny*(i-1);
end
% DOF topology data (connectivity matrix)
Topo_aux1 = kron(TopoN,ones(1,ndof)*ndof);
Topo_aux2 = kron(ones(size(TopoN)),(ndof-1:-1:0));
TopoD = Topo_aux1 - Topo_aux2;

% Global stiffness matrix
K = zeros(nDOF,nDOF);
for i = 1:nElem
    K(TopoD(i,:),TopoD(i,:)) = K(TopoD(i,:),TopoD(i,:)) + Ke;
end
kmax = max(max(abs(K)));
K(abs(K)<kmax*1e-10) = 0;

% Global mass matrix
M = zeros(nDOF,nDOF);
for i = 1:nElem
    M(TopoD(i,:),TopoD(i,:)) = M(TopoD(i,:),TopoD(i,:)) + Me;
end
mmax = max(max(abs(M)));
M(abs(M)<mmax*1e-10) = 0;

% Global force vector
F1 = zeros(nDOF,1);
for i = 1:nElem
    F1(TopoD(i,:),1) = F1(TopoD(i,:),1) + Fe;
end
fmax = max(abs(F1));
F1(abs(F1)<fmax*1e-10) = 0;

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

%% Solution
% Static system
uC1 = K(ActiveDOFs,ActiveDOFs)\F1(ActiveDOFs);

% Displacement field
U1 = zeros(nDOF,1);
U1(ActiveDOFs) = uC1;

%% Visualization
% Visualize the displacement field
figure('Name','Clamped plate: Self-weight','units','normalized',...
       'outerposition',[0 0 1 1])
hold on
rotate3d on

% Local coordinates
xi  = -1:0.05:1;
eta = -1:0.05:1;

% Number of points
nP_1D = numel(xi);

% Grid
[XI,ETA] = meshgrid(xi,eta);
xi_vis  = reshape(XI,1,[]);
eta_vis = reshape(ETA,1,[]);

% Compute shape functions
[N] = SF_BFS(xi_vis,eta_vis);
[Ngeo] = SF_4nodeElem(xi_vis,eta_vis);

for i = 1:nElem
    % Element solution
    uElem = U1(TopoD(i,:));
    % Displacement field
    Uvis = N'*uElem;
    % Nodal coordinates (global)
    xElem = xCoor(TopoN(i,:));
    yElem = yCoor(TopoN(i,:));
    xg = Ngeo'*xElem;
    yg = Ngeo'*yElem;
    % Element i
    Uvis = reshape(Uvis,nP_1D,nP_1D);
    XMat = reshape(xg,nP_1D,nP_1D);
    YMat = reshape(yg,nP_1D,nP_1D);
    surf(XMat,YMat,Uvis,'FaceLighting','gouraud')
end
grid on
hold off
view(3)

%% Example 2: Concentrated force in the middle
% Find node at (lx/2,ly/2)
idx = xCoor == lx/2;
idy = yCoor == ly/2;
id  = and(idx,idy);
% Corresponding displacement degree of freedom
id = find(id)*ndof-3;

% Global force vector
F2 = zeros(nDOF,1);
F2(id) = rho*g*lx*ly*h;        % Same magnitude of the force (self-weight)

%% Solution
% Static system
uC2 = K(ActiveDOFs,ActiveDOFs)\F2(ActiveDOFs);

% Displacement field
U2 = zeros(nDOF,1);
U2(ActiveDOFs) = uC2;

%% Visualization
% Visualize the displacement field
figure('Name','Clamped plate: Point force','units','normalized',...
       'outerposition',[0 0 1 1])
hold on
rotate3d on

% Local coordinates
xi  = -1:0.05:1;
eta = -1:0.05:1;

% Number of points
nP_1D = numel(xi);

% Grid
[XI,ETA] = meshgrid(xi,eta);
xi_vis  = reshape(XI,1,[]);
eta_vis = reshape(ETA,1,[]);

% Compute shape functions
[N] = SF_BFS(xi_vis,eta_vis);
[Ngeo] = SF_4nodeElem(xi_vis,eta_vis);

for i = 1:nElem
    % Element solution
    uElem = U2(TopoD(i,:));
    % Displacement field
    Uvis = N'*uElem;
    % Nodal coordinates (global)
    xElem = xCoor(TopoN(i,:));
    yElem = yCoor(TopoN(i,:));
    xg = Ngeo'*xElem;
    yg = Ngeo'*yElem;
    % Element i
    Uvis = reshape(Uvis,nP_1D,nP_1D);
    XMat = reshape(xg,nP_1D,nP_1D);
    YMat = reshape(yg,nP_1D,nP_1D);
    surf(XMat,YMat,Uvis,'FaceLighting','gouraud')
end
grid on
hold off
view(3)

%% Example 3: Modal analysis
% Number of eigenvalues requested
nEV = 9;
% Solve the eigenvalue problem and compute the nEV smallest eigenvalues
[MS,EV] = eigs(K(ActiveDOFs,ActiveDOFs),M(ActiveDOFs,ActiveDOFs),nEV,'sm');
% Add the fixed degrees of freedom
ModeShapes = zeros(nDOF,nEV);
ModeShapes(ActiveDOFs,:) = MS;
% Eigenfrequencies
EF = sqrt(diag(EV))/(2*pi);

% Visualize the displacement field
figure('Name','Clamped plate: Eigenmodes','units','normalized',...
       'outerposition',[0 0 1 1])

% Local coordinates
xi  = -1:0.1:1;
eta = -1:0.1:1;

% Number of points
nP_1D = numel(xi);

% Grid
[XI,ETA] = meshgrid(xi,eta);
xi_vis  = reshape(XI,1,[]);
eta_vis = reshape(ETA,1,[]);

% Compute shape functions
[N] = SF_BFS(xi_vis,eta_vis);
[Ngeo] = SF_4nodeElem(xi_vis,eta_vis);

% Subplot layout
nSP = ceil(sqrt(nEV));
for j = 1:nEV
    subplot(nSP,nSP,j)
    hold on
    rotate3d on
    for i = 1:nElem
        % Element solution
        uElem = ModeShapes(TopoD(i,:),j);
        % Displacement field
        Uvis = N'*uElem;
        % Nodal coordinates (global)
        xElem = xCoor(TopoN(i,:));
        yElem = yCoor(TopoN(i,:));
        xg = Ngeo'*xElem;
        yg = Ngeo'*yElem;
        % Element i
        Uvis = reshape(Uvis,nP_1D,nP_1D);
        XMat = reshape(xg,nP_1D,nP_1D);
        YMat = reshape(yg,nP_1D,nP_1D);
        surf(XMat,YMat,Uvis,'FaceLighting','gouraud')
    end
    grid on
    hold off
    view(3)
end