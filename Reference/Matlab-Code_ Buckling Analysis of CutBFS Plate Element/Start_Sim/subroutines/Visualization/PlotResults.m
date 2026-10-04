%% Visualization
%-------------------------------------------------------------------------%
%
% (c) 03/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Plot the results of the analysis.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:    Nodal topology (connectivity/coincidence) matrix
% TopoD:    DOF topology (connectivity/coincidence) matrix
% xCoor:    Global x-coordinates of the nodes
% yCoor:    Global y-coordinates of the nodes
% U:        Results vector
% ndiv:     Number of subdivisions for the visualization
%
%-------------------------------------------------------------------------%
function [] = PlotResults(TopoN,TopoD,xCoor,yCoor,U,ndiv)

% Visualize the displacement field
figure('Name','Bogner-Fox-Schmitt plate element: Results','units',...
       'normalized','outerposition',[0 0 1 1])
hold on
rotate3d on

% Number of elements
nElem = size(TopoD,1);

if ndiv <= 1
    ndiv = 1;
end

% Local coordinates
xi  = -1:2/ndiv:1;
eta = -1:2/ndiv:1;

% Number of points
nP_1D = numel(xi);

% Grid
[XI,ETA] = meshgrid(xi,eta);
xi_vis  = reshape(XI,1,[]);
eta_vis = reshape(ETA,1,[]);

% Compute shape functions
[N] = SF_BFS(xi_vis,eta_vis);
[Ngeo] = SF_4nodeElem(xi_vis,eta_vis);

XMatC = cell(nElem,1);
YMatC = cell(nElem,1);
UvisC = cell(nElem,1);
for i = 1:nElem
    % Element solution
    uElem = U(TopoD(i,:));
    % Displacement field
    Uvis = N'*uElem;
    % Nodal coordinates (global)
    xElem = xCoor(TopoN(i,:));
    yElem = yCoor(TopoN(i,:));
    xg = Ngeo'*xElem;
    yg = Ngeo'*yElem;
    % Element i
    UvisC{i} = [reshape(Uvis,nP_1D,nP_1D); NaN(1,nP_1D)];
    XMatC{i} = [reshape(xg,nP_1D,nP_1D); NaN(1,nP_1D)];
    YMatC{i} = [reshape(yg,nP_1D,nP_1D); NaN(1,nP_1D)];
end
XMat = cell2mat(XMatC);
YMat = cell2mat(YMatC);
UMat = cell2mat(UvisC);
surf(XMat,YMat,UMat,'FaceLighting','gouraud')
grid on
hold off
view(3)

end