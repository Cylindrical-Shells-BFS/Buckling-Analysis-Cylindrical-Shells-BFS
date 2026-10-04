%% Visualization - FE-mesh
%-------------------------------------------------------------------------%
%
% (c) 03/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Plot the elements of the FE-model.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:    Nodal topology (connectivity/coincidence) matrix
% xCoor:    Global x-coordinates of the nodes
% yCoor:    Global y-coordinates of the nodes
%
%-------------------------------------------------------------------------%
function FE_Grid(xCoor,yCoor,TopoN)

figure('Name','Finite Element Mesh','units','normalized',...
       'outerposition',[0 0 1 1])
hold on

% Number of elements
nElem = size(TopoN,1);

for count = 1:nElem
    % Nodal coordinates of the finite element
    x = xCoor(TopoN(count,:));
    y = yCoor(TopoN(count,:));
 
    % Plot patch
    faccolor = [210/255 210/255 210/255];
    patch('Faces',[1 2 3 4],'Vertices',[x,y],...
    'FaceColor',faccolor,'EdgeColor','k','LineWidth',1)
end

set(gca,'FontName','Times','FontSize',24);
set(gcf,'PaperPositionMode','auto')
axis equal
hold off

% axis off
% h = gcf;
% iptsetpref('ImshowBorder','tight')
% print(h,'-depsc','-painters','Mesh.eps')
% 
% savefig('Mesh.fig')

end