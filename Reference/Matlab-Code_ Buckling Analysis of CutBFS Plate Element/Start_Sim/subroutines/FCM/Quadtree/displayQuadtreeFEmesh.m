%% Plot FE-mesh including quadtree decomposition
%-------------------------------------------------------------------------%
%
% (c) 07/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Visualization of the finite element mesh including the quadtree
% decomposition.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:    Nodal topology (connectivity/coincidence) matrix
% CutElem:  Elements that are intersected by the boundary of the
%           computational domain
% xCoor:    Global x-coordinates of the nodes
% yCoor:    Global y-coordinates of the nodes
% k:        Number of hierarchic subdivisions of cut elements
% Geometry: Structure that contains geometric- and meta-information about
% 			the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% FCMvec:       Vector containing a flag whether the points are located
%               within the fictitious domain or the physical one.
% IdxDomain:   	Specifies the domain which has the smallest distance to the
%               considered point
% Hole:         Specifies if the ith domain is a void or a different
%               material (material interface/inclusion)
%
%-------------------------------------------------------------------------%
function [] = displayQuadtreeFEmesh(TopoN,CutElem,xCoor,yCoor,k,Geometry)

figure('Name','FE-mesh including quadtree decomposition','units',...
       'normalized','outerposition',[0 0 1 1])
hold on

% Number of finite elements
nElem = size(TopoN,1);

% Number of void regions/inclusions
nD = numel(CutElem);

% Loop over all elements
for l = 1:nElem
    % Nodal coordinates
    xE = xCoor(TopoN(l,:));
    yE = yCoor(TopoN(l,:));
    
    % Only if the element is intersected by the boundary of the domain the
    % quadtree decomposition has to be executed.
    % Determine if the current element is a cut or standard element
    for j = 1:nD
        ID = find(CutElem{j} == l, 1);
        if ~isempty(ID)
            break;
        end
    end
    
    if isempty(ID)  % FEM
        rectangle('Position',[xE(1),yE(1),xE(3)-xE(1),yE(3)-yE(1)],...
                  'FaceColor',[.5 .5 .5])
    else
        % Quadtree decomposition
        [~,~,domains,~,~,XX,YY] = modQuadtree(xE,yE,k,Geometry,true,true);

        for i = 1:size(XX,2)
            x1 = XX(1,i);
            x2 = XX(2,i);
            x3 = XX(3,i);
            x4 = XX(4,i);
            y1 = YY(1,i);
            y2 = YY(2,i);
            y3 = YY(3,i);
            y4 = YY(4,i);

            if domains(i) == 0       % Physical sub-cell
                rectangle('Position',[x1,y1,x3-x1,y3-y1],...
                          'FaceColor',[.5 .5 .5])
            elseif domains(i) == 1   % Fictitious sub-cell
                rectangle('Position',[x1,y1,x3-x1,y3-y1],...
                          'FaceColor',[.9 .9 .9])
            elseif domains(i) == 2 	% Cut sub-cell
                rectangle('Position',[x1,y1,x3-x1,y3-y1],...
                          'FaceColor',[.7 .7 .7])
            elseif domains(i) == 3 	% Sub-cell was not checked
                rectangle('Position',[x1,y1,x3-x1,y3-y1],...
                          'FaceColor',[255/255 204/255 153/255])
            end

            plot([x1 x2], [y1, y2],'k')
            plot([x2 x3], [y2, y3],'k')
            plot([x3 x4], [y3, y4],'k')
            plot([x4 x1], [y4, y1],'k')
        end
    end
end

% Display the polygonal boundary, if the input variable 'geometry' exists.
if exist('geometry','var')
    for i = 1:size(geometry,1)-1
        plot([geometry(i,1), geometry(i+1,1)],[geometry(i,2),...
            geometry(i+1,2)],'b')
    end
end

xlabel('\xi')
ylabel('\eta')

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