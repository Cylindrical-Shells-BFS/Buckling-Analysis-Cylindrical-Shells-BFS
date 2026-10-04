function [] = displayBooleanQuadtree( BQT, geometry )
% Created by: Marton Petoe
% Date: 13.02.2019

%% Description
% The following code displays the sub-cells resulting from a boolean
% quadtree decomposition (booleanQuadtree.m) of a given cell.

% According to the function "booleanQuadtree.m", the resulting sub-cells 
% can have one of the four possible alphaB values: 1, -1, 11, -11. For the
% definition of these values please see the function "booleanQuadtree.m".
% Following color-coding is used for displaying the sub-cells:

% alphaB =   1  --> green
% alphaB =  -1  --> red
% alphaB =  11  --> dark green
% alphaB = -11  --> dark red

%% Input and output arguments
% INPUT -------------------------------------------------------------------
% BQT: Structure, containing the results of a boolean quadtree
%      decomposition of a single cell. For the definition of the required
%      fields XI, ETA, alphaB and tree, see the output arguments of the
%      function "booleanQuadtree.m".
% geometry: [nV+1 x 2] matrix containing the global coordinates of the nV 
%           vertices of the polygonal chain defining the boundary. The 
%           boundary should be given in a closed form, and so the first and
%           last vertices should be identical. The input of 'geometry' is 
%           optional. If the boundary exists only in an implicit form, 
%           geometry cannot and should not be given. In this case the 
%           following code only displays sub-cells resulting from the 
%           boolean quadtree decomposition.

%% Plot
figure('Name','Boolean quadtree decomposition','units','normalized',...
       'outerposition',[0 0 1 1])
hold on
grid on

for i = 1:size(BQT.XI,2)
    [x,y] = meshgrid([min(BQT.XI(:,i)), max(BQT.XI(:,i))],...
        [min(BQT.ETA(:,i)), max(BQT.ETA(:,i))]);
    
    if BQT.alphaB(i) == 11
        col = [0 .7 0]; % Dark green
    elseif BQT.alphaB(i) == 1
        col = [0 1 0]; % Green
    elseif BQT.alphaB(i) == -11
        col = [.7 0 0]; % Dark red
    elseif BQT.alphaB(i) == -1
        col = [1 0 0]; % Red
    end
    surf(x,y, (nnz(BQT.tree(:,i))-1)*ones(2),'FaceColor',col)
end

% Display the polygonal boundary, if the input variable 'geometry' exists.
if exist('geometry','var')
    for i = 1:size(geometry,1)-1
        plot3([geometry(i,1), geometry(i+1,1)],...
            [geometry(i,2),geometry(i+1,2)], [BQT.k, BQT.k],'k')
    end
end

xlabel('\xi')
ylabel('\eta')
zlabel('k')

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