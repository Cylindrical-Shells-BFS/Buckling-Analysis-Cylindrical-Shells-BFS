function [] = displayQuadtree( QT, geometry, treeLabels)
% Created by: Marton Petoe
% Date: 10.01.2019

%% Description
% The following code displays the individual sub-cells resulting from the
% quadtree decomposition (modQuadtree.m) of a single cell. 

% Once the quadtree decomposition has been performed, the resulting sub-
% cells are classified to either physical, fictitious or cut. Additionally,
% sub-cells can also have a 4th kind of classification, if they were not
% checked at all. For further infos about the not checked cells please see
% the section "Checking the leaf-cells" in the function "modQuadtree.m".
% Following color-coding is used:
%
% Light gray: Fictitious sub-cell
% Mid gray: Cut sub-cell
% Dark gray: Physical sub-cell
% Salmon color: Sub-cell was not checked

%% Input and output arguments
% INPUT -------------------------------------------------------------------
% QT: Structure, containing the results of the quadtree decomposition of a 
%     single cell. For the definition of the required fields XI, ETA, tree
%     and domains, see the output arguments of the function "modQuadtree.m".
% geometry: [nV+1 x 2] matrix containing the global coordinates of the nV 
%           vertices of the polygonal chain defining the boundary. The 
%           boundary should be given in a closed form, and so the first and
%           last vertices should be identical. The input of 'geometry' is 
%           optional. If the boundary exists only in an implicit form, 
%           geometry cannot and should not be given. In this case the 
%           following code only displays sub-cells resulting from the 
%           quadtree decomposition.

%% Plot

figure('Name','Quadtree decomposition','units','normalized',...
       'outerposition',[0 0 1 1])
hold on

if ~exist('treeLabels','var')
    treeLabels = false;    
end

for i = 1:size(QT.XI,2)
    x1 = QT.XI(1,i);
    x2 = QT.XI(2,i);
    x3 = QT.XI(3,i);
    x4 = QT.XI(4,i);
    y1 = QT.ETA(1,i);
    y2 = QT.ETA(2,i);
    y3 = QT.ETA(3,i);
    y4 = QT.ETA(4,i);
    
    if QT.domains(i) == 0       % Physical sub-cell
        rectangle('Position',[x1,y1,x3-x1,y3-y1],'FaceColor',[.5 .5 .5])
    elseif QT.domains(i) == 1   % Fictitious sub-cell
        rectangle('Position',[x1,y1,x3-x1,y3-y1],'FaceColor',[.9 .9 .9])
    elseif QT.domains(i) == 2 	% Cut sub-cell
        rectangle('Position',[x1,y1,x3-x1,y3-y1],'FaceColor',[.7 .7 .7])
    elseif QT.domains(i) == 3 	% Sub-cell was not checked
        rectangle('Position',[x1,y1,x3-x1,y3-y1],'FaceColor',...
                  [255/255 204/255 153/255])
    end
    
    plot([x1 x2], [y1, y2],'k')
    plot([x2 x3], [y2, y3],'k')
    plot([x3 x4], [y3, y4],'k')
    plot([x4 x1], [y4, y1],'k')
    
    if treeLabels
        str = '';
        for j = 1:size(QT.tree,1)
            if QT.tree(j,i) ~= 0
                str = [str num2str(QT.tree(j,i))];
            else
                break
            end
        end
        text((x1+x3)/2-0.06,(y1+y3)/2,str,'FontSize',10-QT.levels(i));
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