%-----------------------------------------------------------------
% Decomposition of a rectangle in 4 sub-domains
%-----------------------------------------------------------------
%
% Christoph Neuhaus (07.05.2013)
%
%  
% Input:
%     * x, y: Vertices of the original rectangle
%     
%     
% Output:
%     * X, Y: Matrices containing the vertices of the sub-rectangles
%    
%
function [ X,Y ] = DecompositionQuadrilateral( x,y )
%
%
%              _______________
%             |       |       |
%             |   4   |   3   |
%             |_______C_______|
%             |       |       |
%             |   1   |   2   |
%             |_______|_______|
%
%

Cx = x(1) + (x(2) - x(1))/2;
Cy = y(1) + (y(4) - y(1))/2;

%% 1 %%
X(1,1) = x(1);
X(2,1) = Cx;
X(3,1) = Cx;
X(4,1) = x(1);

Y(1,1) = y(1);
Y(2,1) = y(1);
Y(3,1) = Cy;
Y(4,1) = Cy;

%% 2 %%
X(1,2) = Cx;
X(2,2) = x(2);
X(3,2) = x(2);
X(4,2) = Cx;

Y(1,2) = y(1);
Y(2,2) = y(1);
Y(3,2) = Cy;
Y(4,2) = Cy;

%% 3 %%
X(1,3) = Cx;
X(2,3) = x(2);
X(3,3) = x(2);
X(4,3) = Cx;


Y(1,3) = Cy;
Y(2,3) = Cy;
Y(3,3) = y(3); 
Y(4,3) = y(3); 

%% 4 %%
X(1,4) = x(1);
X(2,4) = Cx;
X(3,4) = Cx;
X(4,4) = x(1);

Y(1,4) = Cy;
Y(2,4) = Cy;
Y(3,4) = y(3); 
Y(4,4) = y(3); 

end