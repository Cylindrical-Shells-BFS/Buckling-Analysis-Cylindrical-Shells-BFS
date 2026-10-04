%% Level-set/Abstandsfunktion für ein Polygon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Die Abstandsfunktion fuer ein Polygon wird aus gewertet. Ist der 
% Funktionswert Distance = 0, so befindet sich der Punkt P auf dem Polygon. 
% Ist der Funktionswert Distance > 0, so befindet sich der Punkt P 
% auszerhalb des Polygons. Ist der Funktionswert Distance < 0, so befindet 
% sich der Punkt P innerhalb des Polygons.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input-Variablen
%
% xp:			vector of points X coordinates (1 X np)
% yp:			vector of points Y coordinates (1 X np)
% xv:			vector of polygon vertices' X coordinates (1 X nv)
% yv:			vector of polygon vertices' Y coordinates (1 X nv)
% flag:			p_poly_dist.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output-Variablen
%
% Distance:		Funktionswert der Abstandsfunktion
% IO:           Im Inneren des physikalischen Gebietes = 0; 
%               Im Inneren des fiktiven Gebietes = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Distance,IO] = DistanceFunctionPolygon(xp,yp,xv,yv,flag)

% Distance to the boundary of the polygon
Distance = p_poly_dist(xp,yp,xv,yv,flag);

% Round the computed distance values
Distance = round(Distance,12);

Dist = Distance;
Dist(Dist <= 0) = 0;
Dist(Dist > 0)  = 1;

% If three values are the same, check if the forth distance is 0. Depending
% on the value of Dist the element is either a standard finite element or a
% fictitious element.
if sum(Dist) == 1
	id = find(Dist == 1);
    if Distance(id) == 1
        Dist(id) = 0;
        Distance(id) = -10*eps;
    end
elseif sum(Dist) == 3
    id = find(Dist == 0);
    if Distance(id) == 0
        Dist(id) = 1;
        Distance(id) = 10*eps;
    end
end

IO = Dist;
end