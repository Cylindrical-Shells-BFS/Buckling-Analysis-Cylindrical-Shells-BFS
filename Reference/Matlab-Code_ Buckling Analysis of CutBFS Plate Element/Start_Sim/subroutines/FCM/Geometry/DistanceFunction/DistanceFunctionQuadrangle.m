%% Level-set/Abstandsfunktion f¸r ein Polygon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Die Abstandsfunktion f¸r ein Viereck wird aus gewertet. Ist der 
% Funktionswert Distance = 0, so befindet sich der Punkt P auf dem Rand des 
% Vierecks. Ist der Funktionswert Distance > 0, so befindet sich der Punkt
% P auﬂerhalb des Vierecks. Ist der Funktionswert Distance < 0, so befindet 
% sich der Punkt P innerhalb des Vierecks. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input-Variablen
%
% P (xp,yp):	Koordinaten des zu ¸berpr¸fenden Punktes				
% PP(xi,yi):	Koordinaten der Eckpunkte des Vierecks (ccw-orientation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output-Variablen
%
% Distance:		Funktionswert der Abstandsfunktion
% IO:           Im Inneren des physikalischen Gebietes = 0; 
%               Im Inneren des fiktiven Gebietes = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Distance,IO] = DistanceFunctionQuadrangle(P,PP)
% Number of points describing the quadrangle
NumberOfPoints = size(PP,1);

% check if the quadrangle line is closed
if PP(1,1) == PP(NumberOfPoints,1) && PP(1,2) == PP(NumberOfPoints,2)
    NumberOfVertices = NumberOfPoints - 1;
    Nodes = PP;
else
    NumberOfVertices = NumberOfPoints;
    Nodes = [PP; PP(1,1) PP(1,2)];
end
    
dVec = cell(1,NumberOfVertices);

for i = 1:NumberOfVertices
    dVec{i} = DistanceFunctionLine(P,Nodes(i,1),Nodes(i,2),...
              Nodes(i+1,1),Nodes(i+1,2));
end

d3 = dVec{1};
for i = 1:NumberOfVertices-1
    d2 = dVec{i+1};
    d1 = d3;
    d3 = DistanceFunctionBooleanOperationIntersection(d2,d1);
end

Distance = d3;

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