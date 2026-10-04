%% Level-set/Abstandsfunktion f¸r ein Rechteck
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Die Abstandsfunktion f¸r ein Rechteck wird aus gewertet. Ist der 
% Funktionswert Distance = 0, so befindet sich der Punkt P auf dem Rand des 
% Rechtecks. Ist der Funktionswert Distance > 0, so befindet sich der Punkt
% P auﬂerhalb des Rechtecks. Ist der Funktionswert Distance < 0, so befindet 
% sich der Punkt P innerhalb des Rechtecks. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input-Variablen
%
% P (xp,yp):	Koordinaten des zu ¸berpr¸fenden Punktes				
% P1(x1,y1):	Koordinaten des ersten Eckpunktes des Rechteck
% P2(x2,y2):	Koordinaten des zweiten Eckpunktes des Rechteck
% P1 und P2 bilden eine Diagonale des Rechtecks.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output-Variablen
%
% Distance:		Funktionswert der Abstandsfunktion
% IO:           Im Inneren des physikalischen Gebietes = 0; 
%               Im Inneren des fiktiven Gebietes = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Distance,IO] = DistanceFunctionRectangle(P,x1,x2,y1,y2)
d = [x1-P(:,1), P(:,1)-x2, y1-P(:,2), P(:,2)-y2];
Distance = max(d,[],2);

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