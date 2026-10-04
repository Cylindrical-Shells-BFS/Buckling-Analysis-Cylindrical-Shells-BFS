%% Level-set/Abstandsfunktion für eine Linie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Die Abstandsfunktion für eine Linie wird aus gewertet. Ist der 
% Funktionswert Distance = 0, so befindet sich der Punkt P auf der Linie. 
% Ist der Funktionswert Distance > 0, so befindet sich der Punkt P auf der 
% rechten Seite der Linie. Ist der Funktionswert Distance < 0, so befindet 
% sich der Punkt P auf der linken Seite der Linie.
%
% Per Definition wird ein Punkt, der sich auf der linken Seite der Linie 
% befindet als innerhalb des Gebietes definiert und erhält somit einen
% negativen Funktionewert.
%
% Ein geschlossener Polygonzug muß entgegen der Uhrzeigerrichtung durch-
% laufen werden. Dann kann mittels Intersection die Schnittmenge der
% Gebiete bestimmt werden und somit die Abstandsfuntion für den Polygon-
% zug berechnet werden.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input-Variablen
%
% P (xp,yp):	Koordinaten des zu überprüfenden Punktes				
% P1(x1,y1):	Koordinaten des ersten Eckpunktes der Linie
% P2(x2,y2):	Koordinaten des zweiten Eckpunktes der Linie
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output-Variablen
%
% Distance:		Funktionswert der Abstandsfunktion
% IO:           Im Inneren des physikalischen Gebietes = 0; 
%               Im Inneren des fiktiven Gebietes = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Distance,IO] = DistanceFunctionLine(P,x1,y1,x2,y2)
% Richtungsvektor der Linie
a = [x2-x1,y2-y1]; 
a = a/norm(a);
% Vektor vom Anfangspunkt der Linie zum Punkt P
b = [P(:,1)-x1,P(:,2)-y1];
% x_3-Komponente des Vektorpruduktes; da a ein Einheitsvektor ist, 
% entspricht dies auch dem senkrechten Abstand von P zur Linie
% |a x b| = |a|*|b|*sin(gamma)
Distance = b(:,1)*a(2) - b(:,2)*a(1);

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