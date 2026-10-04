%% Level-set/Abstandsfunktion f¸r eine Ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Die Abstandsfunktion f¸r eine Ellipse wird aus gewertet. Ist der 
% Funktionswert Distance = 0, so befindet sich der Punkt P auf der Ellipse. 
% Ist der Funktionswert Distance > 0, so befindet sich der Punkt P 
% auﬂerhalb der Ellipse. Ist der Funktionswert Distance < 0, so befindet 
% sich der Punkt P innerhalb der Ellipse.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input-Variablen
%
% P (xp,yp):	Koordinaten des zu ¸berpr¸fenden Punktes				
% Pc(xc,yc):	Koordinaten des Ellipsenmittelpunktes
% a:            Erste Halbachse
% b:			Zweite Halbachse
% gamma:        Winkellage der Ellipse (Angabe in radiant; Drehung im
%               Gegenuhrzeigersinn-mathematisch positive Drehrichtung)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output-Variablen
%
% Distance:		Funktionswert der Abstandsfunktion
% IO:           Im Inneren des physikalischen Gebietes = 0; 
%               Im Inneren des fiktiven Gebietes = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Distance,IO] = DistanceFunctionEllipse(P,xc,yc,a,b,gamma)
% Compute distance from points P to the ellipse centered at the origin with 
% axes=[a,b]. If the ellipse is not centered at the origin the points P
% need to be translated --> P_new.
P_new = [P(:,1)-xc, P(:,2)-yc];
PPrime = RotateCoos(P_new,gamma);
Distance = dellipse(PPrime,[a,b]);

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