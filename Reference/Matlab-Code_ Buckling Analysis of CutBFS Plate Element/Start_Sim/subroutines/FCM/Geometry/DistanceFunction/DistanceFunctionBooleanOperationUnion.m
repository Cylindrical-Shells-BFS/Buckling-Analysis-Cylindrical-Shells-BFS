%% Boolesche Operation: Vereinigung
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Die Abstandsfunktion für die Vereinigung zweier Gebiete wird ausgewertet.
% Ist der Funktionswert Distance = 0, so befindet sich der Punkt P auf dem 
% Rand des entstanden Gebietes. Ist der Funktionswert Distance > 0, so 
% befindet sich der Punkt P außerhalb des entstanden Gebietes. Ist der 
% Funktionswert Distance < 0, so befindet sich der Punkt P innerhalb des 
% entstanden Gebietes.
% 
% Das entstandene Gebiet ist die Vereinigung von Gebiet Omega_1 
% (repräsentiert durch den Funktionswert der Abstandfunktions d1) und 
% Gebiet Omega_2 (repräsentiert durch den Funktionswert der 
% Abstandfunktions d2). D.h. der Punkt muss mindestens innerhalb von 
% Omega_1 oder Omega_2 liegen. Dadurch muss gelten:
% 1) d1 < 0 und d2 > 0 oder
% 2) d1 > 0 und d2 < 0 oder
% 3) d1 < 0 und d2 < 0.
% Es reicht also aus, wenn mindestens ein Abstandswert < 0. Sonst liegt 
% der Punkt P nicht im Gebiet, welches durch die Boolesche Operation 
% entstanden ist.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input-Variablen
%
% d1:		Funktionswert der Absatndfunktion von Gebiet Omega_1
% d2:		Funktionswert der Absatndfunktion von Gebiet Omega_2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output-Variablen
%
% Distance:		Funktionswert der Abstandsfunktion
% Idx:          Der kleinere Abstandswert ist entweder in Vektor d1 (=1) 
%               oder in d2 (=2) zu finden
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Distance,Idx] = DistanceFunctionBooleanOperationUnion(d1,d2)

[Distance,Idx] = min([d1,d2],[],2);

end