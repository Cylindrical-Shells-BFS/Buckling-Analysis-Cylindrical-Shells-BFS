%% Boolesche Operation: Schnittmenge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Die Abstandsfunktion für die Schnittmenge zweier Gebiete wird
% ausgewertet. Ist der Funktionswert Distance = 0, so befindet sich der
% Punkt P auf dem Rand des entstanden Gebietes. Ist der Funktionswert
% Distance > 0, so befindet sich der Punkt P außerhalb des entstanden
% Gebietes. Ist der Funktionswert Distance < 0, so befindet sich der Punkt
% P innerhalb des entstanden Gebietes.
% 
% Das entstandene Gebiet ist die Schnittmenge von Gebiet Omega_1 
% (repräsentiert durch den Funktionswert der Abstandfunktions d1) und 
% Gebiet Omega_2 (repräsentiert durch den Funktionswert der 
% Abstandfunktions d2). D.h. der Punkt liegt innerhalb von Omega_1 und 
% gleichzeitig auch innerhalb von Omega_2. Dadurch muss gelten:
%
% d1 < 0 und d2 < 0.
%
% Das maximum beider Absatndwerte muss also immer < 0 sein. Sonst liegt 
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
function [Distance,Idx] = ...
          DistanceFunctionBooleanOperationIntersection(d1,d2)

[Distance,Idx] = max([d1,d2],[],2);

end