%% Koordinaten von P im gedrehten Koordinatensystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Die Funktion berechnet die Koordinaten der Punkte P in einem um den
% Winkel gamma gedrehten Koordinaten System.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input-Variablen
%
% P (xp,yp):	Koordinaten der Punkte
% gamma:        Drehwinkel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output-Variablen
%
% PPrime:		Koordinaten der Punkte im gedrehten Koordinatensystem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PPrime = RotateCoos(P,gamma)

PPrime = ([cos(gamma),sin(gamma);-sin(gamma),cos(gamma)]*P')';

end