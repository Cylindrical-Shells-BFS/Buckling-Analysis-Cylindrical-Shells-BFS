%% Amplitude function
%-------------------------------------------------------------------------%
%
% (c) 06/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the amplitude a specific locations. Can be used for the
% application of spatially varying loads and Dirichlet boundary conditions.
%
%-------------------------------------------------------------------------%
%% Input
%
% (x,y):	Global coordinates of the evaluation point
%
%-------------------------------------------------------------------------%
%% Output
%
% Amp:      Amplitude of the function
%
%-------------------------------------------------------------------------%
function [Amp] = AmplitudeFunc(x,y)

% Function
Amp = sin(pi*x).*sin(pi*y);

end