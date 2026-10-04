%% FCM - check geometry
%-------------------------------------------------------------------------%
%
% (c) 16.08.2016, Sascha Duczek, OvGU Magdeburg, Germany
%
% Determine whether a point is located in the fictitious domain --> FCM.
%
%-------------------------------------------------------------------------%
%% Input
%
% [x,y]:            Coordinates of points that are tested
% Geometry:         Structure that contains geometric- and meta-
%                   information about the given boundary.
%
%-------------------------------------------------------------------------%
%% Output
%
% FCMvec:       Vector containing a flag whether the points are located
%               within the geometric domains defined by the variable
%               Geometry. A value of 1 means that the points (x,y) are
%               located in one of the domains and accordingly 0 signifies
%               that they are not in any of those domains.
% IdxDomain:   	Specifies the domain which has the smallest distance to the
%               considered point
% Hole:         Specifies if the ith domain is an inclusion (0), void (1),
%               or the physical domain (2); physical domain means that the
%               outer boundary is described by the geometric entity.
%
%-------------------------------------------------------------------------%
function [FCMvec,IdxDomain,Hole] = PointInMembershipTest(x,y,Geometry)

% Level set function
[LevelSetValue,IdxDomain,Hole] = LevelSetFunctionGeometry(x,y,Geometry);

% Division by zero is not defined
LevelSetValue(LevelSetValue == 0) = -eps;

% Components are only allowed to have the values 1 or 0, i.e. either they
% are inside (1) or outside of the geometric domain(s) defined in the 
% variable Geometry
LevelSetValue = -LevelSetValue./abs(LevelSetValue);
LevelSetValue(LevelSetValue < 1) = 0;

% Status variable: Point is within the fictitious domain (1) or within the
% physical domain (0)
FCMvec = LevelSetValue;

end