%% Jacobian matrix
%-------------------------------------------------------------------------%
%
% (c) 03/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% Compute the Jacobian matrix for a 4-node element (geometry mapping).
%
%-------------------------------------------------------------------------%
%% Input
%
% (xi,eta):     Local coordinates of the evaluation points (have to be
%               provided as vectors)
% (x,y):        Global coordinates of the nodes
%
%-------------------------------------------------------------------------%
%% Output
%
% J:            Matrix of shape functions evaluated at (xi,eta)
%
%-------------------------------------------------------------------------%
function [J] = JacobianMatrix(xi,eta,x,y)

% Jacobian matrix (2D)
J = zeros(2,2);

% Compute components
J_11 = dot(SF_4nodeElem_dxi(eta),x);
J_12 = dot(SF_4nodeElem_dxi(eta),y);
J_21 = dot(SF_4nodeElem_deta(xi),x);
J_22 = dot(SF_4nodeElem_deta(xi),y);

% Assign values
J(1,1) = J_11;
J(1,2) = J_12;
J(2,1) = J_21;
J(2,2) = J_22;

end