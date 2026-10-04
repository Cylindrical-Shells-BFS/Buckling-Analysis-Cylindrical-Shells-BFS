%% Diagonal preconditioning
%-------------------------------------------------------------------------%
%
% (c) 11/04/2022, Sascha Eisentraeger, UNSW Sydney, Australia
%
% Diagonal preconditioning procedure for linear systems of equations.
%
% Original system of equations:
% A*x = b
%
% Preconditioned system of equations:
% D*A*D'*y = D*b, with x = D'*y.
% To obtain the solution to the original system of equations x, an
% additional operation (matrix multiplication) is required: x = D'*y.
%
%-------------------------------------------------------------------------%
%% Input
%
% A:	System matrix
% b:	Right-hand-side vector
%
%-------------------------------------------------------------------------%
%% Output
%
% Acond:	Preconditioned system matrix
% bcond:	Preconditioned right-hand-side vector
% D			Diagonal preconditioning matrix
%
%-------------------------------------------------------------------------%
function [Acond,bcond,D] = PreconditioningLSE(A,b)

% Number of rows/columns
nrc = size(A,1);

% Generate the diagonal preconditioning matrix
D = spdiags(1./sqrt(diag(A)),0,nrc,nrc);

% Compute preconditioned system matrix (transformation is symmetry
% preserving)
Acond = D*A*(D');
bcond = D*b;

end