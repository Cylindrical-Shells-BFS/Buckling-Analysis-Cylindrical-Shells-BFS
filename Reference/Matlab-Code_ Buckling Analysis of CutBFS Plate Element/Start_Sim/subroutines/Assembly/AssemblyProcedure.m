%% Assembly of the system matrices and load vector
%-------------------------------------------------------------------------%
%
% (c) 03/05/2019, Sascha Duczek, UNSW Sydney, Australia
%
% The overall (global) stiffness, mass matrices and load vector are
% assembled.
%
%-------------------------------------------------------------------------%
%% Input
%
% TopoN:    Nodal topology (connectivity/coincidence) matrix
% TopoD:    DOF topology (connectivity/coincidence) matrix
% xCoor:    Global x-coordinates of the nodes
% yCoor:    Global y-coordinates of the nodes
% h:       	Plate thickness
% D:      	Plate elasticity matrix
% rho:    	Mass density
% g:    	Acceleration due to gravity
% flagE:    Flag if elements are identical
%
%-------------------------------------------------------------------------%
%% Output
%
% K:        Global stiffness matrix
% M:        Global mass matrix
% F:        Global load vector (due to gravity load)
%
%-------------------------------------------------------------------------%
function [K,M,F] = AssemblyProcedure(TopoN,TopoD,xCoor,yCoor,h,D,rho,g,...
                   flagE)

% Start time measurement
tic;

% Number of elements
nElem = size(TopoD,1);

% Number of degrees of freedom
nDOF = max(max(TopoD));

% Save results in cell arrays
Kcell     = cell(nElem,1);
Mcell     = cell(nElem,1);
Fcell     = cell(nElem,1);
topo_row  = cell(nElem,1);
topo_col  = cell(nElem,1);

% Initialize the waitbar
anz = waitbar(0,'Computing FE element matrices');

% Compute the element matrices
for i = 1:nElem
    % Nodal coordinates
    xElem = xCoor(TopoN(i,:));
    yElem = yCoor(TopoN(i,:));
    
    % Elemental stiffness matrix
    [Ke] = StiffMat(xElem,yElem,h,D);
    
    % Elemental mass matrix
    [Me] = MassMat(xElem,yElem,h,rho);
    
    % Elemental force vector due to gravity
    [Fe] = Gravity(xElem,yElem,h,rho,g);
    
    % If the element are identical leave the loop
    if flagE == 1
        for j = 1:nElem
            % Save elemental matrices
            Kcell{j} = reshape(Ke,[],1);
            Mcell{j} = reshape(Me,[],1);
            Fcell{j} = Fe;
            topo_row{j}  = reshape(repmat(TopoD(j,:),16,1),[],1);
            topo_col{j}  = reshape(repmat(TopoD(j,:)',16,1),[],1);
			
			% Waitbar
			waitbar(j/nElem,anz,['Computing FE element matrices: ' ...
			num2str(j/nElem*100,'%5.1f') '%'])
        end
        break
    else
        % Save elemental matrices
        Kcell{i}    = reshape(Ke,[],1);
        Mcell{i}    = reshape(Me,[],1);
        Fcell{i}    = Fe;
        topo_row{i} = reshape(repmat(TopoD(i,:),16,1),[],1);
        topo_col{i} = reshape(repmat(TopoD(i,:)',16,1),[],1);
    end
	
	% Waitbar
	waitbar(i/nElem,anz,['Computing FE element matrices: ' ...
    num2str(i/nElem*100,'%5.1f') '%'])
end

% Concatenate all the cell arrays into a matrix
k   = cell2mat(Kcell);
m   = cell2mat(Mcell);
row = cell2mat(topo_row);
col = cell2mat(topo_col);

% Save the result vectors in sparse matrices
% Global stiffness matrix
K = sparse(row,col,k,nDOF,nDOF);
% Delete small components
kmax = max(max(abs(K)));
K(abs(K)<kmax*1e-10) = 0;

% Global mass matrix
M = sparse(row,col,m,nDOF,nDOF); 
% Delete small components
mmax = max(max(abs(M)));
M(abs(M)<mmax*1e-10) = 0;

% Global force vector
F = zeros(nDOF,1);
for i = 1:nElem
    F(TopoD(i,:),1) = F(TopoD(i,:),1) + Fcell{i};
end
fmax = max(abs(F));
F(abs(F)<fmax*1e-10) = 0;

% Close the waitbar
close(anz)

% Stop time measurement
toc

end