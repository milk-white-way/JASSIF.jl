
 function x = precLU(PRE, rhs) 
%% function rhs = precLU(PRE, rhs) 
%% arms preconditioning operation
%% PRE = struct for preconditioner
%%-------------------------------------------------
fprintf('Apply ILU preconditioner...\n');
L = PRE.L;
U = PRE.U; 
x = U \ (L \ rhs);

