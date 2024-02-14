function J=Jacobian_Assemble(x)

J_t = J_Time_Derivative(x);
J_D = J_Viscous_Tri_Diag(x);
J_C = J_Convection_Tri_Diag(x);
J = J_t + J_D + J_C;
