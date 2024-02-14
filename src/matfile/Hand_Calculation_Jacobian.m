function J = Hand_Calculation_Jacobian(x)
M = length(x.Ucont_x(:,1));
N = length(x.Ucont_x(1,:));

J = spalloc(2*M*N,2*M*N, 22 * M*N);

Ucont_x_0 = x.Ucont_x;
Ucont_y_0 = x.Ucont_y;
dx      = x.dx;
dy      = x.dy;
dt      = x.dt;
Ucat_x_0  = x.Ucat_x;
Ucat_y_0  = x.Ucat_y;

M4 = M+4;
N4 = N+4;

Ucat_x(1:M4,1:N4) = 0;
Ucat_y(1:M4,1:N4) = 0;
Ucont_x(1:M4,1:N4) = 0;
Ucont_y(1:M4,1:N4) = 0;
%Artificial Ucat_x, Ucat_y, Ucont_x Ucont_y
counter = 0;
for i = 1:M4
    for j = 1:N4
        
        if (i >= 5 && j >= 5)
        Ucat_x(i,j) = Ucat_x_0(i-4,j-4);   
        Ucat_y(i,j) = Ucat_y_0(i-4,j-4);
        Ucont_x(i,j) = Ucont_x_0(i-4,j-4);
        Ucont_y(i,j) = Ucont_y_0(i-4,j-4);   
        end
    end
end

%---------------  First M*N values of F ----------------------------------
for i = 1:M4
    for j =1:N4

             i0 = i-4;
             j0 = j-4;

        
         if (i >=6  && i<=M4-5 && j>=6 && j<=N4-5)

             
        % Set up all the indexes
        index_N_N_U = glidx(i0,j0-2,M,N); %1
        index_N_U   = glidx(i0,j0-1,M,N); %2
        index_P_U   = glidx(i0,j0,M,N)  ; %3
        index_S_U   = glidx(i0,j0+1,M,N); %4
        index_S_S_U = glidx(i0,j0+2,M,N); %5
        
        index_W_U   = glidx(i0+1,j0,M,N); %6
        index_W_W_U = glidx(i0+2,j0,M,N); %7       
        index_E_U   = glidx(i0-1,j0,M,N); %8
        index_E_E_U = glidx(i0-2,j0,M,N); %9 
        
     
        index_N_N_V = glidx(i0,j0-2,M,N) + M*N; %10
        index_N_V   = glidx(i0,j0-1,M,N) + M*N; %11
        index_P_V   = glidx(i0,j0,M,N)   + M*N; %12
        index_S_V   = glidx(i0,j0+1,M,N) + M*N; %13
        index_S_S_V = glidx(i0,j0+2,M,N) + M*N; %14
        
        index_W_V   = glidx(i0+1,j0,M,N) + M*N; %15
        index_W_W_V = glidx(i0+2,j0,M,N) + M*N; %16       
        index_E_V   = glidx(i0-1,j0,M,N) + M*N; %17
        index_E_E_V = glidx(i0-2,j0,M,N) + M*N; %18   
        
       

         for k=1:22
            ii(k) = 0;
            jj(k) = 0;
            vals(k)  = 0;
         end
      
         %% ------------ Diagonal \partial F \partial u_i
        
        % Time derivative
        ii(1) = index_P_U;
        jj(1) = index_P_U;
        vals(1) = -1.5/dt;
        
        % Viscous_Flux
        vals(1) = vals(1) - 2/(dx).^2 - 2/(dy).^2;

      
        % \partial F partial u_p
        %Convective flux        
        Fpx1_i     =  Eval_Fp(i-4,M, Ucont_x(i,j),Ucat_x(i+2,j), Ucat_x(i+1,j),Ucat_x(i,j),Ucat_x(i-1,j),1,1);
        Fpx1_i_m_1 =  Eval_Fp(i-1-4,M, Ucont_x(i-1,j),Ucat_x(i+1,j), Ucat_x(i,j),Ucat_x(i-1,j),Ucat_x(i-2,j),1,1);
      
        Fpx2_j     =  Eval_Fp(j-4,N, Ucont_y(i,j),Ucat_x(i,j+2), Ucat_x(i,j+1),Ucat_x(i,j),Ucat_x(i,j-1),0,1);
        Fpx2_j_m_1 =  Eval_Fp(j-1-4,N, Ucont_y(i,j-1),Ucat_x(i,j+1), Ucat_x(i,j),Ucat_x(i,j-1),Ucat_x(i,j-2),0,1);
        
        Conv_Flux_x = (Fpx1_i - Fpx1_i_m_1) / dx  +  (Fpx2_j - Fpx2_j_m_1) / dy;        
        
        vals(1) = vals(1) - Conv_Flux_x;
        
        %%-------------------- Point E ----------------------
        % Time derivative
        ii(2) = index_P_U;
        jj(2) = index_E_U;
        vals(2) = 0;
        
        % Viscous_Flux
        vals(2) = vals(2) + 1/(dx).^2;
        
        % Convective Flux
        Fpx1_i     =  Eval_Fe(i-4,M, Ucont_x(i,j),Ucat_x(i+2,j), Ucat_x(i+1,j),Ucat_x(i,j),Ucat_x(i-1,j),1,1);
        Fpx1_i_m_1 =  Eval_Fe(i-1-4,M, Ucont_x(i-1,j),Ucat_x(i+1,j), Ucat_x(i,j),Ucat_x(i-1,j),Ucat_x(i-2,j),1,1);
      
        Fpx2_j     =  Eval_Fe(j-4,N, Ucont_y(i,j),Ucat_x(i,j+2), Ucat_x(i,j+1),Ucat_x(i,j),Ucat_x(i,j-1),0,1);
        Fpx2_j_m_1 =  Eval_Fe(j-1-4,N, Ucont_y(i,j-1),Ucat_x(i,j+1), Ucat_x(i,j),Ucat_x(i,j-1),Ucat_x(i,j-2),0,1);
        
        Conv_Flux_x = (Fpx1_i - Fpx1_i_m_1) / dx  +  (Fpx2_j - Fpx2_j_m_1) / dy;        
        
        vals(2) = vals(2) - Conv_Flux_x;
        
        %%-------------------- Point W ----------------------
        % Time derivative
        ii(3) = index_P_U;
        jj(3) = index_W_U;
        vals(3) = 0;
        
        % Viscous_Flux
        vals(3) = vals(3) + 1/(dx).^2;
        
        % Convective Flux
        Fpx1_i     =  Eval_Fw(i-4,M, Ucont_x(i,j),Ucat_x(i+2,j), Ucat_x(i+1,j),Ucat_x(i,j),Ucat_x(i-1,j),1,1);
        Fpx1_i_m_1 =  Eval_Fw(i-1-4,M, Ucont_x(i-1,j),Ucat_x(i+1,j), Ucat_x(i,j),Ucat_x(i-1,j),Ucat_x(i-2,j),1,1);
      
        Fpx2_j     =  Eval_Fw(j-4,N, Ucont_y(i,j),Ucat_x(i,j+2), Ucat_x(i,j+1),Ucat_x(i,j),Ucat_x(i,j-1),0,1);
        Fpx2_j_m_1 =  Eval_Fw(j-1-4,N, Ucont_y(i,j-1),Ucat_x(i,j+1), Ucat_x(i,j),Ucat_x(i,j-1),Ucat_x(i,j-2),0,1);
        
        Conv_Flux_x = (Fpx1_i - Fpx1_i_m_1) / dx  +  (Fpx2_j - Fpx2_j_m_1) / dy;        
        
        vals(3) = vals(3) - Conv_Flux_x;
        
        %%-------------------- Point N ----------------------
        % Time derivative
        ii(4) = index_P_U;
        jj(4) = index_N_U;
        vals(4) = 0;
        
        % Viscous_Flux
        vals(4) = vals(4) + 1/(dy).^2;
        
        % Convective Flux
        
        Fpx2_j     =  Eval_Fe(j-4,N, Ucont_y(i,j),Ucat_x(i,j+2), Ucat_x(i,j+1),Ucat_x(i,j),Ucat_x(i,j-1),0,1);
        Fpx2_j_m_1 =  Eval_Fe(j-1-4,N, Ucont_y(i,j-1),Ucat_x(i,j+1), Ucat_x(i,j),Ucat_x(i,j-1),Ucat_x(i,j-2),0,1);
        
        Conv_Flux_x = (Fpx2_j - Fpx2_j_m_1) / dy;        
        
        vals(4) = vals(4) - Conv_Flux_x;
        %%-------------------- Point S ----------------------
        % Time derivative
        ii(5) = index_P_U;
        jj(5) = index_S_U;
        vals(5) = 0;
        
        % Viscous_Flux
        vals(5) = vals(5) + 1/(dy).^2;
        
        % Convective Flux
        
        Fpx2_j     =  Eval_Fw(j-4,N, Ucont_y(i,j),Ucat_x(i,j+2), Ucat_x(i,j+1),Ucat_x(i,j),Ucat_x(i,j-1),0,1);
        Fpx2_j_m_1 =  Eval_Fw(j-1-4,N, Ucont_y(i,j-1),Ucat_x(i,j+1), Ucat_x(i,j),Ucat_x(i,j-1),Ucat_x(i,j-2),0,1);
        
        Conv_Flux_x = (Fpx2_j - Fpx2_j_m_1) / dy;        
        
        vals(5) = vals(5) - Conv_Flux_x;

        %%-------------------- Point E_E ----------------------
        % Time derivative
        ii(6) = index_P_U;
        jj(6) = index_E_E_U;
        vals(6) = 0;
        
        % Convective Flux
        Fpx1_i     =  Eval_Fe_e(i-4,M, Ucont_x(i,j),Ucat_x(i+2,j), Ucat_x(i+1,j),Ucat_x(i,j),Ucat_x(i-1,j),0,1);
        Fpx1_i_m_1 =  Eval_Fe_e(i-1-4,M, Ucont_x(i-1,j),Ucat_x(i+1,j), Ucat_x(i,j),Ucat_x(i-1,j),Ucat_x(i-2,j),0,1);
      
        Conv_Flux_x = (Fpx1_i - Fpx1_i_m_1) / dx;        
        
        vals(6) = vals(6) - Conv_Flux_x;
        %%-------------------- Point W_W ----------------------
        % Time derivative
        ii(7) = index_P_U;
        jj(7) = index_W_W_U;
        vals(7) = 0;
        
        % Convective Flux
        Fpx1_i     =  Eval_Fw_w(i-4,M, Ucont_x(i,j),Ucat_x(i+2,j), Ucat_x(i+1,j),Ucat_x(i,j),Ucat_x(i-1,j),0,1);
        Fpx1_i_m_1 =  Eval_Fw_w(i-1-4,M, Ucont_x(i-1,j),Ucat_x(i+1,j), Ucat_x(i,j),Ucat_x(i-1,j),Ucat_x(i-2,j),0,1);
      
        Conv_Flux_x = (Fpx1_i - Fpx1_i_m_1) / dx;        
        
        vals(7) = vals(7) - Conv_Flux_x;
        %%-------------------- Point N_N ----------------------
        % Time derivative
        ii(8) = index_P_U;
        jj(8) = index_N_N_U;
        vals(8) = 0;
        
        % Convective Flux
        Fpx2_j     =  Eval_Fe_e(j-4,N, Ucont_y(i,j),Ucat_x(i,j+2), Ucat_x(i,j+1),Ucat_x(i,j),Ucat_x(i,j-1),0,1);
        Fpx2_j_m_1 =  Eval_Fe_e(j-1-4,N, Ucont_y(i,j-1),Ucat_x(i,j+1), Ucat_x(i,j),Ucat_x(i,j-1),Ucat_x(i,j-2),0,1);
        
        Conv_Flux_x = (Fpx2_j - Fpx2_j_m_1) / dy;       
        
        vals(8) = vals(8) - Conv_Flux_x;
        %%-------------------- Point S_S ----------------------
        % Time derivative
        ii(9) = index_P_U;
        jj(9) = index_S_S_U;
        vals(9) = 0;
        
        % Convective Flux
        Fpx2_j     =  Eval_Fw_w(j-4,N, Ucont_y(i,j),Ucat_x(i,j+2), Ucat_x(i,j+1),Ucat_x(i,j),Ucat_x(i,j-1),0,1);
        Fpx2_j_m_1 =  Eval_Fw_w(j-1-4,N, Ucont_y(i,j-1),Ucat_x(i,j+1), Ucat_x(i,j),Ucat_x(i,j-1),Ucat_x(i,j-2),0,1);
        
        Conv_Flux_x = (Fpx2_j - Fpx2_j_m_1) / dy;       
        
        vals(9) = vals(9) - Conv_Flux_x;
        
        %%-------------------- Point P_V ----------------------
        % Time derivative
        ii(10) = index_P_U;
        jj(10) = index_P_V;
        vals(10) = 0;
        
        % Convective Flux
        Fpx2_j     =  Eval_Fp(j-4,N, Ucont_y(i,j),Ucat_x(i,j+2), Ucat_x(i,j+1),Ucat_x(i,j),Ucat_x(i,j-1),1,0);
        Fpx2_j_m_1 =  Eval_Fp(j-1-4,N, Ucont_y(i,j-1),Ucat_x(i,j+1), Ucat_x(i,j),Ucat_x(i,j-1),Ucat_x(i,j-2),1,0);
        
        Conv_Flux_x = (Fpx2_j - Fpx2_j_m_1) / dy;       
        
        vals(10) = vals(10) - Conv_Flux_x;
       
        %%-------------------- Point N_V ----------------------
        % Time derivative
        ii(11) = index_P_U;
        jj(11) = index_N_V;
        vals(11) = 0;
        
        % Convective Flux
        Fpx2_j     =  Eval_Fp(j-4,N, Ucont_y(i,j),Ucat_x(i,j+2), Ucat_x(i,j+1),Ucat_x(i,j),Ucat_x(i,j-1),1,0);
        Fpx2_j_m_1 =  Eval_Fp(j-1-4,N, Ucont_y(i,j-1),Ucat_x(i,j+1), Ucat_x(i,j),Ucat_x(i,j-1),Ucat_x(i,j-2),1,0);
        
        Conv_Flux_x = (Fpx2_j - Fpx2_j_m_1) / dy;       
        
        vals(11) = vals(11) - Conv_Flux_x;

       %% End of first M*N components 
       %%-------------------------- Specify for M*N ->>> 2*M*N -----------------------------------
       
        %% ------------ Diagonal \partial F \partial u_i
        
        % Time derivative
        ii(12) = index_P_V;
        jj(12) = index_P_V;
        vals(12) = -1.5/dt;
        
        % Viscous_Flux
        vals(12) = vals(12) - 2/(dx).^2 - 2/(dy).^2;

      
        % \partial F partial u_p
        %Convective flux        
        Fpy1_i     =  Eval_Fp(i-4,M, Ucont_x(i,j),Ucat_y(i+2,j), Ucat_y(i+1,j),Ucat_y(i,j),Ucat_y(i-1,j),1,1);
        Fpy1_i_m_1 =  Eval_Fp(i-1-4,M, Ucont_x(i-1,j),Ucat_y(i+1,j), Ucat_y(i,j),Ucat_y(i-1,j),Ucat_y(i-2,j),1,1);
      
        Fpy2_j     =  Eval_Fp(j-4,N, Ucont_y(i,j),Ucat_y(i,j+2), Ucat_y(i,j+1),Ucat_y(i,j),Ucat_y(i,j-1),0,1);
        Fpy2_j_m_1 =  Eval_Fp(j-1-4,N, Ucont_y(i,j-1),Ucat_y(i,j+1), Ucat_y(i,j),Ucat_y(i,j-1),Ucat_y(i,j-2),0,1);
        
        Conv_Flux_y = (Fpy1_i - Fpy1_i_m_1) / dx  +  (Fpy2_j - Fpy2_j_m_1) / dy;        
        
        vals(12) = vals(12) - Conv_Flux_y;        
        
        
        %%-------------------- Point E_V ----------------------
        % Time derivative
        ii(13) = index_P_V;
        jj(13) = index_E_V;
        vals(13) = 0;
        
        % Viscous_Flux
        vals(13) = vals(13) + 1/(dx).^2;
        
        % Convective Flux
        Fpy1_i     =  Eval_Fe(i-4,M, Ucont_x(i,j),Ucat_y(i+2,j), Ucat_y(i+1,j),Ucat_y(i,j),Ucat_y(i-1,j),1,1);
        Fpy1_i_m_1 =  Eval_Fe(i-1-4,M, Ucont_x(i-1,j),Ucat_y(i+1,j), Ucat_y(i,j),Ucat_y(i-1,j),Ucat_y(i-2,j),1,1);
      
        Fpy2_j     =  Eval_Fe(j-4,N, Ucont_y(i,j),Ucat_y(i,j+2), Ucat_y(i,j+1),Ucat_y(i,j),Ucat_y(i,j-1),0,1);
        Fpy2_j_m_1 =  Eval_Fe(j-1-4,N, Ucont_y(i,j-1),Ucat_y(i,j+1), Ucat_y(i,j),Ucat_y(i,j-1),Ucat_y(i,j-2),0,1);
        
        Conv_Flux_y = (Fpy1_i - Fpy1_i_m_1) / dx  +  (Fpy2_j - Fpy2_j_m_1) / dy;        
        
        vals(13) = vals(13) - Conv_Flux_y;
        
        %%-------------------- Point W_V ----------------------
        % Time derivative
        ii(14) = index_P_V;
        jj(14) = index_W_V;
        vals(14) = 0;
        
        % Viscous_Flux
        vals(14) = vals(14) + 1/(dx).^2;
        
        % Convective Flux
        Fpy1_i     =  Eval_Fw(i-4,M, Ucont_x(i,j),Ucat_y(i+2,j), Ucat_y(i+1,j),Ucat_y(i,j),Ucat_y(i-1,j),1,1);
        Fpy1_i_m_1 =  Eval_Fw(i-1-4,M, Ucont_x(i-1,j),Ucat_y(i+1,j), Ucat_y(i,j),Ucat_y(i-1,j),Ucat_y(i-2,j),1,1);
      
        Fpy2_j     =  Eval_Fw(j-4,N, Ucont_y(i,j),Ucat_y(i,j+2), Ucat_y(i,j+1),Ucat_y(i,j),Ucat_y(i,j-1),0,1);
        Fpy2_j_m_1 =  Eval_Fw(j-1-4,N, Ucont_y(i,j-1),Ucat_y(i,j+1), Ucat_y(i,j),Ucat_y(i,j-1),Ucat_y(i,j-2),0,1);
        
        Conv_Flux_y = (Fpy1_i - Fpy1_i_m_1) / dx  +  (Fpy2_j - Fpy2_j_m_1) / dy;        
        
        vals(14) = vals(14) - Conv_Flux_y;
        
        %%-------------------- Point N_V ----------------------
        % Time derivative
        ii(15) = index_P_V;
        jj(15) = index_N_V;
        vals(15) = 0;
        
        % Viscous_Flux
        vals(15) = vals(11) + 1/(dy).^2;
        
        % Convective Flux
        
        Fpy2_j     =  Eval_Fe(j-4,N, Ucont_y(i,j),Ucat_y(i,j+2), Ucat_y(i,j+1),Ucat_y(i,j),Ucat_y(i,j-1),0,1);
        Fpy2_j_m_1 =  Eval_Fe(j-1-4,N, Ucont_y(i,j-1),Ucat_y(i,j+1), Ucat_y(i,j),Ucat_y(i,j-1),Ucat_y(i,j-2),0,1);
        
        Conv_Flux_y = (Fpy2_j - Fpy2_j_m_1) / dy;        
        
        vals(15) = vals(11) - Conv_Flux_y;
        %%-------------------- Point S_V ----------------------
        % Time derivative
        ii(16) = index_P_V;
        jj(16) = index_S_V;
        vals(16) = 0;
        
        % Viscous_Flux
        vals(16) = vals(16) + 1/(dy).^2;
        
        % Convective Flux
        
        Fpy2_j     =  Eval_Fw(j-4,N, Ucont_y(i,j),Ucat_y(i,j+2), Ucat_y(i,j+1),Ucat_y(i,j),Ucat_y(i,j-1),0,1);
        Fpy2_j_m_1 =  Eval_Fw(j-1-4,N, Ucont_y(i,j-1),Ucat_y(i,j+1), Ucat_y(i,j),Ucat_y(i,j-1),Ucat_y(i,j-2),0,1);
        
        Conv_Flux_y = (Fpy2_j - Fpy2_j_m_1) / dy;        
        
        vals(16) = vals(16) - Conv_Flux_y;

        %%-------------------- Point E_E_V ----------------------
        % Time derivative
        ii(17) = index_P_V;
        jj(17) = index_E_E_V;
        vals(17) = 0;
        
        % Convective Flux
        Fpy1_i     =  Eval_Fe_e(i-4,M, Ucont_x(i,j),Ucat_y(i+2,j), Ucat_y(i+1,j),Ucat_y(i,j),Ucat_y(i-1,j),0,1);
        Fpy1_i_m_1 =  Eval_Fe_e(i-1-4,M, Ucont_x(i-1,j),Ucat_y(i+1,j), Ucat_y(i,j),Ucat_y(i-1,j),Ucat_y(i-2,j),0,1);
      
        Conv_Flux_y = (Fpy1_i - Fpy1_i_m_1) / dx;        
        
        vals(17) = vals(17) - Conv_Flux_y;
        %%-------------------- Point W_W ----------------------
        % Time derivative
        ii(18) = index_P_V;
        jj(18) = index_W_W_V;
        vals(18) = 0;
        
        % Convective Flux
        Fpy1_i     =  Eval_Fw_w(i-4,M, Ucont_x(i,j),Ucat_y(i+2,j), Ucat_y(i+1,j),Ucat_y(i,j),Ucat_y(i-1,j),0,1);
        Fpy1_i_m_1 =  Eval_Fw_w(i-1-4,M, Ucont_x(i-1,j),Ucat_y(i+1,j), Ucat_y(i,j),Ucat_y(i-1,j),Ucat_y(i-2,j),0,1);
      
        Conv_Flux_y = (Fpy1_i - Fpy1_i_m_1) / dx;        
        
        vals(18) = vals(18) - Conv_Flux_y;
        %%-------------------- Point N_N ----------------------
        % Time derivative
        ii(19) = index_P_V;
        jj(19) = index_N_N_V;
        vals(19) = 0;
        
        % Convective Flux
        Fpy2_j     =  Eval_Fe_e(j-4,N, Ucont_y(i,j),Ucat_y(i,j+2), Ucat_y(i,j+1),Ucat_y(i,j),Ucat_y(i,j-1),0,1);
        Fpy2_j_m_1 =  Eval_Fe_e(j-1-4,N, Ucont_y(i,j-1),Ucat_y(i,j+1), Ucat_y(i,j),Ucat_y(i,j-1),Ucat_y(i,j-2),0,1);
        
        Conv_Flux_y = (Fpy2_j - Fpy2_j_m_1) / dy;       
        
        vals(19) = vals(19) - Conv_Flux_y;
        %%-------------------- Point S_S ----------------------
        % Time derivative
        ii(20) = index_P_V;
        jj(20) = index_S_S_V;
        vals(20) = 0;
        
        % Convective Flux
        Fpy2_j     =  Eval_Fw_w(j-4,N, Ucont_y(i,j),Ucat_y(i,j+2), Ucat_y(i,j+1),Ucat_y(i,j),Ucat_y(i,j-1),0,1);
        Fpy2_j_m_1 =  Eval_Fw_w(j-1-4,N, Ucont_y(i,j-1),Ucat_y(i,j+1), Ucat_y(i,j),Ucat_y(i,j-1),Ucat_y(i,j-2),0,1);
        
        Conv_Flux_y = (Fpy2_j - Fpy2_j_m_1) / dy;       
        
        vals(20) = vals(20) - Conv_Flux_y;
        
        %%-------------------- Point P_U ----------------------
        % Time derivative
        ii(21) = index_P_V;
        jj(21) = index_P_U;
        vals(21) = 0;
        
        % Convective Flux
        Fpy1_j     =  Eval_Fp(i-4,M, Ucont_x(i,j),Ucat_y(i,j+2), Ucat_y(i,j+1),Ucat_y(i,j),Ucat_y(i,j-1),1,0);
        Fpy1_j_m_1 =  Eval_Fp(i-1-4,M, Ucont_x(i-1,j),Ucat_y(i,j+1), Ucat_y(i,j),Ucat_y(i,j-1),Ucat_y(i,j-2),1,0);
        
        Conv_Flux_y = (Fpy1_j - Fpy1_j_m_1) / dx;       
        
        vals(21) = vals(21) - Conv_Flux_y;
       
        %%-------------------- Point N_U ----------------------
        % Time derivative
        ii(22) = index_P_V;
        jj(22) = index_N_U;
        vals(22) = 0;

        % Convective Flux
        Fpy1_j     =  Eval_Fp(i-4,M, Ucont_x(i,j),Ucat_y(i,j+2), Ucat_y(i,j+1),Ucat_y(i,j),Ucat_y(i,j-1),1,0);
        Fpy1_j_m_1 =  Eval_Fp(i-1-4,M, Ucont_x(i-1,j),Ucat_y(i,j+1), Ucat_y(i,j),Ucat_y(i,j-1),Ucat_y(i,j-2),1,0);
        
        Conv_Flux_y = (Fpy1_j - Fpy1_j_m_1) / dx;       
        
        vals(22) = vals(22) - Conv_Flux_y;
        
        
        % Assemble ii , jj , vals
        ct = 0;
        for k = 1:22
            if (ii(k) > 0 && jj(k) >0)
                ct = ct + 1;
             ii_0(counter+ct) = ii(k);
             jj_0(counter+ct) = jj(k);
             vals_0(counter+ct) = vals(k);
            end
        end
        counter = counter + ct;
         else
            index_P_U   = glidx(i0,j0,M,N)      ; 
            index_P_V   = glidx(i0,j0,M,N) +M*N ; 
            
            if index_P_U > 0
            counter = counter + 1;
            ii_0(counter) = index_P_U; 
            jj_0(counter) = index_P_U;
            vals_0(counter) = 1;
            end
            
            if index_P_V > 0
            counter = counter + 1;
            ii_0(counter) = index_P_V; 
            jj_0(counter) = index_P_V;
            vals_0(counter) = 1;
            end
        end % End of inside the domain       
       

        
    end % End of j
end % End of i
J = sparse(ii_0,jj_0,vals_0);