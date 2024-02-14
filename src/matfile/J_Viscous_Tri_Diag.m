function J = J_Viscous_Contribution(x)
M = length(x.Ucont_x(:,1));
N = length(x.Ucont_x(1,:));

J = spalloc(2*M*N,2*M*N, 3 * M*N);

dx      = x.dx;
dy      = x.dy;
dt      = x.dt;
Re      = x.Re;
Ucat_x  = x.Ucat_x;
Ucat_y  = x.Ucat_y;

Ucont_x  = x.Ucont_x;
Ucont_y  = x.Ucont_y;

counter = 0;

%---------------  First M*N values of F ----------------------------------
for i = 1:M
    for j =1:N             
        
         

             
        % Set up all the indexes
        index_N_N_U = glidx(i,j-2,M,N); %1
        index_N_U   = glidx(i,j-1,M,N); %2
        index_P_U   = glidx(i,j,M,N)  ; %3
        index_S_U   = glidx(i,j+1,M,N); %4
        index_S_S_U = glidx(i,j+2,M,N); %5
        
        index_W_U   = glidx(i+1,j,M,N); %6
        index_W_W_U = glidx(i+2,j,M,N); %7       
        index_E_U   = glidx(i-1,j,M,N); %8
        index_E_E_U = glidx(i-2,j,M,N); %9 
        
     
        index_N_N_V = glidx(i,j-2,M,N) + M*N; %10
        index_N_V   = glidx(i,j-1,M,N) + M*N; %11
        index_P_V   = glidx(i,j,M,N)   + M*N; %12
        index_S_V   = glidx(i,j+1,M,N) + M*N; %13
        index_S_S_V = glidx(i,j+2,M,N) + M*N; %14
        
        index_W_V   = glidx(i+1,j,M,N) + M*N; %15
        index_W_W_V = glidx(i+2,j,M,N) + M*N; %16       
        index_E_V   = glidx(i-1,j,M,N) + M*N; %17
        index_E_E_V = glidx(i-2,j,M,N) + M*N; %18   
        
       
         
        Viscous_coeff_x = dD_dui(i,j,M,N,Re,dx,dy,0);
        Viscous_coeff_y = dD_dui(i,j,M,N,Re,dx,dy,1);
        
        %% ------------ Diagonal \partial F \partial u_i
        % Time derivative
        counter= counter + 1;
        ii(counter) = index_P_U;
        jj(counter) = index_P_U;
        
        % Viscous contribution
        if (j ~= 1 && j~=N && i ~=1 && i ~=M && i~=M-1)
            vals(counter) = Viscous_coeff_x(5);
        else
            vals(counter) = 0;
        end
        
        %P
        counter = counter + 1;
        ii(counter) = index_P_V;
        jj(counter) = index_P_V;
        
        % Viscous contribution
        if (i ~= 1 && i~=M && j ~=1 && j~=N && j~=N-1)
            vals(counter) = Viscous_coeff_y(5);
        else
            vals(counter) = 0;
        end
        %% -------------- x direction components -------------------------
        % Off diagonal component % i - 1-----------------------------------
        % Viscous contribution
        
        %if (i >2 && i<M-1 && j ~= 1 && j~=N)
           
        %    counter = counter + 1;        
            
        %    ii(counter) = index_P_U;
        %    jj(counter) = index_E_U;
            
        %    vals(counter) =  Viscous_coeff_x(2);
            
        %end      
        % Off diagonal component % i + 1-----------------------------------
        % Viscous contribution
        
        %if (i >1 && i<M-2 && j ~= 1 && j~=N)
        
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_U;
        %    jj(counter) = index_W_U;
            
        %    vals(counter) =  Viscous_coeff_x(3);
            
        %end      
         
        % Off diagonal component % j-1-------------------------------------
        % Viscous contribution
        %if (i >1 && i<M-1 && j >2 && j~=N)
            
        %    counter = counter + 1;
            
        %    ii(counter) = index_P_U;
        %    jj(counter) = index_N_U;
            
        %    vals(counter) =  Viscous_coeff_x(7);
            
        %end      
        
        % Off diagonal component % j+1-------------------------------------
        % Viscous contribution
        %if (i >1 && i<M-1 && j >1 && j<N-1)
            
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_U;
        %    jj(counter) = index_S_U;
            
        %    vals(counter) =  Viscous_coeff_x(8);
            
        %end      
        
        % Off diagonal component % i + 1,j+1-----------------------------------
        % Viscous contribution
        
        %if (i >1 && i<M-2 && j >1 && j<N-1)
        
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_U;
        %    jj(counter) = index_W_U+1;
            
        %    vals(counter) = Viscous_coeff_x(11);
            
        %end      
        % Off diagonal component % i + 1,j-1-----------------------------------
        % Viscous contribution
        
        %if (i >2 && i<M-1 && j >2 && j<N)
            
        %   counter = counter + 1;
        
        %    ii(counter) = index_P_U;
        %    jj(counter) = index_E_U-1;
            
        %    vals(counter) = Viscous_coeff_x(10);
            
        %end      
        % Off diagonal component % i - 1,j+1-----------------------------------
        % Viscous contribution
        
        %if (i >1 && i<M-2 && j >2 && j<N)
        
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_U;
        %    jj(counter) = index_W_U-1;
            
        %    vals(counter) = Viscous_coeff_x(10);
            
        %end      
        
        % Off diagonal component % i + 1,j+1-------------------------------
        %if (i >2 && i<M-1 && j >1 && j<N-1)
        
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_U;
        %    jj(counter) = index_E_U+1;
            
        %    vals(counter) = Viscous_coeff_x(10);
            
        %end      
        % Off diagonal component % i - 2,j-----------------------------------
        % Viscous contribution
        
        % if (i >3 && i<M-1 && j >1 && j<N)
        % 
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_U;
        %    jj(counter) = index_P_U-2*N;
            
        %    vals(counter) = Viscous_coeff_x(12);
            
        %end      
        % Off diagonal component % i + 2,j-----------------------------------
        % Viscous contribution
        
        %if (i >1 && i<M-3 && j >1 && j<N)
        % 
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_U;
        %   jj(counter) = index_P_U+2*N;
            
        %    vals(counter) = Viscous_coeff_x(12);
            
        %end      
        
        
        %% -------------- y direction components-------------------------
        %% --------------------------------------------------------------
        % Off diagonal component % i - 1-----------------------------------
        % Viscous contribution
        
        %if (i >2 && i<M && j >1 && j<N-1)
        
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_V;
        %    jj(counter) = index_E_V;
        %    
        %    vals(counter) = Viscous_coeff_y(13);
            
        %end      
        
        % Off diagonal component % i + 1-----------------------------------
        % Viscous contribution
        
        %if (i >1 && i<M-1 && j >1 && j<N-1)
            
        %    counter = counter + 1;
         
        %    ii(counter) = index_P_V;
        %    jj(counter) = index_W_V;
            
        %    vals(counter) = Viscous_coeff_y(13);
            
        %end      
        % Off diagonal component % j - 1-----------------------------------
        % Viscous contribution
        
        %if (i >1 && i<M && j >2 && j<N-1)
            
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_V;
        %    jj(counter) = index_N_V;
            
        %    vals(counter) = Viscous_coeff_y(14);
            
        %end      
        % Off diagonal component % j + 1-----------------------------------
        % Viscous contribution
        
        %if (i >1 && i<M && j >1 && j<N-2)
        
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_V;
        %    jj(counter) = index_S_V;
           
        %    vals(counter) = Viscous_coeff_y(14);
            
        %end
        %----------------- i - 1, j -1
        % Viscous contribution        
        %if (i >2 && i<M && j >2 && j<N-1)
        
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_V;
        %    jj(counter) = index_E_V -1;
            
        %    vals(counter) = Viscous_coeff_y(14);
            
        %end      
        %----------------- i + 1, j -1
        % Viscous contribution        
        %if (i >1 && i<M-1 && j >1 && j<N-2)
        
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_V;
        %    jj(counter) = index_W_V +1;
            
        %    vals(counter) = Viscous_coeff_y(14);
            
        %end      
        %----------------- i - 1, j+1
        % Viscous contribution        
        %if (i >2 && i<M && j >1 && j<N-2)
            
        %    counter = counter + 1;
        
        %   ii(counter) = index_P_V;
        %    jj(counter) = index_E_V +1;
            
        %    vals(counter) = Viscous_coeff_y(14);
            
        %end      
        
        %----------------- i + 1, j-1
        % Viscous contribution        
        %if (i >1 && i<M-1 && j >2 && j<N-1)
        
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_V;
        %    jj(counter) = index_W_V -1;
            
        %    vals(counter) = Viscous_coeff_y(15);
            
        %end      
        %----------------- i, j-2
        % Viscous contribution        
        %if (i >1 && i<M && j >3 && j<N-1)
        
        %    counter = counter + 1;
        
        %    ii(counter) = index_P_V;
        %    jj(counter) = index_N_V - 1;
            
        %    vals(counter) =  Viscous_coeff_y(15);
            
        %end      
        
        %----------------- i, j+2
        % Viscous contribution        
        %if (i >1 && i<M && j >1 && j<N-3)
        
        %   counter = counter + 1;
        
        %    ii(counter) = index_P_V;
        %    jj(counter) = index_S_V + 1;
            
        %    vals(counter) =  Viscous_coeff_y(15);
            
        %end      
        
        
    end % End of j
end % End of i

J = sparse(ii,jj,vals);