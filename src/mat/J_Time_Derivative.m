function J = Time_Derivative(x)
M = length(x.Ucont_x(:,1));
N = length(x.Ucont_x(1,:));

J = spalloc(2*M*N,2*M*N, 22 * M*N);

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
        
       

         for k=1:22
            ii(k) = 0;
            jj(k) = 0;
            vals(k)  = 0;
         end
      
        %% ------------ Diagonal \partial F \partial u_i
        % Time derivative
        ii(3) = index_P_U;
        jj(3) = index_P_U;
                       
        ii(12) = index_P_V;
        jj(12) = index_P_V;

        
        % The the boundary parts
        if (i ~= 1 && i ~= M-1 && i~=M && j ~=1 && j~=N)
            vals(3) = -1.5/dt;
        else
            vals(3) = 1;
        end
              
        if (j ~= 1 && j ~= N-1 && j~=N && i ~=1 && i~=M)
            vals(12) = -1.5/dt;
        else
            vals(12) = 1;        
        end

        %%----------------------------------------------------------------
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
    end % End of j
end % End of i
J = sparse(ii_0,jj_0,vals_0);