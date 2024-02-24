function Main()

clear all
MAXTIME =1000;
M = 20;
N = 20 ;
% Add two more ghost points
M2 = M+2;
N2 = N+2;

dt = 0.01;
Re = 1;
% Remember dU for RK
dU_x = 0;
dU_y = 0;
% Initialization process
[Ucont_x Ucont_y Ucat_x Ucat_y Ubcs_x Ubcs_y Pressure dx dy] = Init(M2,N2);%Init_Taylor_Green(M2,N2);

% It is the time integration scheme
for time_step = 1:MAXTIME       
       
    tic
    time_step
    t = 0;% time_step * dt;
    
    %[U_im_x U_im_y] =  Ex_Runge_Kutta(dU_x , dU_y, Ucont_x, Ucont_y, Ucat_x, Ucat_y, Ubcs_x, Ubcs_y, Pressure, Re, dx, dy, dt,t);            
     
     Vel_Star0 = [Vectorize(Ucont_x)  Vectorize(Ucont_y)];
    
     Vel_Star = fsolve(@myfun,Vel_Star0);
     
     U_im_x = Ucont_x;
     U_im_y = Ucont_y;
     
    [phi]  = Poisson_Solver(U_im_x, U_im_y, dx, dy, dt);
    
    [U_x_new U_y_new P_new] = Update_Solution(U_im_x, U_im_y, Pressure, phi,dx,dy,dt);
    
    %Update dU
    dU_x = U_x_new - Ucont_x;
    dU_y = U_y_new - Ucont_y;
    
    % Go next time step
    [Ucat_x Ucat_y Ucont_x Ucont_y] = FormBCS(U_x_new, U_y_new, Ubcs_x, Ubcs_y, dx, dy,Re,t);
    
    Pressure = P_new;          
        
    % Assess the divergence of flow field
    Div = Divergence(U_x_new, U_y_new,dx ,dy);
    MaxDiv = norm(Div,inf)
   
    toc             
end

    % Display flow field
    figure(1);        
    quiver(Ucat_x,Ucat_y,1);
    
    % Error 
    [Exact_x Exact_y] = Taylor_Green(M2,N2,t);    
    
    for i = 1 :M2
        for j = 1:N2
            
            index = glidx(i,j,M2,N2);
            
            if (i>=2 && i<=M-1 && j >=2 && j<=N-1)
            e(index) = Ucat_x(i,j) - Exact_x(i,j);
            else
                e(index) = 0;
            end
        end
    end
       
    norm(e,2) 
    
    figure(2);
    plot(e);
   
  