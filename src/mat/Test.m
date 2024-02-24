function Main()

MAXTIME = 10;
M = 20;
N = 20;
M2 = M + 2;
N2 = N + 2;
dt = 0.01;
Re = 1;
% Remember dU for RK
dU_x = 0;
dU_y = 0;
% Initialization process
[Ucont_x Ucont_y Ucat_x Ucat_y Ubcs_x Ubcs_y Pressure dx dy] = Init(M,N);

quiver(Ucat_x, Ucat_y);

for i = 1:M2-1
    for j = 1:N2
        
        x = (i-2+0.5) * dx;
        y = (j-2) * dy;                
        
        Exact_x(i,j) = - 2 *  sin(x) * cos(y);
        
    end
end


for i = 1:M2
    for j = 1:N2-1
        
        x = (i-2) * dx;
        y = (j-2+0.5) * dy;        
        
        Exact_y(i,j) = -2* (- cos(x)*sin(y));
    end
end

% Zero out the boundaries
  for i = 1 :M2
      Exact_x(i,1) = 0;
      Exact_x(i,N2) = 0;
      Exact_x(i,N2-1) = 0;
  end
  for j=1:N2
      Exact_x(M2,j) = 0;
      Exact_x(M2-1,j) = 0;
      Exact_x(1,j) = 0;
  end
    
  for i=1:M2
      Exact_y(i,N2-1) = 0;
      Exact_y(i,1) = 0;
      Exact_y(i,N2) = 0;
  end
  
  for j=1:N2
      Exact_y(1,j) =0;
      Exact_y(M2,j) = 0;
      Exact_y(M2-1,j) = 0;
  end
%% -- Test the convergence of RHS

[RHS_x RHS_y] = RHS_Calculation(Ucont_x, Ucont_y, Ucat_x, Ucat_y, Pressure,Re,dx,dy);

plot(Exact_y','b x');
%mesh(Exact_x - Convective_Flux_x);

hold on
plot(RHS_y','r o');
%mesh(F_x);


error =0;
for i=1:M2
    for j = 1:N2
       error = error + ( Exact_y(i,j) - RHS_y(i,j)).^2;
    end
 end

norm(sqrt(error),2)