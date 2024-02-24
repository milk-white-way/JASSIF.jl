function [U_x_new U_y_new P_new] = Update_Solution(Ucont_x, Ucont_y, Pressure, phi,dx,dy,dt)

% Return to the old things
M = length(Ucont_x(:,1));
N = length(Ucont_x(1,:));

% Return the phi at integer node
for im = 1:M*N
    [i j] = lidx(im,M,N);
    a_Phi(i,j) = phi(im);
end

% Calculate pressure gradient at half node
P_G_x(1:M,1:N) = 0;

for i = 2:M-2
    for j = 2:N-1
   
        
       P_G_x(i,j) = ( a_Phi(i+1,j) - a_Phi(i,j)) / dx;       
       
    end
end

P_G_y(1:M,1:N) = 0;
for i = 2:M-1
    for j = 2:N-2
        
       P_G_y(i,j) =  (a_Phi(i,j+1) - a_Phi(i,j)) / dy;
       
    end
end


% Update contravariant x & y
U_x_new = Ucont_x;
U_y_new = Ucont_y;

for i = 1:M-1
    for j = 1:N
        
       U_x_new(i,j) = Ucont_x(i,j) - P_G_x(i,j)*dt /1.5;
       
    end
end

for i = 1:M
    for j = 1:N-1
        
    U_y_new(i,j) = Ucont_y(i,j) -    P_G_y(i,j)*dt /1.5;
    end
end

% Update pressure
P_new = Pressure + a_Phi;
