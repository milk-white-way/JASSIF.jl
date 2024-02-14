function [E_x E_y] = Taylor_Green(M,N,time)

dx = pi / (M-3);
dy = pi / (N-3);


for i=1:M
    for j = 1:N
        
        x = (i-2) * dx;
        y = (j-2) * dy;
        
    E_x(i,j) = -cos(x) * sin(y) *exp(-2*time);
    E_y(i,j) = sin(x) * cos(y)  * exp(-2*time);
    end
end
 