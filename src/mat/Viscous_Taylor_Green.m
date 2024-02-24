function [Viscous_x Viscous_y] = Viscous_Taylor_Green(M,N, time)

dx = pi / M;
dy = pi / N;


for i =1:M
 for j = 1:N     
     
        x = (i-2) * dx;
        y = (j-2) * dy;                
        
        Viscous_x(i,j) =  -2*(-cos(x) * sin(y) * exp(-2 *time));     
        
        if (i==M ||  i==1 ||j==1 ||j==N)
            Viscous_x(i,j) = 0;
        end
 end
end


    for i = 1:M    
        for j=1:N        
        
        x = (i-2) * dx;
        y = (j-2) * dy;        
        
        Viscous_y(i,j) = -2*(sin(x) * cos(y) * exp(-2*time));
        
        if (j==N ||j==1||i==1||i==M)
            Viscous_y(i,j) = 0;
        end
        
        end
    end