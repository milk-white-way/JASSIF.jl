function [RHS_x RHS_y] = RHS_Taylor_Green(M,N, time)

dx = pi / M;
dy = pi / N;


for i =1:M
 for j = 1:N     
     
        x = (i-2+0.5) * dx;
        y = (j-2) * dy;                
        
        RHS_x(i,j) =  -2*(-cos(x) * sin(y) * exp(-2 *time));     
        
        if (i==M ||  i==1 ||j==1 ||j==N ||i==M-1)
            RHS_x(i,j) = 0;
        end
 end
end


    for i = 1:M    
        for j=1:N        
        
        x = (i-2) * dx;
        y = (j-2+0.5) * dy;        
        
        RHS_y(i,j) = -2*(sin(x) * cos(y) * exp(-2*time));
        
        if (j==N ||j==1||i==1||i==M ||j==N-1)
            RHS_y(i,j) = 0;
        end
        
        end
    end