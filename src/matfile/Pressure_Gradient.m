function [P_Gradient_x P_Gradient_y] = Pressure_Gradient(Pressure,dx,dy)

M = length(Pressure(:,1));
N = length(Pressure(1,:));

for i = 1:M
    for j = 1:N
        
        % ------ x direction ------------------------------
        if (i >=2 && i<=M-1)
            P_Gradient_x(i,j) = ( Pressure(i+1,j) - Pressure(i-1,j)) / (2*dx);
        else
            if (i==2)
                P_Gradient_x(i,j) = (Pressure(i+1,j) - Pressure(i,j)) / (dx);
            else
                if (i== M-1)
                    P_Gradient_x(i,j) = (Pressure(i,j) - Pressure(i-1,j)) / (dx);
                end
            end            
        end
        
        % ---------- y direction --------------------------
        if (j >=2 && j<=N-1)
            P_Gradient_y(i,j) = ( Pressure(i,j+1) - Pressure(i,j-1)) / (2*dy);
        else
            if (j==2)
                P_Gradient_y(i,j) = (Pressure(i,j+1) - Pressure(i,j)) / (dy);
            else
                if (j== N-1)
                    P_Gradient_y(i,j) = (Pressure(i,j) - Pressure(i,j-1)) / (dy);
                end
            end            
        end
        
        if (i==1 || i == M || j==1 || j==N) 
            P_Gradient_x(i,j) = 0;          
            P_Gradient_y(i,j) = 0;
        end    
        
    end
end
