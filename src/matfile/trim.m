function [E_x E_y] = trim(Ucat_x, Ucat_y)

M2 = length(Ucat_x(:,1));
N2 = length(Ucat_x(1,:));

 for i= 1:M2
        for j =1:N2   
            
            if (i>1 && i<M2 && j>1 && j<N2)
               
                E_x(i-1,j-1) = Ucat_x(i,j);
                E_y(i-1,j-1) = Ucat_y(i,j);
            end
            
            
        end
    end