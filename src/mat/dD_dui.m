function [coef]  = dD_dui(i,j,M,N,Re,dx,dy,direction);
%coef(1) = i-2;
%coef(2) = i-1;
%coef(3) = i+1;
%coef(4) = i+2;
% --  coef(5) = i,j
%coef(6) = j-2;
%coef(7) = j-1;
%coef(8) = j+1;
%coef(9) = j+2;
%-----cross term
%coef(10) = i-2 j =
%coef(11) = i-1, j -1
%coef(12) = i-2,j i.e - 2*N
% direction = 0 for x
% direction = 1 for y

for k = 1:15
    coef(k) = 0;
end

% i,j --- depending on the average of the direction
if direction == 0
  if (j ~= 1 && j~=N)
            coef(5) = -(1/Re) * ( 0.5/(dx.^2) + 1/(dy.^2)); % approximated for u(i,j)
   end
        
end

if direction == 1    
        if (i ~= 1 && i~=M)
            coef(5) =  - (1/Re) * ( 1/(dx.^2) + 0.5/(dy.^2));            
        end
        
end

%------------- i direction  ----------
% i - 1 ,j 
if direction == 0
  if (j ~= 1 && j~=N)
            coef(2) = (1/Re) * ( -0.5/(dy.^2));
   end
        
end

% i + 1,j
if direction == 0
  if (j ~= 1 && j~=N)
            coef(3) = (1/Re) * ( -0.5/(dy.^2));
   end
        
end
% i  ,j -1
if direction == 0
  if (i ~= 1 && i~=M)
            coef(7) = (1/Re) * ( 0.5/(dy.^2)); 
   end
        
end

% i  ,j +1
if direction == 0
  if (i ~= 1 && i~=M)
            coef(8) = (1/Re) * ( 0.5/(dy.^2)); 
   end
        
end

% i+1  ,j - 1
if direction == 0
  if (i ~= 1 && i~=M)
            coef(10) = (1/Re) * ( 0.25/(dy.^2)); 
   end
        
end

% i+1  ,j +1
if direction == 0
  if (i ~= 1 && i~=M)
            coef(11) = (1/Re) * ( 0.25/(dy.^2)); 
   end
        
end

% i-2  ,j
if direction == 0
  if (i ~= 1 && i~=M)
            coef(12) = (1/Re) * ( 0.25/(dx.^2)); 
   end
        
end

%---------- for j direction - ------------------
% i-2  ,j
if direction == 1
  if (i ~= 1 && i~=M)
            coef(13) = (1/Re) * ( 0.5/(dx.^2)); 
   end
        
end

% i  ,j-1
if direction == 1
  if (i ~= 1 && i~=M)
            coef(14) = -(1/Re) * ( 0.5/(dx.^2)); 
   end
        
end
% i-1  ,j-1
if direction == 1
  if (i ~= 1 && i~=M)
            coef(15) = (1/Re) * ( 0.25/(dx.^2)); 
   end
        
end




