function [Output_x Output_y] = Un_Vectorize(vector, M,N)


for counter = 1:length(vector) 
    
  if counter <= M*N
      index = counter;
      [i j] = lidx(index,M,N);       
      Output_x(i,j) = vector(counter);    
  else
      index = counter - M*N;
      [i j] = lidx(index,M,N);       
      Output_y(i,j) = vector(counter);    
  end
      
end
