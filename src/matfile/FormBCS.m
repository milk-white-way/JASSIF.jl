function [Ucat_x Ucat_y Ucont_new_x Ucont_new_y] = FormBCS(Ucont_x, Ucont_y, Ubcs_x, Ubcs_y, dx, dy,Re,time)

M = length(Ucont_x(:,1));
N = length(Ucont_x(1,:));

% Apply boundary condition
Ucont_new_x = Ucont_x;
Ucont_new_y = Ucont_y;

% Setup contravariant component ---------------- U cont
    for j = 1:N
       
        %i=1; 
        %x = (i-2 +0.5) * dx;
        %y = (j-2) * dy;                
        %Ucont_new_x(i,j) =  -cos(x) * sin(y) * exp(-2 *time);
        
        
        %i=M-1; 
        %x = (i-2+0.5) * dx;
        %y = (j-2) * dy;                
        %Ucont_new_x(i,j) =  -cos(x) * sin(y) * exp(-2*time);
                
        Ucont_new_x(1,j) = 0;
        Ucont_new_x(M-1,j) = 0;        
    end



    for i = 1:M    
        
        %j = 1;
        %x = (i-2) * dx;
        %y = (j-2+0.5) * dy;        
        %Ucont_new_y(i,j) = sin(x) * cos(y) * exp(-2*time);
        
        %j = N-1;
        %x = (i-2) * dx;
        %y = (j-2+0.5) * dy;        
        %Ucont_new_y(i,j) = sin(x) * cos(y) * exp(-2*time);
        
        Ucont_new_y(i,1) = 0;
        Ucont_new_y(i,N-1) = 0;
    end


% Set boundary condition - Ubcs ----------------------
    for j=1:N      
        
        i=1; 
        x = (i-2 +0.5) * dx;
        y = (j-2) * dy;                
        Ubcs_x_3(j) =  0;%-cos(x) * sin(y) * exp(-2 *time);
        Ubcs_y_3(j) =  0;% sin(x) * cos(y) * exp(-2*time);    
        
        i=M-1; 
        x = (i-2 +0.5) * dx;
        y = (j-2) * dy;                
        Ubcs_x_4(j) = 0;%-cos(x) * sin(y) * exp(-2 *time);
        Ubcs_y_4(j) = 0;% sin(x) * cos(y) * exp(-2*time);    
        
                
    end


    for i=1:M                
     
        
        j=1; 
        x = (i-2) * dx;
        y = (j-2 +0.5) * dy;                
        Ubcs_x_1(i) =  0;%-cos(x) * sin(y) * exp(-2 *time);
        Ubcs_y_1(i) =  0;% sin(x) * cos(y) * exp(-2*time);    
        
        j=N-1; 
        x = (i-2) * dx;
        y = (j-2 +0.5) * dy;                
        Ubcs_x_2(i) =  1;%-cos(x) * sin(y) * exp(-2 *time);
        Ubcs_y_2(i) =  0;% sin(x) * cos(y) * exp(-2*time);    
               
    end



% Convert to Cartesian components
[Ucat_x Ucat_y] = Contra_To_Cart(Ucont_new_x, Ucont_new_y);

% Mirror to the ghost nodes
for i = 1:M
    Ucat_x(i,1) = 2 * Ubcs_x_1(i) - Ucat_x(i,2);
    Ucat_y(i,1) = 2 * Ubcs_y_1(i) - Ucat_y(i,2);
    
    Ucat_x(i,N) = 2 * Ubcs_x_2(i) - Ucat_x(i,N-1);
    Ucat_y(i,N) = 2 * Ubcs_y_2(i) - Ucat_y(i,N-1);
    
end


for j = 1:N
    Ucat_x(1,j) = 2 * Ubcs_x_3(j) - Ucat_x(2,j);
    Ucat_y(1,j) = 2 * Ubcs_y_3(j) - Ucat_y(2,j);
    
    Ucat_x(M,j) = 2 * Ubcs_x_4(j) - Ucat_x(M-1,j);
    Ucat_y(M,j) = 2 * Ubcs_y_4(j) - Ucat_y(M-1,j);
    
end
% Apply directly the boundaries
%for i=1:M
%    for j = 1:N
        
%        x = (i-2) * dx;
%        y = (j-2) * dy;
        
%        if (i==1 || i==M ||j==1 ||j==N)
%            Ucat_x(i,j) =  -cos(x) * sin(y)* exp(-2 *time);
%            Ucat_y(i,j) =   sin(x) * cos(y) * exp(-2 *time);
%        end
%    end
%end


% Zero corner points
Ucat_x(1,1) = 0;
Ucat_y(1,1) = 0;

Ucat_x(1,N) = 0;
Ucat_y(1,N) = 0;

Ucat_x(M,1) = 0;
Ucat_y(M,1) = 0;

Ucat_x(M,N) = 0;
Ucat_y(M,N) = 0;