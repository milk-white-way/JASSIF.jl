function Solution_Comparison();

[Ucat_x Ucat_y Pressure dx dy M N t] = Main();

%
M2 = M+ 2;
N2 = N+ 2;
% Find the Solution
    Solution_x(1:M,1:N) = 0;
    Solution_y(1:M,1:N) = 0;
    
    
    for i= 1:M2
        for j =1:N2   
            
                x(i,j) = (i-2)*dx;
                y(i,j) = (j-2)*dy;  
            if (i>1 && i<M2 && j>1 && j<N2)
               
                x1(i-1,j-1) = (i-2)*dx;
                y1(i-1,j-1) = (j-2)*dy;   
                Solution_x(i-1,j-1) = Ucat_x(i,j);
                Solution_y(i-1,j-1) = Ucat_y(i,j);
            end
            
            
        end
    end
    
    grid on
    figure(1)
    quiver(x1,y1,Solution_x,Solution_y,2);
    title(' Velocity field');
    
    [omega_z, cav] = curl(Solution_x, Solution_y);
    
    % Plot vorticity 
    figure(2)
    level = [-4 -3 -2 -1  0 1 2 3 4 5 6];
    [c h] = contour(x1,y1,omega_z,10);
    clabel(c,h);
    title('Vorticity field');
    
    % --------------- Ghia's data --------------
    y_Ghia         = [0  0.0547   0.0625   0.0703   0.1016   0.1719   0.2813   0.4531   0.5      0.6172   0.7344   0.8516   0.9531  0.9609  0.9688  0.9766  1];
    Ghia_U_100     = [0  -0.03717 -0.04192 -0.04775 -0.06434 -0.10150 -0.15662 -0.21090 -0.20581 -0.13641 0.00332  0.23151  0.68717 0.73722 0.78871 0.84123 1];
    
    x_Ghia         = [0 0.0625  0.0703  0.0781  0.0938  0.1563  0.2266  0.2344  0.5      0.8047   0.8594   0.9063   0.9453   0.9531   0.9609   0.9688   1];
    Ghia_V_100     = [0 0.09233 0.10091 0.10890 0.12317 0.16077 0.17507 0.17527 0.05454  -0.24533 -0.22445 -0.16914 -0.10313 -0.08864 -0.07391 -0.05906 0];
   
    %Plot u profile
    figure(3)
    center = floor(M/2);
    u_center = Ucat_x(center,:);
    y_center = y(center,:);
    plot(y_center,u_center);
    hold on
    plot(y_Ghia, Ghia_U_100,'r .');
    
    %Plot u profile
    figure(4)
    center = floor(N/2);
    v_center = Ucat_y(:,center);
    x_center = x(:,center);
    plot(x_center,v_center);
    hold on
    plot(x_Ghia, Ghia_V_100,'r .');
    
    % --------------  Error checking ------------------- 
    %[Exact_x Exact_y] = Taylor_Green(M2,N2,t);  
    
    %Error = Ucat_x - Exact_x;
    %for i=1:M2
    %    for j= 1:N2
    
    %       if (i==1 ||j==1 ||i==M2 ||j==N2)
    %           Error(i,j) = 0;
    %        end
    %    end
    %end

    %figure(3)
    %mesh(Error)
    %plot(Vectorize(Exact_x));
    %hold on
    %plot(Vectorize(Ucat_x),'r o');
    
    %figure(4)
    %plot(Vectorize(Error)); 
   
   