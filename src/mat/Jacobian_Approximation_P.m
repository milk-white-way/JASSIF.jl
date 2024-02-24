function J = Jacobian_Approximation_P(RHS, x)
M = length(x.Ucont_x(:,1));
N = length(x.Ucont_x(1,:));


n= 100;

V = speye(2*M*N);


J = [];

% Current evalutation of function
counter = 0;
P   = feval(RHS,x);

for i = 1:n

x_not = x;
h = V(:,i);
epsilon = 1e-14;

% auxillary vector - to evaluate the contravariant vector
if i <= M*N
    cter = i;
    [id_x id_y] = lidx(cter,M,N);
    x_not.U_im_x(id_x,id_y) =  x.U_im_x(id_x,id_y) + epsilon;     
else
    cter = i - M*N;
    [id_x id_y] = lidx(cter,M,N);
    x_not.U_im_y(id_x,id_y) =  x.U_im_y(id_x,id_y) + epsilon;
end


P_1 = feval(RHS,x_not);


%% -------------  Exclude the boundary terms
if i <=M*N % Ucont_x
        if (id_x > 1 && id_x < M-1)
            flag =0;
        else
            flag = 1;
        end
        
        if (id_y == 1 || id_y == N)
            flag = 1;
        end
else % Ucont_y
        if (id_y > 1 && id_y < N-1)
            flag = 0;
        else
            flag = 1;
        end
        
        if (id_x == 1 || id_x == M)
            flag = 1;
        end
end

%% -------------  
if (flag == 0)
    j_s = (P_1 - P) / epsilon;
else
    j_s = h;
end
% Pick only three terms

if (flag == 0)
    
    if (i >1)
        counter= counter+1;
        ii(counter)  = i-1;
        jj(counter)  = i;
        vals(counter) = j_s(i-1);
    end

    counter= counter+1;
    ii(counter)  = i;
    jj(counter)  = i;
    vals(counter) = j_s(i);

    if (i<n)
        counter= counter+1;
        ii(counter)  = i+1;
     jj(counter)  = i;
     vals(counter) = j_s(i+1);
    end
else
    

    counter= counter+1;
    ii(counter)  = i;
    jj(counter)  = i;
    vals(counter) = 1;
    
end

J = [J,j_s];
end


for i = n+1:2*M*N

    h = zeros(2*M*N,1);
    J = [J,h];
end


% If use only 3 diagonal terms
%J = sparse(ii,jj,vals);
