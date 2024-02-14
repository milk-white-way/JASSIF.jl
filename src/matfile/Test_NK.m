function Test_NK


 
tic
[x1 iter] = Newton_Krylov('myfunc',2,rand(100,1),30,1e-6)
%[x1 iter] = Newton_Krylov_Frechet('myfunc',2,0.5*ones(100,1),30,1e-6)
fx = norm(feval('myfunc',x1),inf)
toc