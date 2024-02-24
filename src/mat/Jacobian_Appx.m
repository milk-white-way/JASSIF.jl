function J = Jacobian_Appx(RHS, x)
n = length(x);
V = speye(n);

J = [];

% Current evalutation of function
counter = 0;
P   = feval(RHS,x);

for i = 1:n
h = V(:,i);
epsilon = 1e-9;

P_1 = feval(RHS,x+epsilon*h);
j_s = (P_1 - P) / epsilon;

% Pick only three terms

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

end

% If use only 3 diagonal terms
J = sparse(ii,jj,vals);
