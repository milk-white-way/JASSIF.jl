function f  = myfunc(x)
n = length(x);

for i =1 :n
f(i) = cos(x(i)) - x(i)*sin(x(i));
end
f=f';