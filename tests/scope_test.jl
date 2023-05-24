

# Test that the scope of a variable is correct
a = zeros(Float64, 1, 7)
function myfunc(a)
   let b = deepcopy(a)
       c = deepcopy(a)
      for ii = 1:7
         b[ii] = 1
         c[ii] = 0
      end
   println(b)
   return b,c
   end
end

myfunc(a)
println(b)
println(a)