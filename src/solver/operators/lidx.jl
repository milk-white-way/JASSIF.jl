"""
# module glidx

- Julia version: 1.8.5
- Author: tamnt
- Date: 2023-04-03

# Examples

```jldoctest
julia>
```
"""
module LIDX
   export lidx

   function lidx(IM, M, N)
      i::UInt128 = floor(IM/N)
      j::UInt128 = IM - ( N*i )
      if (j != 0)
         i += 1
      else
         j = N
      end
      return i, j
   end
end