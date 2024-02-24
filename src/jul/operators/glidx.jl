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
module GLIDX
   export glidx

   function glidx(ii, jj, mm, nn)
      idx = ( ii-1 )*nn + jj
      if ( ii<1 || ii>mm || jj<1 || jj>nn )
         idx = 0
      end
      return idx
   end
end