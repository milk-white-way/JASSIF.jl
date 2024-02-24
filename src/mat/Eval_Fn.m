function Fpx1_i = Eval_Fp(i,M, Ucont,Ucat_2, Ucat_1,Ucat_0,Ucat_m_1,d_Cont,d_Cat)

        coef = 1/8;

        ucon = Ucont /  2;
         
        up = ucon + abs(ucon);
        um = ucon - abs(ucon);
        
        Fpx1_i = 0;
         
   if (d_Cat == 1) 
       % Derivative w.r.t U_cat  
       if (i > 1 && i < M-1)
            Fpx1_i =  um * ( coef * 3) +  up * ( coef * (-2)  + 1);                          
       else
          if (i == 1)
              Fpx1_i =um * ( coef * (3)) +  up * ( coef * ( -1-2 ) + 1);                                                      
          else
                 if (i == M-1)
                    Fpx1_i =  um * ( coef * (3)) + up * ( coef * (-2) + 1);                                             
                 end % End of i == M-1
          end % End of i == 1
       end % End of internal nodes  

   end
       %Derivative w.r.t Ucont_x
       if ucon > 0            
            ap = 0.5;
            am = 0;
       else
            ap = 0;
            um = -0.5;
       end
       
       
  if (d_Cont == 1)     
       if (i > 1 && i < M-1)
            Fpx1_i = Fpx1_i +  am * ( coef * ( - Ucat_2 - 2 * Ucat_1 + 3 * Ucat_0  ) + Ucat_1) +  ap * ( coef * ( - Ucat_m_1 - 2 * Ucat_0   + 3 * Ucat_1) + Ucat_0);                          
       else
          if (i == 1)
              Fpx1_i  = Fpx1_i + am * ( coef * ( - Ucat_2 - 2 * Ucat_1 + 3 * Ucat_0  ) + Ucat_1) +  ap * ( coef * ( - Ucat_0   - 2 * Ucat_0   + 3 * Ucat_1) + Ucat_0);                                                      
          else
                 if (i == M-1)
                    Fpx1_i =  Fpx1_i + am * ( coef * ( - Ucat_1 - 2 * Ucat_1 + 3 * Ucat_0  ) + Ucat_1) + ap * ( coef * ( - Ucat_m_1 - 2 * Ucat_0   + 3 * Ucat_1) + Ucat_0);                                             
                 end % End of i == M-1
          end % End of i == 1
      end % End of internal nodes       
  end
  
  
      % Zero-out the boundary      
      if i==M
      Fpx1_i = 0;    
      end
      
      