  function [a] = ValidAssemblageOneRound(sigma)

 %This function checks if the assemblage created after the first
 %measurement round is valid, i.e. if it is non-signalling, and the
 %elements are positive, semidefinite matrices. The output a is a flag, if
 %a = 0, the assemblage is valid, if a = 1, the assemblage is not valid.

%Numerical Tolerance
tol = 1e-10; 

    [dA(1),dA(1),ob(1),mb(1)] = size(sigma);
   
    
  for i1 = 1:dA(1)
     for i2 = 1:dA(1)
       for y1 = 1:mb(1)
                   %Checks that the assemblage is non-signalling to within the tolerance                   
                   if ~(abs(sum(sigma(i1,i2,:,y1),3) - sum(sigma(i1,i2,:,y1),3)) <=tol)
                     disp(['assemblage ',num2str(1),' is signalling'])
                     a = 1;
                   return
                   else 
                     signal1 = 0; 
                   end
            
          for b1 = 1:ob(1)
              %Checks if the assemblage elements are positive semidefinite
                         if ~all(eig(sigma(:,:,b1,y1))>= -tol)
                                 disp(['assemblage ',num2str(1),' not positive semidefinite'])
                                 a = 1;
                         return
                         else
                           signal2 = 0;
                         end 
            
          end
       end
     end
  end
        
  %If all requirements are met for the assemblage to be valid, output a =
  %0, If not output a = 1.
  if signal1 == 0 && signal2 == 0 
      a = 0;
  else 
      a = 1;
  end

  end