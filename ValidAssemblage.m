  function [a] = ValidAssemblage(sigma)

%Numerical Tolerance
tol = 1e-10; 

    [dA(1),dA(1),ob(1),mb(1)] = size(sigma{1});
    [dA(2),dA(2),~,ob(2),~,mb(2)] = size(sigma{2});



for y1 = 2:mb(1)
        for b1 = 1:ob(1)
              for y2 = 1:mb(2)
                  for b2 = 1:ob(2)
                      
                   if ~all(eig(sigma{1}(:,:,b1,y1))>= -tol)
                        
                      disp(['assemblage ',num2str(1),' not positive semidefinite'])
                      a = 1;
                         return
                   else
                       sig1 = 0;
                   end 
                  
                   if ~all(eig(sigma{2}(:,:,b1,b2,y1,y2))>= -tol)
                         disp(['assemblage ',num2str(2),' not positive semidefinite'])
                         a = 1;
                         return
                   else
                       sig2 = 0;
                   end
                  
                  end
              end
        end
  end
  for n1 = 1:dA(1)
     for n2 = 1:dA(1)
        for a1 = 1:dA(2)
          for a2 = 1:dA(2)
            for y1 = 1:mb(1)
                for b1 = 1:ob(1)
                   for y2 = 1:mb(2)
                       for b2 = 1:ob(2)
                           
                      
if ~(abs(sum(sigma{1}(n1,n2,:,2),3) - sum(sigma{1}(n1,n2,:,1),3)) <=tol)
      
           disp(['assemblage ',num2str(1),' is signalling'])
         a = 1;
   return
else 
    sig3 = 0;
end

if (~(abs(squeeze(sum(sum(sigma{2}(a1,a2,:,:,2,2),3),4))-squeeze(sum(sum(sigma{2}(a1,a2,:,:,2,1),3),4))) <= tol) ||... 
        ~(abs(squeeze(sum(sum(sigma{2}(a1,a2,:,:,1,2),3),4))-squeeze(sum(sum(sigma{2}(a1,a2,:,:,1,1),3),4))) <= tol))
      disp(['assemblage ',num2str(2),' is signalling'])
        a = 1;
           return
else
    sig4 = 0;
end
                       end
                   end
                end
            end
          end
        end
     end
  end
  
  %If all requirements are met for the assemblage to be valid, output a =
  %0, If not output a = 1.
  if sig1 == 0 && sig2 == 0 && sig3 == 0 && sig4 == 0
      a = 0;
  else 
      a = 1;
  end

  
  end