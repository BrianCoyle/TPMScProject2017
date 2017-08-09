 function [a] = ValidAssemblageTwoRounds(sigma)

 %This function checks if the assemblage created after the three
 %measurement rounds is valid, i.e. if it is non-signalling, and the
 %elements are positive, semidefinite matrices. The output a is a flag, if
 %a = 0, the assemblage is valid, if a = 1, the assemblage is not valid.

%Numerical Tolerance
tol = 1e-10; 

    [dA(1),dA(1),ob(1),mb(1)] = size(sigma{1});
    [dA(2),dA(2),~,ob(2),~,mb(2)] = size(sigma{2});


        
  for i1 = 1:dA(1)
     for i2 = 1:dA(1)
       for y1 = 1:mb(1)
           
            %Checks that the assemblage, sigma_{b1|y1} is non-signalling to within the tolerance                   
           if ~(abs(sum(sigma{1}(i1,i2,:,y1),3) - sum(sigma{1}(i1,i2,:,y1),3)) <=tol)
            disp(['assemblage ',num2str(1),' is signalling'])
           a = 1;
           return
           else 
             signal1 = 0; 
           end
           
          for b1 = 1:ob(1)
              
           %Checks if the assemblage elements are positive semidefinite
                         if ~all(eig(sigma{1}(:,:,b1,y1))>= -tol)
                                 disp(['assemblage ',num2str(1),' not positive semidefinite'])
                                 a = 1;
                         return
                         else
                           signal2 = 0;
                         end 
                                    
            for j1 = 1:dA(2)
                for j2 = 1:dA(2)
                   for y2 = 1:mb(2)
                       
                 %Checks that the assemblage, sigma_{b1b2|y1y2} is non-signalling to within the tolerance
                               if ~(abs(squeeze(sum(sum(sigma{2}(j1,j2,:,:,y1,y2),3),4))-...
                                        squeeze(sum(sum(sigma{2}(j1,j2,:,:,y1,y2),3),4))) <= tol) 
                                     disp(['assemblage ',num2str(2),' is signalling'])
                                       a = 1;
                               return
                               else
                                 signal3 = 0;
                               end
                                  
                       for b2 = 1:ob(2)
                           
                %Checks if the assemblage elements are positive semidefinite
                              if ~all(eig(sigma{2}(:,:,b1,b2,y1,y2))>= -tol)
                                 disp(['assemblage ',num2str(2),' not positive semidefinite'])
                                  a = 1;
                                return
                              else
                                  signal4 = 0;
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
  if signal1 == 0 && signal2 == 0 && signal3 == 0 && signal4 == 0
      a = 0;
  else 
      a = 1;
  end

  end