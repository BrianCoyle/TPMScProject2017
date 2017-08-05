function [F,violation] = GenerateFunctionalThreeRounds(sigma)


[dA(1),dA(1),ob(1),mb(1)] = size(sigma{1});
[dA(2),dA(2),~,ob(2), ~ ,mb(2)] = size(sigma{2});
[dA(3),dA(3),~,~,ob(3),~, ~ ,mb(3)] = size(sigma{3});


[Fprime{1}] = GenerateFunctionalOneRoundSW(sigma{1});

F{1} = real(Fprime{1});
tr1 = zeros(ob(1),mb(1));


[~, Fprime{2}] = GenerateSteeringFunctionalTwoRounds(sigma);

 F{2} = real(Fprime{2});
 tr2 = zeros(ob(1),ob(2),mb(1),mb(2));

  F{3} = zeros(size(sigma{3}));
  tr3 = zeros(ob(1),ob(2),ob(3),mb(1),mb(2),mb(3));
 
 alpha = 100;
 
for b1=1:ob(1)
    for y1=1:mb(1)
        
        tr1(b1,y1) = trace(sigma{1}(:,:,b1,y1)*F{1}(:,:,b1,y1));
       
        for b2 = 1:ob(2)
            for y2 = 1:mb(2)
                
        tr2(b1,b2,y1,y2) =...
            trace(sigma{2}(:,:,b1,b2,y1,y2)*F{2}(:,:,b1,b2,y1,y2));
          
        
                for b3 = 1:ob(3)
                  for y3=1:mb(3)
                     if trace(sigma{3}(:,:,b1,b2,b3,y1,y2,y3)) ==0
                             F{3}(:,:,b1,b2,b3,y1,y2,y3) = alpha*(eye(dA(3)) -[0,0;0,0]); 
                     else
                          F{3}(:,:,b1,b2,b3,y1,y2,y3) = alpha*(eye(dA(3)) - ...
                            (sigma{3}(:,:,b1,b2,b3,y1,y2,y3)/trace(sigma{3}(:,:,b1,b2,b3,y1,y2,y3))));
                     end
                        tr3(b1,b2,b3,y1,y2,y3) = ...
                          trace(sigma{3}(:,:,b1,b2,b3,y1,y2,y3)*F{3}(:,:,b1,b2,b3,y1,y2,y3));
                  end
                end
             end
        end
    end
end
        
 violation(1)  = sum(sum(tr1,1),2);

 violation(2)  = sum(sum(sum(sum(tr2,1),2),3),4);

 violation(3) = sum(sum(sum(sum(sum(sum(tr3,1),2),3),4),5),6);

 
 
 
 






