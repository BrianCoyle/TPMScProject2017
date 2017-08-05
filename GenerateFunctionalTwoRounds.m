function [F,violation] = GenerateFunctionalTwoRounds(sigma)

[dA(1),dA(1),ob(1),mb(1)] = size(sigma{1});
[dA(2),dA(2),~,ob(2),~,mb(2)] = size(sigma{2});

[Fprime{1}] = GenerateFunctionalOneRoundSW(sigma{1});

F{1} = real(Fprime{1});

tr1 = zeros(ob(1),mb(1));

for b1=1:ob(1)
    for y1=1:mb(1)
        
        tr1(b1,y1) = trace(sigma{1}(:,:,b1,y1)*F{1}(:,:,b1,y1));
       
    end
end

violation(1)  = sum(sum(tr1,1),2);

F{2} = zeros(size(sigma{2}));
alpha = 100;
tr2 = zeros(ob(1),ob(2),mb(1),mb(2));

 for y1 = 1:mb(1)
     for y2 = 1:mb(2)
         for b1 = 1:ob(1)
             for b2 = 1:ob(2)
                 if trace(sigma{2}(:,:,b1,b2,y1,y2)) ==0
                   F{2}(:,:,b1,b2,y1,y2) = alpha*(eye(dA(2)) -[0,0;0,0]); 
                 else
    F{2}(:,:,b1,b2,y1,y2) = alpha*(eye(dA(2)) - ...
                (sigma{2}(:,:,b1,b2,y1,y2)/trace(sigma{2}(:,:,b1,b2,y1,y2))));
                 end
    tr2(b1,b2,y1,y2) = trace(sigma{2}(:,:,b1,b2,y1,y2)*F{2}(:,:,b1,b2,y1,y2));
    
             end
         end
     end
 end
 violation(2)  = sum(sum(sum(sum(tr2,1),2),3),4);

 
 
 






