function [F,violation] = GenerateFunctionalOneRound(sigma)

[~,~,ob(1),mb(1)] = size(sigma);

 [Fprime] = GenerateFunctionalOneRoundSW(sigma);
 F = real(Fprime);
  
 tr1 = zeros(ob(1),mb(1));
 for b1=1:ob(1)
    for y1=1:mb(1)
        
        tr1(b1,y1) = trace(sigma(:,:,b1,y1)*F(:,:,b1,y1));
       
    end
 end
        
 violation  = sum(sum(tr1,1),2);

 
 
 






