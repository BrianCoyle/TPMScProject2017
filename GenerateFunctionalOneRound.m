function [F,v] = GenerateFunctionalOneRound(sigma)

%This function generates the Steering inequality that is maximally violated
%by the assemblage, sigma, using GenerateSteeringFunctionalOneRound. It also
%calculates the observed violation of that steering inequality, v,
%for one round of measurements.

[~,~,ob(1),mb(1)] = size(sigma);

%Calculate steering functional, F
 [Fprime] = GenerateSteeringFunctionalOneRound(sigma);
 F = real(Fprime);
  
 %Calculate steering inequality violation, v
 tr1 = zeros(ob(1),mb(1));
 for b1=1:ob(1)
    for y1=1:mb(1)
        
        tr1(b1,y1) = trace(sigma(:,:,b1,y1)*F(:,:,b1,y1));
       
    end
 end
        
 v  = sum(sum(tr1,1),2);

 
 
 






