function [Pg] = DualGPOneRoundOneOutcome(a1,violation1,sigma,F1)

%This function calculates the guessing probability for a particular 
%outcome in one round of measurements, a_1, for measurement choice y_1=2,
%i.e. a measurement using non-projective X operators. It takes as inputs
%the assemblage, sigma, the outcome choice, a1, the steering inequality
%violation, violation1, that is achieved by the steering inequality F1.
%The output is the guessing probability for a single outcome.
[dA(1),~,ob(1),mb(1)] = size(sigma);

%Start cvx to solve optimisation problem
cvx_begin sdp quiet

%Declare cvx variable which is the assemblages that Eve could create
variable sigma1e(dA(1),dA(1),ob(1),mb(1)) hermitian semidefinite


maximise trace(sigma1e(:,:,a1,2))

subject to 

%Enforcing the constraint that the assemblages that Eve can create must
%satisfy the observed steering inequality violation.
for b1 =1:ob(1)
    for y1 = 1:mb(1)
            
      temp1(:,:,b1,y1) =  F1(:,:,b1,y1)*sigma1e(:,:,b1,y1);
    
    end
end

 violation1 == trace(sum(sum(temp1,3),4));
        

 for y1 = 1:mb(1)
     %No-signalling constraints on the assemblages
     squeeze(sum(sigma1e(:,:,:,y1),3)) == squeeze(sum(sigma1e(:,:,:,1),3));     
     %Constraint that Eve's assembalges must create a valid quantum state.
     trace(squeeze(sum(sigma1e(:,:,:,y1),3))) == 1; 
     
 end

  
cvx_end
Pg = cvx_optval;


end






