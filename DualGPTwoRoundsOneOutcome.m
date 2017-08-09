function [Pg] = DualGPTwoRoundsOneOutcome(a1,a2,v,sigma,F)

%This function calculates the GP for a particular choice of b1 = a1, b2=a2,
%b3=a3. To get the full guessing probability, all possible values of
%b1, b2, b3 must be maximised over. Measurement choices y_1=2, y_2 =1,
%i.e. a measurement using non-projective X operators in the first round,
%and projective Z in the second round, are used. It takes as inputs
%the assemblage, sigma, the outcome choices, a1, a2 the steering inequality
%violations, v, that is achieved by the steering inequalites F{1}, F{2}.
%The output is the guessing probability for a single outcome.


%Input observed assemblages, sigma_{b1|y1}, sigma_{b1b2|y1y2}, sigma_{b1b2b3|y1y2y3}
[dA(1),~,ob(1),mb(1)] = size(sigma{1});
[dA(2),~,~,ob(2), ~ ,mb(2)] = size(sigma{2});

%Start cvx to solve optimisation problem
cvx_begin sdp quiet

%Declare cvx variables: sigma^E_{b1|y1}, sigma^E_{b1b2|y1y2},
% the assemblages the Eve can create to maximise
%her GP, possible with some LHS model
variable sigma1e(dA(1),dA(1),ob(1),mb(1)) hermitian semidefinite
variable sigma2e(dA(2),dA(2),ob(1),ob(2),mb(1),mb(2)) hermitian semidefinite

%Maximise the trace of sigma^E_{b1b2|y1y2} for particular
%b1=a1, b2=a2, b3=a3, and given y1 = y1*,  y2 = y2*,  
maximise trace(sigma2(:,:,a1,a2,2,1))

subject to 

%Causality constraints on the assemblages that Eve can create, sum_{b2}
%sigma^E_{b1b2|y1y2} = sigma^E_{b1|y1} for all y2,... etc.
for y2 = 1:ob(2)
    sigma1e == squeeze(sum(sigma2e(:,:,:,:,:,y2),4));
end

%Define temporary variables temp = F*sigma, for each round. All values in
%trace(temp) will be summed over for all values of b and y.
for b1 =1:ob(1)
    for y1 = 1:mb(1)
      temp1(:,:,b1,y1) =  F{1}(:,:,b1,y1)*sigma1e(:,:,b1,y1);
      
        for b2 =1:ob(2)
             for y2 = 1:mb(2)
                 
      temp2(:,:,b1,b2,y1,y2) =  F{2}(:,:,b1,b2,y1,y2)*sigma2e(:,:,b1,b2,y1,y2);
      
             end
        end
    end
end
%Violation constraints: the observed Steering inequality violations for each
%measurement round. The assemblages Eve creates must satisfy these observed
%violations
 v(1) == trace(sum(sum(temp1,3),4));
        
 v(2) == trace(sum(sum(sum(sum(temp2,3),4),5),6));
       
 %No signalling constraints on the assemblages Eve can create, i.e. the
 %assemblages cannot signal faster than the speed of light, and also they
 %must sum over all values of b and y to give a valid quantum state on
 %Alice's side
 for y1 = 1:mb(1)
     squeeze(sum(sigma1e(:,:,:,y1),3)) == squeeze(sum(sigma1e(:,:,:,1),3));     
     trace(squeeze(sum(sigma1e(:,:,:,y1),3))) == 1; 
     
     for y2=1:mb(2)
  trace(squeeze(sum(sum(sigma2e(:,:,:,:,y1,y2),3),4))) == 1 ;
  squeeze(sum(sum(sigma2e(:,:,:,:,y1,y2),3),4)) == squeeze(sum(sum(sigma2e(:,:,:,:,1,1),3),4));

     end
 end

  
cvx_end
Pg = cvx_optval;


end






