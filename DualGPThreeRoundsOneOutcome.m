function outcome = DualGPThreeRoundsOneOutcome(a1,a2,a3,sigma,F,violation)

%This function calculates the GP for a particular choice of b1 = a1, b2=a2,
%b3=a3. To get the full guessing probability, all possible values of
%b1, b2, b3 must be maximised over.



%Input observed assemblages, sigma_{b1|y1}, sigma_{b1b2|y1y2}, sigma_{b1b2b3|y1y2y3}
[dA(1),dA(1),ob(1),mb(1)] = size(sigma{1});
[dA(2),dA(2),~,ob(2), ~ ,mb(2)] = size(sigma{2});
[dA(3),dA(3),~,~,ob(3),~, ~ ,mb(3)] = size(sigma{3});
 

cvx_begin sdp quiet
%Declare cvx variables: sigma^E_{b1|y1}, sigma^E_{b1b2|y1y2},
%sigma^E_{b1b2b3|y1y2y3}, the assemblages the Eve can create to maximise
%her GP, possible with some LHS model
variable sigma1(dA(1),dA(1),ob(1),mb(1)) hermitian semidefinite
variable sigma2(dA(2),dA(2),ob(1),ob(2),mb(1),mb(2)) hermitian semidefinite
variable sigma3(dA(3),dA(3),ob(1),ob(2),ob(3),mb(1),mb(2),mb(3)) hermitian semidefinite

%Maximise the trace of sigma^E_{b1b2b3|y1y2y3} for particular
%b1=a1, b2=a2, b3=a3, and given y1 = y1*,  y2 = y2*,  y3 = y3*
maximise trace(sigma3(:,:,a1,a2,a3,2,1,2))

subject to 

%Causality constraints on the assemblages that Eve can create, sum_{b2}
%sigma^E_{b1b2|y1y2} = sigma^E_{b1|y1} for all y2,... etc.
for y2 = 1:mb(2)
    sigma1 == squeeze(sum(sigma2(:,:,:,:,:,y2),4));
end
for y3 = 1:mb(3)
    sigma2 == squeeze(sum(sigma3(:,:,:,:,:,:,:,y3),5));
end

%Define temporary variables temp = F*sigma, for each round. All values in
%trace(temp) will be summed over for all values of b and y.

for b1 =1:ob(1)
    for y1 = 1:mb(1)
         temp1(:,:,b1,y1) =  F{1}(:,:,b1,y1)*sigma1(:,:,b1,y1);

            for b2 =1:ob(2)
                for y2 = 1:mb(2)
                    temp2(:,:,b1,b2,y1,y2) = ...
                         F{2}(:,:,b1,b2,y1,y2)*sigma2(:,:,b1,b2,y1,y2);
                    
                    for b3 = 1:ob(3)
                         for y3 = 1:mb(3)
                          temp3(:,:,b1,b2,b3,y1,y2,y3) =...
                         F{3}(:,:,b1,b2,b3,y1,y2,y3)*sigma3(:,:,b1,b2,b3,y1,y2,y3);
                         end
                    end
                    
                end
            end
            
    end
end
  
%Violation constraints: the observed Steering inequality violations for each
%measurement round. The assemblages Eve creates must satisfy these observed
%violations
 violation(1) == trace(sum(sum(temp1,3),4));
        
 violation(2) == trace(sum(sum(sum(sum(temp2,3),4),5),6));
  
 violation(3) == trace(sum(sum(sum(sum(sum(sum(temp3,3),4),5),6),7),8));
 
 %No signalling constraints on the assemblages Eve can create, i.e. the
 %assemblages cannot signal faster than the speed of light, and also they
 %must sum over all values of b and y to give a valid quantum state on
 %Alice's side
 for y1 = 1:mb(1)
     squeeze(sum(sigma1(:,:,:,y1),3)) == squeeze(sum(sigma1(:,:,:,1),3));     
     trace(squeeze(sum(sigma1(:,:,:,y1),3))) == 1; 
     
     for y2=1:mb(2)
  trace(squeeze(sum(sum(sigma2(:,:,:,:,y1,y2),3),4))) == 1 ;
  
  squeeze(sum(sum(sigma2(:,:,:,:,y1,y2),3),4)) ==...
               squeeze(sum(sum(sigma2(:,:,:,:,1,1),3),4));
           
           for y3 = 1:mb(3)
   trace(squeeze(sum(sum(sum(sigma3(:,:,:,:,:,y1,y2,y3),3),4),5))) == 1;
   
   squeeze(sum(sum(sum(sigma3(:,:,:,:,:,y1,y2,y3),3),4),5)) ==...
               squeeze(sum(sum(sum(sigma3(:,:,:,:,:,1,1,1),3),4),5));
           end
     end
 end

  
cvx_end
 outcome = cvx_optval;
