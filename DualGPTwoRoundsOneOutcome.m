function [Pg] = DualGPTwoRoundsOneOutcome(a1,a2,violation,sigma,F)

[dA(1),~,ob(1),mb(1)] = size(sigma{1});
[dA(2),~,~,ob(2), ~ ,mb(2)] = size(sigma{2});

cvx_begin sdp quiet

variable sigma1(dA(1),dA(1),ob(1),mb(1)) hermitian semidefinite
variable sigma2(dA(2),dA(2),ob(1),ob(2),mb(1),mb(2)) hermitian semidefinite

maximise trace(sigma2(:,:,a1,a2,2,1))

subject to 


%Causality constraints
for y2 = 1:ob(2)
    sigma1 == squeeze(sum(sigma2(:,:,:,:,:,y2),4));
end

for b1 =1:ob(1)
    for y1 = 1:mb(1)
      temp1(:,:,b1,y1) =  F{1}(:,:,b1,y1)*sigma1(:,:,b1,y1);
      
        for b2 =1:ob(2)
             for y2 = 1:mb(2)
                 
     
      temp2(:,:,b1,b2,y1,y2) =  F{2}(:,:,b1,b2,y1,y2)*sigma2(:,:,b1,b2,y1,y2);
      
             end
        end
    end
end

 violation(1) == trace(sum(sum(temp1,3),4));
        
 violation(2) == trace(sum(sum(sum(sum(temp2,3),4),5),6));
       
 
 for y1 = 1:mb(1)
     squeeze(sum(sigma1(:,:,:,y1),3)) == squeeze(sum(sigma1(:,:,:,1),3));     
     trace(squeeze(sum(sigma1(:,:,:,y1),3))) == 1; 
     
     for y2=1:mb(2)
  trace(squeeze(sum(sum(sigma2(:,:,:,:,y1,y2),3),4))) == 1 ;
  squeeze(sum(sum(sigma2(:,:,:,:,y1,y2),3),4)) == squeeze(sum(sum(sigma2(:,:,:,:,1,1),3),4));

     end
 end

  
cvx_end
Pg = cvx_optval;


end






