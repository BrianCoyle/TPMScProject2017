function [Pg] = DualGPOneRoundOneOutcome(a1,violation1,sigma,F1)

[dA(1),~,ob(1),mb(1)] = size(sigma);


cvx_begin sdp quiet

variable sigma1(dA(1),dA(1),ob(1),mb(1)) hermitian semidefinite

maximise trace(sigma1(:,:,a1,2))

subject to 

for b1 =1:ob(1)
    for y1 = 1:mb(1)
            
      temp1(:,:,b1,y1) =  F1(:,:,b1,y1)*sigma1(:,:,b1,y1);
    
    end
end

 violation1 == trace(sum(sum(temp1,3),4));
        

 for y1 = 1:mb(1)
     
     squeeze(sum(sigma1(:,:,:,y1),3)) == squeeze(sum(sigma1(:,:,:,1),3));     
     trace(squeeze(sum(sigma1(:,:,:,y1),3))) == 1; 
     
 end

  
cvx_end
Pg = cvx_optval;


end






