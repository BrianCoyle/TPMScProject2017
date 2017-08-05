function [SW, F] = GenerateSteeringFunctionalTwoRounds(sigma)


[dA(1),dA(1),ob(1),mb(1)] = size(sigma{1});
[dA(2),dA(2),~,ob(2),~, mb(2)] = size(sigma{2});


Ndetone = ob(1)^mb(1);
Ndettwo = ob(2)^mb(2);

DeterministicTwoRound = GenerateDeterministicTwoRounds(ob,mb,dA);

% NOTE: Here we use the dual formulation of the steering weight.

cvx_begin sdp quiet
    
    variable F(dA(2),dA(2),ob(1),ob(2),mb(1),mb(2)) hermitian semidefinite
    % members of the steering functional


    maximise 1 - real(sum(reshape(F.*conj(sigma{2}),1,[])))
    % max 1 - sum_b1b2y1y2 trace(F_b1b2|y1y2*sigma_b1b2|y1y2)
    
    subject to
    
    for i = 1:Ndetone
        for j = Ndettwo
            
        
    sum(sum(sum(sum(DeterministicTwoRound(:,:,:,:,:,:,i,j).*F,3),4),5),6) ...
        -eye(dA(2)) == hermitian_semidefinite(dA(2)); 
      
        % standard SW: sum_ax F_ax D(a|x,lam) - eye(dB) >= 0, forall lam
        end
    end

cvx_end

SW = cvx_optval;
    
end