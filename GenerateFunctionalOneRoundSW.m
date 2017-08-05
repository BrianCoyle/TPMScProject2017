function [F] = GenerateFunctionalOneRoundSW(sigma)


% if unspecified, it is assumed the standard steering weight is required.

[dA(1),~,ob(1),mb(1)] = size(sigma);
% dB = dim. of Bob, oa = # outcomes for Alice, ma = # inputs for Alice

Ndet = ob(1)^mb(1);
% number of determinstic probability distributions for Alice
DeterministicOneRound = GenerateDeterministicOneRound(ob,mb);

% generate array containing the single party distributions
% NOTE: Here we use the dual formulation of the steering weight.

cvx_begin sdp quiet
    
    variable F(dA(1),dA(1),ob(1),mb(1)) hermitian semidefinite
    % members of the steering functional

    maximise 1 - real(sum(reshape(F.*conj(sigma),1,[])))
    % max 1 - sum_b1y1 trace(F_b1y1*sigma_b1|y1)
    
    subject to
    
    for i = 1:Ndet
        sum(sum(permute(repmat(DeterministicOneRound(:,:,i),[1,1,dA,dA]),[3,4,1,2]).*F,3),4) ...
             - eye(dA) == hermitian_semidefinite(dA);
        % standard SW: sum_ax F_ax D(a|x,lam) - eye(dB) >= 0, forall lam
        
        % eye(dB) >= 0, forall lam
    end

cvx_end


    
end
