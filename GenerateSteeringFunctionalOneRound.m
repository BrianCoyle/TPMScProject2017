function [F] = GenerateSteeringFunctionalOneRound(sigma)

%This function finds the steering weight of an assemblage via the dual
%problem, which also calculates the steering inequality, F which is maximally
%violated by the assemblage sigma. It is based from steeringWeight.m in
%''Steering: a review with focus on semidefinite programming''. This works
%for one round of measurements

[dA(1),~,ob(1),mb(1)] = size(sigma);


Ndet = ob(1)^mb(1);
% Number of determinstic probability distributions for Bob
 
%GenerateDeterministicOneRound generates an array of the deterministic
%probability distributions that Eve can create
DeterministicOneRound = GenerateDeterministicOneRound(ob,mb);

%Solve SDP to find steering weight and Dual variables, F
cvx_begin sdp quiet
    
    variable F(dA(1),dA(1),ob(1),mb(1)) hermitian semidefinite
    % members of the steering functional

    maximise 1 - real(sum(reshape(F.*conj(sigma),1,[])))
    % max 1 - sum_b1y1 trace(F_b1y1*sigma_b1|y1)
    
    subject to
    
    %Constraint in SDP that checks given assemblage, sigma, against all
    %possible deterministic strategies.
    for i = 1:Ndet
        sum(sum(permute(repmat(DeterministicOneRound(:,:,i),[1,1,dA,dA]),[3,4,1,2]).*F,3),4) ...
             - eye(dA) == hermitian_semidefinite(dA);
      
    end

cvx_end
end
