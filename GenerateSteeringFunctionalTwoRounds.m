function [F] = GenerateSteeringFunctionalTwoRounds(sigma)


%This function finds the steering weight of an assemblage via the dual
%problem, which also calculates the steering inequality, F which is maximally
%violated by the assemblage sigma{2}. It is based from steeringWeight.m in
%''Steering: a review with focus on semidefinite programming'', but adapted for two rounds of measurements


[dA(1),dA(1),ob(1),mb(1)] = size(sigma{1});
[dA(2),dA(2),~,ob(2),~, mb(2)] = size(sigma{2});

% Number of determinstic probability distributions is Ndetone x Ndettwo
Ndetone = ob(1)^mb(1);
Ndettwo = ob(2)^mb(2);

%GenerateDeterministicTwoRounds generates an array of the deterministic
%probability distributions that Eve can create
DeterministicTwoRound = GenerateDeterministicTwoRounds(ob,mb,dA);

%Solve SDP to find steering weight and Dual variables, F
cvx_begin sdp quiet
    
    variable F(dA(2),dA(2),ob(1),ob(2),mb(1),mb(2)) hermitian semidefinite
    % members of the steering functional


    maximise 1 - real(sum(reshape(F.*conj(sigma{2}),1,[])))
    % max 1 - sum_b1b2y1y2 trace(F_b1b2|y1y2*sigma_b1b2|y1y2)
    
    subject to
    
    %Constraint in SDP that checks given assemblage, sigma, against all
    %possible deterministic strategies.
    for i = 1:Ndetone
        for j = Ndettwo
               
    sum(sum(sum(sum(DeterministicTwoRound(:,:,:,:,:,:,i,j).*F,3),4),5),6) ...
        -eye(dA(2)) == hermitian_semidefinite(dA(2)); 
      
     
        end
    end

cvx_end
end