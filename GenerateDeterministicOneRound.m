function DeterministicOneRound = GenerateDeterministicOneRound(ob,mb)
%This function calculates the array of possible deterministic strategies
%that Eve can produce. It is based from genSinglePartyArray.m in
%''Steering: a review with focus on semidefinite programming''. This works
%for one round of measurements.

%Total number of deterministic strategies is Ndet for 1 round
Ndet = ob(1)^mb(1); 

DeterministicOneRound = zeros(ob(1),mb(1),Ndet); 

for lam1 = 0:Ndet-1
    lamvar1 = dec2base(lam1,ob(1),mb(1))-'0'; 
    %generates the string of outcomes b1 (for each y1), for the given variable lam1
    
    for y1 = 1:mb(1)
        for b1 = 0:ob(1)-1
            DeterministicOneRound(b1+1,y1,lam1+1) = (lamvar1(y1) == b1);
            % probability = 1 if b1 = lamvar1(y1), 0 otherwise
        end
    end
    
end

end