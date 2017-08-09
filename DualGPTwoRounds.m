function [Pg] = DualGPTwoRounds(v,sigma,F)

 %This function calculates the full guessing probability, which is the
 %maximum over all outcomes, b_1, b_2. It takes as inputs
%the assemblage, sigma,  the steering inequality
%violation, v, that is achieved by the steering inequalities F.
%These inputs are passed to DualGPOneRoundTwoOutcome to calculate the
%guessing probability for a single outcome, and then takes the maximum over
%all outcomes to be the guessing probability, which is outputted.

[~,~,ob(1),~] = size(sigma{1});
[~,~,~,ob(2), ~ ,~] = size(sigma{2});

outcome = zeros(ob(1),ob(2));

for b1 = 1:ob(1)
    for b2 =1:ob(2)
        
outcome(b1,b2)= ...
    DualGPTwoRoundsOneOutcome(b1,b2,v,sigma,F);

    end
end

%Maximization over all outcomes.
Pg = max(max(outcome));

end

