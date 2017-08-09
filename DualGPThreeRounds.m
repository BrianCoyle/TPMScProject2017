function [Pg] = DualGPThreeRounds(sigma,F,v)

%This function calculates the full guessing probability, which is the
%maximum over all outcomes, b_1, b_2, b_3. It takes as inputs
%the assemblage, sigma,  the steering inequality
%violation, v, that is achieved by the steering inequalities F.
%These inputs are passed to DualGPOneRoundThreeOutcome to calculate the
%guessing probability for a single outcome, and then takes the maximum over
%all outcomes to be the guessing probability, which is outputted.

[~,~,ob(1),~] = size(sigma{1});
[~,~,~,ob(2), ~ ,~] = size(sigma{2});
[~,~,~,~,ob(3),~, ~ ,~] = size(sigma{3});
 
outcome = zeros(ob(1),ob(2),ob(3));

for b1 = 1:ob(1)
    for b2 = 1:ob(2)
        for b3  = 1:ob(3)
            outcome(b1,b2,b3) = DualGPThreeRoundsOneOutcome(b1,b2,b3,sigma,F,v);
        end
    end
end

%Maximization over all outcomes.
Pg = max(max(max(outcome)));

end






