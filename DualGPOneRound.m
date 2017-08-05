function [Pg] = DualGPOneRound(violation1,sigma,F1)

[~,~,ob(1),~] = size(sigma);

outcome = zeros(ob(1),1);

for b1 = 1:ob(1)
         
outcome(b1)= DualGPOneRoundOneOutcome(b1,violation1,sigma,F1);

end

Pg = max(outcome);

