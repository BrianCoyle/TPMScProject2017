function [Pg] = DualGPTwoRounds(violation,sigma,F)

[~,~,ob(1),~] = size(sigma{1});
[~,~,~,ob(2), ~ ,~] = size(sigma{2});

outcome = zeros(ob(1),ob(2));

for b1 = 1:ob(1)
    for b2 =1:ob(2)
        
outcome(b1,b2)= ...
    DualGPTwoRoundsOneOutcome(b1,b2,violation,sigma,F);

    end
end

Pg = max(max(outcome));

