function [Pg] = DualGPThreeRounds(sigma,F,violation)

%This function caluculates and plots the GP which is the maximum of the
%single outcome GP, over all outcomes b1, b2, b3.
[~,~,ob(1),~] = size(sigma{1});
[~,~,~,ob(2), ~ ,~] = size(sigma{2});
[~,~,~,~,ob(3),~, ~ ,~] = size(sigma{3});
 
outcome = zeros(ob(1),ob(2),ob(3));

for b1 = 1:ob(1)
    for b2 = 1:ob(2)
        for b3  = 1:ob(3)
            outcome(b1,b2,b3) = DualGPThreeRoundsOneOutcome(b1,b2,b3,sigma,F,violation);
        end
    end
end

Pg = max(max(max(outcome)));

end






