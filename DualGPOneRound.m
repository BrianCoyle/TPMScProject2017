function [Pg] = DualGPOneRound(v,sigma,F1)
 %This function calculates the full guessing probability, which is the
 %maximum over all outcomes, b_1. It takes as inputs
%the assemblage, sigma,  the steering inequality
%violation, vthat is achieved by the steering inequality F1.
%These inputs are passed to DualGPOneRoundOneOutcome to calculate the
%guessing probability for a single outcome, and then takes the maximum over
%all outcomes to be the guessing probability, which is outputted.
[~,~,ob(1),~] = size(sigma);

outcome = zeros(ob(1),1);

for b1 = 1:ob(1)
         
outcome(b1)= DualGPOneRoundOneOutcome(b1,v,sigma,F1);

end
%Maximization over all outcomes.
Pg = max(outcome);

end