function [Hmin] = PlotDualGpTwoRounds(rho1)

[sigma,theta,~] = GenAssemblagesTwoRounds(rho1,0);
r = length(theta);

sigmatemp{1} = zeros(size(sigma{1,1}));
sigmatemp{2} = zeros(size(sigma{1,2}));

for j=1:r
    sigmatemp{1} = sigma{j,1};
    sigmatemp{2} = sigma{j,2};
    
[Ftemp,violationtemp] = ...
    GenerateFunctionalTwoRounds(sigmatemp);

[Pg(j)] = DualGPTwoRounds(violationtemp,sigmatemp,Ftemp);

Hmin(j) = -log2(Pg(j));
end
subplot(1,2,1); plot(theta/pi,Pg,'linewidth',1.5); grid on;...
xlabel('\theta_1 (units of \pi)');...
ylabel('Guessing Probability');...
title('Guessing Probability as a function of \theta_1 for two rounds');...
hold on;

subplot(1,2,2); plot(theta/pi,Hmin,'linewidth',1.5);...
grid on;...
xlabel('\theta_1 (units of \pi)');...
ylabel('H_{min} (Number of Random Bits)');...
title('Min Entropy as a function of \theta_1 for two rounds');...
hold on;







