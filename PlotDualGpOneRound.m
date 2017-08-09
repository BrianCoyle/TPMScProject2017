function [Pg, Hmin] = PlotDualGpOneRound(rho1)
%This function plots the guessing probability, and the number of random
%bits produced, H_{min}, for the given input state rho. The measurement
%angles in the first round are chosen in GenAssemblagesOneRound and noisy X
%and projective Z measurements are used 

[sigma,theta,~] = GenAssemblagesOneRound(rho1);
r = length(theta);

a = zeros(r,1);
Pg = zeros(r,1);
Hmin = zeros(r,1);

for j=1:r
    sigmatemp = sigma{j,1};
    %Checks if the assemblage for angle j is valid (a = 0) or not (a = 1)
    [a(j)] = ValidAssemblageOneRound(sigmatemp);
   
    if a(j) == 0
 %Calculates stering inequality, Ftemp, and violation, violationtemp, for 
 %each measurement angle, theta(j).
[Ftemp,violationtemp] = ...
    GenerateFunctionalOneRound(sigmatemp);

%Finds the guessing probability for measurement angle theta(j)
[Pg(j)] = DualGPOneRound(violationtemp,sigmatemp,Ftemp);

%Calculates the Min entropy from the guessing probability
Hmin(j) = -log2(Pg(j));
    else
        return
    end
end

%Only plots if a(j)=0  for all j, i.e. no errors with the assemblages
if all(~a) 
%Plot both the guessing probability and min entropy on seperate subplots.
subplot(1,2,1); plot(theta/pi,Pg,'linewidth',1.5); grid on;...
xlabel('\theta_1 (units of \pi)');...
ylabel('Guessing Probability');...
title('Guessing Probability as a function of \theta_1 for one round');...
hold on;

subplot(1,2,2); plot(theta/pi,Hmin,'linewidth',1.5);...
grid on;...
xlabel('\theta_1 (units of \pi)');...
ylabel('H_{min} (Number of Random Bits)');...
title('Min Entropy as a function of \theta_1 for one round');...
hold on;

else 
    disp('Cannot plot, there is an issue with an assemblage')
end




