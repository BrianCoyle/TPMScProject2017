function [Pg,Hmin] = PlotDualGpThreeRounds(rho1,phi2)

%This function plots the guessing probability, and the number of random
%bits produced, H_{min}, for the given input state rho. The measurement
%angles in the first round are chosen in GenAssemblagesThreeRounds and noisy X
%and projective Z measurements are used. in the second round, X measurement
%is projective, and the Z measurement angle is defined by input phi2. In th
%third round, both X and Z are projective, phi3 = 0, theta3 = 0.

%Generate Assemblages after 3 rounds, sigma{1} = sigma_{b1|y1} and
%sigma{2}= sigma_{b1b2|y1y2}, sigma{3}= sigma_{b1b2b3|y1y2y3}
[sigma,theta,~] = GenAssemblagesThreeRounds(rho1,0,phi2,0,0,0);
r = length(theta);

a = zeros(r,1);
Pg = zeros(r,1);
Hmin = zeros(r,1);

for j=1:r
    sigmatemp{1} = sigma{j,1};
    sigmatemp{2} = sigma{j,2};
    sigmatemp{3} = sigma{j,3};
   
    %Checks if the assemblage for angle j is valid (a = 0) or not (a = 1)
   [a(j)] = ValidAssemblageThreeRounds(sigmatemp);
   
    if a(j) == 0
    
 %Calculates stering inequalities, Ftemp{1}, Ftemp{2} and Ftemp{3}, for 3 measurement rounds
 %and their violations, violationtemp(1), violationtemp(2), and violationtemp(3) for 
 %each measurement angle, theta(j).
[Ftemp,violationtemp] = ...
    GenerateFunctionalThreeRounds(sigmatemp);

%Finds the guessing probability for measurement angle theta(j)
[Pg(j)] = DualGPThreeRounds(sigmatemp,Ftemp,violationtemp);

%Calculates the Min entropy from the guessing probability
Hmin(j) = -log2(Pg(j));
    end
end

%Only plots if a(j)=0  for all j, i.e. no errors with the assemblages
if all(~a) 
%Plot both the guessing probability and min entropy on seperate subplots.

subplot(1,2,1); plot(theta/pi,Pg,'linewidth',1.5); grid on;...
xlabel('\theta_1 (units of \pi)');...
ylabel('Guessing Probability');...
title('Guessing Probability as a function of \theta_1 for three rounds');...
hold on;

subplot(1,2,2); plot(theta/pi,Hmin,'linewidth',1.5);...
grid on;...
xlabel('\theta_1 (units of \pi)');...
ylabel('H_{min} (Number of Random Bits)');...
title('Min Entropy as a function of \theta_1 for three rounds');...
hold on;

else 
    disp('Cannot plot, there is an issue with an assemblage')
end
