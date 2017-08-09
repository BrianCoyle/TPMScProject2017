function [sigma,theta,rho] = GenAssemblagesOneRound(rho1)
%This function calculates the assemblages for 1 round of measurements, and
%assigns them to a cell array, sigma. Theses assemblages are the result of
%non-projective X measurements, and projective Z measurements, and the
%assemblages are calculated for a range of X measurement angles, theta_1
%The assemblages sigma{1} = sigma_{b1|y1},
%the first round measurement angles, theta, and the post measurement
%states, rho{2} = rho_{b1|y1} are outputted
%Choose initial measurement angle, theta1, for the noisy X measurements

i=1;
thetavar = 0;
rho{1,1} = rho1;

%Calculates set of measurement operators for a range of measurement 
%angles, thetavar \in [0,\pi/4]
while thetavar <= pi/4
%GenerateMeas generates a set of measurement operators M_{b_1y_1},
%based on a set of Kraus Operators, B_{b_1y_1},
% s.t. B_{b_1y_1} = B_{b_1y_1}^{\dagger}B_{b_1y_1}
[M{i,1},B{i,1},ob(1),mb(1),dB(1)] = GenerateMeas(0,thetavar);
dA(1) = length(rho{1,1})/dB(1);
theta(i) = thetavar;
thetavar = thetavar + pi/32;
i=i+1;
end

r = length(theta);

for j=1:r
    
for b1 = 1:ob(1)
    for y1 = 1:mb(1)
     
%Generate assemblage for measurements M{j,1} on state rho, 
%to get assemblage sigma{j,1} for each angle theta(j).

 sigma{j,1}(:,:,b1,y1)= ...
     PartialTrace(kron(eye(dA(1)),M{j,1}(:,:,b1,y1))*(rho{1,1}),2);
 %Generates (normalized) post measurement states for measurements M{j,1} on state rho, 
 %for each angle theta(j).

 rho{j,2}(:,:,b1,y1) = GeneratePostMeasState(B{j,1}(:,:,b1,y1),rho{1,1},1);
    end
end


end

