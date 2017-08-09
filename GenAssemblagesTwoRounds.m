function [sigma,theta,rho] = GenAssemblagesTwoRounds(rho1,phi1,phi2,theta2)
%This function calculates the assemblages for two rounds of measurements, and
%assigns them to a cell array, sigma. Theses assemblages are the result of
%non-projective X measurements, and non-projective Z measurements with angle phi1,
%which is always set to zero and the assemblages are calculated for a range of 
%X measurement angles, theta_1.
%The measurements in the second round are defined by measurment angles,
%phi2 and theta2, but in this case they will be inputted as zero
%The assemblages sigma{1} = sigma_{b1|y1}, sigma{2} = sigma_{b1b2|y1y2},
%the initial round measurement angles, theta, and the post measurment
%states, rho{2} = rho_{b1|y1}, rho{3} = rho_{b1b2|y1y2} are outputted


%Choose initial measurement angle, theta1, for the noisy X measurements
i=1;
thetavar = 0.01;
rho{1,1} = rho1;

%Calculates set of measurement operators for a range of measurement 
%angles, thetavar \in (0,\pi/4)
while thetavar < pi/4
%GenerateMeas generates a set of measurement operators M_{b_1y_1},
%based on a set of Kraus Operators, B_{b_1y_1},
% s.t. B_{b_1y_1} = B_{b_1y_1}^{\dagger}B_{b_1y_1}
[M{i,1},B{i,1},ob(1),mb(1),dB(1)] = GenerateMeas(phi1,thetavar);
dA(1) = length(rho{1,1})/dB(1);
theta(i) = thetavar;
thetavar = thetavar + pi/32;
i=i+1;
end

r = length(theta);

for j=1:r
    
%Generate Measurement set for second round of measurements, with angle
%theta2, phi2
[M{j,2},B{j,2},ob(2),mb(2),dB(2)] = GenerateMeas(phi2,theta2);

dA(2) = dB(2);  
for b1 = 1:ob(1)
    for y1 = 1:mb(1)
     %Generate assemblage for measurements M{j,1} on state rho, 
     %to get assemblage sigma{j,1} for each angle theta(j), after one measurement round.

    sigma{j,1}(:,:,b1,y1)= ...
          PartialTrace(kron(eye(dA(1)),M{j,1}(:,:,b1,y1))*(rho{1,1}),2);
            
       %Generate (unnormalized) post measurement states for measurements M{j,1} on state rho, 
       %for each angle theta(j), after one measurement round.
     rho{j,2}(:,:,b1,y1) = GeneratePostMeasState(B{j,1}(:,:,b1,y1),rho{1,1},0);
     
       for b2=  1:ob(2)
            for y2 = 1:mb(2)
          %Generate assemblage for measurements M{j,2} on state rho{j,2}, 
          %to get assemblage sigma{j,2} for each angle theta(j), after two measurement rounds.

     sigma{j,2}(:,:,b1,b2,y1,y2)= ...
         PartialTrace(kron(eye(dA(2)),M{j,2}(:,:,b2,y2))*(rho{j,2}(:,:,b1,y1)),2);
     
       %Generate (normalized) post measurement states for measurements M{j,2} on state rho{j,2}, 
       %for each angle theta(j), after two measurement rounds.
       
            rho{j,3}(:,:,b1,b2,y1,y2) =...
                 GeneratePostMeasState(B{j,2}(:,:,b2,y2),rho{j,2}(:,:,b1,y1),1); 
 
            end
       end
    end
end

    

end
end
