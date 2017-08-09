function [sigma,theta,rho] = GenAssemblagesThreeRounds(rho1,phi1,phi2,theta2,phi3,theta3)
%This function calculates the assemblages for three rounds of measurements, and
%assigns them to a cell array, sigma. Theses assemblages are the result of
%non-projective X measurements, and non-projective Z measurements with angle phi1,
%which is always set to zero and the assemblages are calculated for a range of 
%X measurement angles, theta_1 in the first round.
%The measurements in the second round are defined by measurement angles,
%phi2 and theta2
%The measurements in the third round are defined by measurement angles,
%phi3 and theta3, but in this case they will be inputted as zero.
%The assemblages sigma{1} = sigma_{b1|y1}, sigma{2} = sigma_{b1b2|y1y2},
%sigma{3} = sigma_{b1b2b3|y1y2y3} the initial round measurement angles, theta, and the post measurment
%states, rho{2} = rho_{b1|y1}, rho{3} = rho_{b1b2|y1y2},  rho{4} = rho_{b1b2b3|y1y2y3} are outputted


%Choose initial measurement angle, theta1, for the non-projective X measurements in
%the first round
i=1;
thetavar =0.05*pi;
rho{1,1} = rho1;

%Calculates set of measurement operators for a range of measurement 
%angles, thetavar \in (0.05\pi,\pi/4)
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
% Generate assemblage for measurements M{1} on state rho, to get assemblage
% sigma{j,1} for each angle theta(j).

for j=1:r
    
%Generate Measurement set for second round of measurements, with angles
%theta2, phi2
[M{j,2},B{j,2},ob(2),mb(2),dB(2)] = GenerateMeas(phi2,theta2);
dA(2) = dB(2);

%Generate Measurement set for third round of measurements, with angles
%theta3, phi3
 [M{j,3},B{j,3},ob(3),mb(3),dB(3)] = GenerateMeas(phi3,theta3);
 dA(3) = dB(3);


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
     
          %Generate (unnormalized) post measurement states for measurements M{j,2} on state rho{j,2}, 
          %for each angle theta(j), after two measurement rounds.
             rho{j,3}(:,:,b1,b2,y1,y2) =...
               GeneratePostMeasState(B{j,2}(:,:,b2,y2),rho{j,2}(:,:,b1,y1),0);
       
              
                  for b3 = 1:ob(3)
                     for y3 = 1:mb(3)
          %Generate assemblage for measurements M{j,3} on state rho{j,3}, 
          %to get assemblage sigma{j,3} for each angle theta(j), after three measurement rounds.
                sigma{j,3}(:,:,b1,b2,b3,y1,y2,y3)= ...
                     PartialTrace(kron(eye(dA(3)),M{j,3}(:,:,b3,y3))*...
                                        (rho{j,3}(:,:,b1,b2,y1,y2)),2);
                                    
          %Generate (normalized) post measurement states for measurements M{j,3} on state rho{j,3}, 
          %for each angle theta(j), after three measurement rounds.
                     rho{j,4}(:,:,b1,b2,b3,y1,y2,y3) =...
                     GeneratePostMeasState(B{j,2}(:,:,b3,y3),rho{j,3}(:,:,b1,b2,y1,y2),1);  
  
                     end
                  end 
            end
         end
    end
end
    
    
end
 
 
    










