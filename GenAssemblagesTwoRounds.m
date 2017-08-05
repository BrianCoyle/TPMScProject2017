function [sigma,theta,rho] = GenAssemblagesTwoRounds(rho1,theta2)
%This function calculates the assemblages for 2 rounds of measurements, and
%assigns them to a cell array, sigma. The first part calculates the
%assemblage for only 1 measurement angle (theta1 needs to be specified),
%and the second calculates it for a range of initial angles. Only the
%assemblage elements for the second round are calculated if the state
%passes through X measurements in the first round.

%---------------------------------------------------------------------------
%Uncomment the following to use only one measurement angle for the initial
%measurement set
%---------------------------------------------------------------------------

% 
% %Choose initial measurement angle, theta1, for the noisy X measurements
% 
% rho{1} = rho1;
% 
% %GenerateMeas generates a measurement Mb1y1, based on a set of Kraus
% %Operators, Bb1y1, s.t. Mb1y1 = Bb1y1^{\dagger}Bb1y1
% [M{1},B{1},ob(1),mb(1),dB(1)] = GenerateMeas(theta1);
% dA(1) = length(rho{1})/dB(1);
% 
% %Generate assemblage for measurements M{1} on state rho, to get assemblage
% %sigma{j,1} for each angle theta(j).
% 
% 
% for b1 = 1:ob(1)
%     for y1 = 1:mb(1)
%      
%  sigma{1}(:,:,b1,y1)= ...
%      PartialTrace(kron(eye(dA),M{1}(:,:,b1,y1))*(rho{1}),2);
%             
% 
%     end
% end
% %Calculate guessing probability for round 1, Pg1, and associated Steering
% %inequality, F1.
% %[Pg1(j), F1(:,:,:,:,j)] = localSteeringGuessProb(sigma_b1y1(:,:,:,:,j));
% 
% %Generate (unnormalised) Post measurement states after first round of measurements,
% %rho^2_{b1|y1}
% rho{2} = post_meas(B{1},rho{1},0);
% 
% %Generate Measurement set for second round of measurments, with angle
% %theta2e
% [M{2},B{2},ob(2),mb(2),dB(2)] = GenerateMeas(theta2);
% 
% dA(2) = length(rho{2})/dB(2);
% %Generate assemblage (sigma_{b_2|y_2,b_1|y_1=0}) for new measurement set if
% %Bob measured in X basis in first round (y_1 = 0)
% for b1 = 1:ob(1)
%     for b2 = 1:ob(2)
%         for y1=1:mb(1)
%          for y2 = 1:mb(2)
%           
%  sigma{2}(:,:,b1,b2,y1,y2)= ...
%      PartialTrace(Tensor(eye(dA),M{2}(:,:,b2,y2))*(rho{2}(:,:,b1,y1)),2);
%             
%             
%          end
%         end
%     end
% end
%     
% 
%  rho{3}(:,:,1,:,1,:) = post_meas(B{2},rho{2}(:,:,1,1),1);
%  rho{3}(:,:,2,:,1,:) = post_meas(B{2},rho{2}(:,:,2,1),1);
%  rho{3}(:,:,1,:,2,:) = post_meas(B{2},rho{2}(:,:,1,2),1);
%  rho{3}(:,:,2,:,2,:) = post_meas(B{2},rho{2}(:,:,2,2),1);
% % ---------------------------------------------------------------------------
% % Uncomment the following to use a range of angles for the initial
% % measurement set
% % ---------------------------------------------------------------------------

%Choose initial measurement angle, theta1, for the noisy X measurements
i=1;
theta1 = 0.01;
rho{1,1} = rho1;
while theta1 < pi/4
%GenerateMeas generates a measurement Mb1y1, based on a set of Kraus
%Operators, Bb1y1, s.t. Mb1y1 = Bb1y1^{\dagger}Bb1y1
[M{i,1},B{i,1},ob(1),mb(1),dB(1)] = GenerateMeas(0,theta1);
dA(1) = length(rho{1,1})/dB(1);
theta(i) = theta1;
theta1 = theta1 + pi/64;
i=i+1;
end

r = length(theta);
%Generate assemblage for measurements M{1} on state rho, to get assemblage
%sigma{j,1} for each angle theta(j).


for j=1:r
    
%Generate Measurement set for second round of measurments, with angle
%theta2e
[M{j,2},B{j,2},ob(2),mb(2),dB(2)] = GenerateMeas(0,theta2);

dA(2) = dB(2);
%Generate assemblage (sigma_{b_2|y_2,b_1|y_1=0}) for new measurement set if
%Bob measured in X basis in first round (y_1 = 0)    
 
%Generate (unnormalised) Post measurement states after first round of measurements,
%rho^2_{b1|y1}

for b1 = 1:ob(1)
    for y1 = 1:mb(1)
     
    sigma{j,1}(:,:,b1,y1)= ...
          PartialTrace(kron(eye(dA(1)),M{j,1}(:,:,b1,y1))*(rho{1,1}),2);
            
     rho{j,2}(:,:,b1,y1) = post_meas_new(B{j,1}(:,:,b1,y1),rho{1,1},0);
     
       for b2=  1:ob(2)
            for y2 = 1:mb(2)
          
     sigma{j,2}(:,:,b1,b2,y1,y2)= ...
         PartialTrace(kron(eye(dA(2)),M{j,2}(:,:,b2,y2))*(rho{j,2}(:,:,b1,y1)),2);
            
            rho{j,3}(:,:,b1,b2,y1,y2) =...
                 post_meas_new(B{j,2}(:,:,b2,y2),rho{j,2}(:,:,b1,y1),1); 
 
            end
       end
    end
end

    

end
