function [sigma,theta,rho] = GenAssemblagesThreeRounds(rho1,phi2)
%This function calculates the assemblages for 3 rounds of measurements, and
%assigns them to a cell array, sigma. The first section calculates the
%assemblage for only 1 measurement angle in the first round
%(theta1 and phi1 needs to be specified),
%and the second calculates it for a range of initial angles. 

%---------------------------------------------------------------------------
%Uncomment the following to use only one measurement angle for the initial
%measurement set
%---------------------------------------------------------------------------

% ------------------------------------------------------------------------
% % Uncomment the following to use a range of angles for the initial
% % measurement set
% % ---------------------------------------------------------------------------

% Genertate assemblage sigma_{b1|y1} for a range of measurement angles  in
% the first round. Assign the initial state to cell array

i=1;
theta1 =0.05*pi;
rho{1,1} = rho1;
while theta1 < pi/4
theta1 = theta1 + pi/64;
[M{i,1},B{i,1},ob(1),mb(1),dB(1)] = GenerateMeas(0,theta1);
dA(1) = length(rho{1,1})/dB(1);
theta(i) = theta1;
i=i+1;
end

r = length(theta);
% Generate assemblage for measurements M{1} on state rho, to get assemblage
% sigma{j,1} for each angle theta(j).

for j=1:r
    
% Generate Measurement set for second round of measurments, with angle
% theta2
[M{j,2},B{j,2},ob(2),mb(2),dB(2)] = GenerateMeas(phi2,0);

dA(2) = dB(2);

 [M{j,3},B{j,3},ob(3),mb(3),dB(3)] = GenerateMeas(0,0);
 dA(3) = dB(3);


for b1 = 1:ob(1)
    for y1 = 1:mb(1)
     
 sigma{j,1}(:,:,b1,y1)= ...
     PartialTrace(kron(eye(dA(1)),M{j,1}(:,:,b1,y1))*(rho{1,1}),2);
 rho{j,2}(:,:,b1,y1) = GeneratePostMeasState(B{j,1}(:,:,b1,y1),rho{1,1},0);     

         for b2=  1:ob(2)
            for y2 = 1:mb(2)
          
        sigma{j,2}(:,:,b1,b2,y1,y2)= ...
         PartialTrace(kron(eye(dA(2)),M{j,2}(:,:,b2,y2))*(rho{j,2}(:,:,b1,y1)),2);
            
         rho{j,3}(:,:,b1,b2,y1,y2) =...
           GeneratePostMeasState(B{j,2}(:,:,b2,y2),rho{j,2}(:,:,b1,y1),0);
       
              
                  for b3 = 1:ob(3)
                     for y3 = 1:mb(3)
          
                sigma{j,3}(:,:,b1,b2,b3,y1,y2,y3)= ...
                     PartialTrace(kron(eye(dA(3)),M{j,3}(:,:,b3,y3))*...
                                        (rho{j,3}(:,:,b1,b2,y1,y2)),2);
            
                     rho{j,4}(:,:,b1,b2,b3,y1,y2,y3) =...
                     GeneratePostMeasState(B{j,2}(:,:,b3,y3),rho{j,3}(:,:,b1,b2,y1,y2),1);  
  
                     end
                  end 
            end
         end
    end
end
    
    
end
 
 
    










