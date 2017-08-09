function [F,v] = GenerateFunctionalThreeRounds(sigma)


%This function generates the Steering inequalities that are maximally violated
%by the assemblage, sigma, after two rounds of measurements
%GenerateFunctionalOneRoundSW generates the steering inequality for the first round
%It also calculates the observed violation of that steering inequality, v(1),
%by the assemablage sigma{1} generated for the first round.
%Similarly, for the second round GenerateSteeringFunctionalTwoRounds
%generates the steering funcional maximally violated by the assemblage
%sigma{2}, and its observed violation, v(3)
%For the third round, the steering inequality defined in 
%''Quantifying Einstein-Podolsky-Rosen steering'' is calculated and again
%the observed violation, v(3), is calculated, which is acheived by the
%assemblage generated after the third round sigma{3}

[dA(1),dA(1),ob(1),mb(1)] = size(sigma{1});
[dA(2),dA(2),~,ob(2), ~ ,mb(2)] = size(sigma{2});
[dA(3),dA(3),~,~,ob(3),~, ~ ,mb(3)] = size(sigma{3});

%Generates Steering functional after Round one
[Fprime{1}] = GenerateSteeringFunctionalOneRound(sigma{1});

F{1} = real(Fprime{1});
tr1 = zeros(ob(1),mb(1));

%Generates Steering functional after Round two
[~, Fprime{2}] = GenerateSteeringFunctionalTwoRounds(sigma);

 F{2} = real(Fprime{2});
 tr2 = zeros(ob(1),ob(2),mb(1),mb(2));

 
 %Generates Steering functional after Round three
  F{3} = zeros(size(sigma{3}));
  tr3 = zeros(ob(1),ob(2),ob(3),mb(1),mb(2),mb(3));
 
 alpha = 100;
 
for b1=1:ob(1)
    for y1=1:mb(1)
        
        tr1(b1,y1) = trace(sigma{1}(:,:,b1,y1)*F{1}(:,:,b1,y1));
       
        for b2 = 1:ob(2)
            for y2 = 1:mb(2)
                
        tr2(b1,b2,y1,y2) =...
            trace(sigma{2}(:,:,b1,b2,y1,y2)*F{2}(:,:,b1,b2,y1,y2));
        
        
                for b3 = 1:ob(3)
                  for y3=1:mb(3)
                       %If the assemblage element has 0 probability of being
                       %produced, given the initial state state
                     if trace(sigma{3}(:,:,b1,b2,b3,y1,y2,y3)) ==0
                             F{3}(:,:,b1,b2,b3,y1,y2,y3) = alpha*(eye(dA(3)) -[0,0;0,0]); 
                     else
                          %If the assemblage element has nonzero probability of being
                          %produced, given the initial state state
                          F{3}(:,:,b1,b2,b3,y1,y2,y3) = alpha*(eye(dA(3)) - ...
                            (sigma{3}(:,:,b1,b2,b3,y1,y2,y3)/trace(sigma{3}(:,:,b1,b2,b3,y1,y2,y3))));
                     end
                        tr3(b1,b2,b3,y1,y2,y3) = ...
                          trace(sigma{3}(:,:,b1,b2,b3,y1,y2,y3)*F{3}(:,:,b1,b2,b3,y1,y2,y3));
                  end
                end
             end
        end
    end
end
  %Calculate observed violations, v(1), v(2), v(3)
 v(1)  = sum(sum(tr1,1),2);

 v(2)  = sum(sum(sum(sum(tr2,1),2),3),4);

 v(3) = sum(sum(sum(sum(sum(sum(tr3,1),2),3),4),5),6);

 
 
 
 






