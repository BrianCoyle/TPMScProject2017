function [F,v] = GenerateFunctionalTwoRounds(sigma)

%This function generates the Steering inequalities that are maximally violated
%by the assemblage, sigma, after two rounds of measurements
%GenerateFunctionalOneRoundSW generates the steering inequality for the first round
%It also calculates the observed violation of that steering inequality, v(1),
%by the assemablage sigma{1} generated for the first round. 
%For the second round, the steering inequality defined in 
%''Quantifying Einstein-Podolsky-Rosen steering'' is calculated and again
%the observed violation, v(2), is calculated, which is acheived by the
%assemblage generated after the second round sigma{2}

[dA(1),dA(1),ob(1),mb(1)] = size(sigma{1});
[dA(2),dA(2),~,ob(2),~,mb(2)] = size(sigma{2});

%Generates Steering functional after Round one
[Fprime{1}] = GenerateSteeringFunctionalOneRound(sigma{1});

F{1} = real(Fprime{1});

%Calculates observed steering inequality violation after one round.
tr1 = zeros(ob(1),mb(1));

for b1=1:ob(1)
    for y1=1:mb(1)
        
        tr1(b1,y1) = trace(sigma{1}(:,:,b1,y1)*F{1}(:,:,b1,y1));
       
    end
end

v(1)  = sum(sum(tr1,1),2);

%Calculates steering inequality and violation after Round 2.
F{2} = zeros(size(sigma{2}));
alpha = 100;
tr2 = zeros(ob(1),ob(2),mb(1),mb(2));

 for y1 = 1:mb(1)
     for y2 = 1:mb(2)
         for b1 = 1:ob(1)
             for b2 = 1:ob(2)
                 %If the assemblage element has 0 probability of being
                 %produced, given the initial state state
                 if trace(sigma{2}(:,:,b1,b2,y1,y2)) ==0
                   F{2}(:,:,b1,b2,y1,y2) = alpha*(eye(dA(2)) -[0,0;0,0]); 
                 else
                  %Otherwise, if the assemblage element has nonzero
                  %probability of being produced, given the initial state state
    F{2}(:,:,b1,b2,y1,y2) = alpha*(eye(dA(2)) - ...
                (sigma{2}(:,:,b1,b2,y1,y2)/trace(sigma{2}(:,:,b1,b2,y1,y2))));
                 end
    tr2(b1,b2,y1,y2) = trace(sigma{2}(:,:,b1,b2,y1,y2)*F{2}(:,:,b1,b2,y1,y2));
    
             end
         end
     end
 end
 v(2)  = sum(sum(sum(sum(tr2,1),2),3),4);

 
 
 






