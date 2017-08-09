function  DeterministicTwoRound = GenerateDeterministicTwoRounds(ob,mb,dA)

%This function calculates the array of possible deterministic strategies
%that Eve can produce. It is based from genSinglePartyArray.m in
%''Steering: a review with focus on semidefinite programming'', but adapted for two measurement rounds.
%This works for two round of measurements.

%Total number of deterministic strategies is Ndetone x Ndettwo for 2 rounds
Ndetone = ob(1)^mb(1);
Ndettwo = ob(2)^mb(2);

DeterministicTwoRound  = zeros(dA(2),dA(2),ob(1),ob(2),mb(1),mb(2),Ndetone,Ndettwo); 

for lamone = 0:Ndetone-1
    for lamtwo = 0:Ndettwo-1
    %generates the string of outcomes b1 (for each y1), for the given
    %variable lamone
    lamdecone = dec2base(lamone,ob(1),mb(1))-'0'; 
   
    %generates the string of outcomes b2 (for each y2), for the given
    %variable lamtwo
    lamdectwo = dec2base(lamtwo,ob(2),mb(2))-'0';
   

        for i = 1:dA(2)
           for j = 1:dA(2)
              for y1 = 1:mb(1)
                 for y2 = 1:mb(2)
                     for b1 = 0:ob(1)-1
                         for b2 = 0:ob(2)-1
              DeterministicTwoRound(i,j,b1+1,b2+1,y1,y2,lamone+1,lamtwo+1) ...
                              = (lamdecone(y1) == b1 && lamdectwo(y2) == b2);
                         % probability = 1 if b1 = lamdec(y1) AND b2 = lamdec(y2), 0 otherwise
                         end
                     end
                 end
              end
           end
         end
        
    end
end
    