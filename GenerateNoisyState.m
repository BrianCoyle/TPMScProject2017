function [rhonoisy] = GenerateNoisyState(epsilon)

%This function generates the maximally entangled phi^+ state with
%probability 1-epsilon, with probability epsilon/3 to generate any of the
%other Bell states as a raw entangled state rhonoisy{1}. The other states
%rhonoisy{2,3,4} are the states obtained using the purification protocol in
%''Minimally complex ion traps as modules for quantum communication and
%computing''

phiplus = (1/2)*[1,0,0,1 ; 0,0,0,0 ; 0,0,0,0 ;1,0,0,1];
phiminus = (1/2)*[1,0,0,-1 ; 0,0,0,0 ; 0,0,0,0 ;-1,0,0,1];
psiplus = (1/2)*[0,0,0,0; 0,1,1,0 ; 0,1,1,0;0,0,0,0];
psiminus = (1/2)*[0,0,0,0; 0,1,-1,0 ; 0,-1,1,0;0,0,0,0];

r1{1} = epsilon/3;
r2{1} = epsilon/3;
r3{1} = epsilon/3;
%Raw entangled state
rhonoisy{1} = (1-r1{1}-r2{1}-r3{1})*phiplus+...
            r1{1}*phiminus +r2{1}*psiplus+r3{1}*psiminus;
        
r1{2} = (2/3)*epsilon+(2/9)*epsilon^2;
r2{2} = (2/9)*epsilon^2;
r3{2} = (2/9)*epsilon^2;
%State after one round of purification protocol
rhonoisy{2} = (1-r1{2}-r2{2}-r3{2})*phiplus+...
            r1{2}*phiminus +r2{2}*psiplus+r3{2}*psiminus;
      
r1{3} = (4/9)*epsilon^2;
r2{3} = (4/9)*epsilon^2;
r3{3} = (8/27)*epsilon^3;

%State after two rounds of purification protocol
rhonoisy{3} = (1-r1{3}-r2{3}-r3{3})*phiplus+...
            r1{3}*phiminus +r2{3}*psiplus+r3{3}*psiminus;
        
r1{4} = (2/9)*epsilon^2;
r2{4} =  (8/27)*epsilon^3;
r3{4} = (8/27)*epsilon^3;

%State after three rounds of purification protocol
rhonoisy{4} = (1-r1{4}-r2{4}-r3{4})*phiplus+...
            r1{4}*phiminus +r2{4}*psiplus+r3{4}*psiminus;
end