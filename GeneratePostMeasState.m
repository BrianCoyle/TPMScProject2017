function [rhonew] = GeneratePostMeasState(Bby,rho,arg)

%This function calculates the unnormalised (arg = 0) or normalised (arg  = 1)
%post measurement states, for a given state rho and a set of measurments Mby,
%to be fed into the next round of measurements.

[n,m] = size(Bby);
iden = eye(n,m);


%Computes tensor product M = I \otimes Mby to calculate general measurement
%for both parties, M, from identity on Alice's side and a measurement
%B_{b|y} on Bob's side

genmeas = kron(iden,Bby);


%Calculates post measurement state on initial state \rho to give final
%state \rho_{b|y}, normalised if arg == 1, unnormalised if arg == 0;
if arg == 0
rhonew =  genmeas*rho*genmeas';
                             
else
    if arg == 1
        if trace(genmeas*rho*genmeas') == 0 
rhonew =  [0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
        else
rhonew =  genmeas*rho*genmeas'...
                  /trace(genmeas*rho*genmeas');
        end
    
    else
        disp('Please enter 0 or 1 for arg, post_meas(~,~,arg)')
    end
end
end

