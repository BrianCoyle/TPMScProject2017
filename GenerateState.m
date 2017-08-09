function [rho] = GenerateState(zeta)

%This function generates the state cos(zeta)|00> + sin(zeta)|11> by
%taking an input angle, zeta, in the form of a density matrix rho
rho = [cos(zeta)^2,0,0,cos(zeta)*sin(zeta);...
          0,0,0,0;...
          0,0,0,0;...
       cos(zeta)*sin(zeta),0,0,sin(zeta)^2];

end