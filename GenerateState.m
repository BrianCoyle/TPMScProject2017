function [rho] = GenerateState(angle)

rho = [cos(angle)^2,0,0,cos(angle)*sin(angle); 0,0,0,0 ; 0,0,0,0;...
       cos(angle)*sin(angle),0,0,sin(angle)^2];

end