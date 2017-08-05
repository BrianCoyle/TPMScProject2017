function [Mby,Bby,ob,mb,dB] = GenerateMeas(phi,theta)

%This function makes a set of measurements on Bob_i's side 
%for an input angle theta in the noisy X basis and an input angle 
%phi for a measurement in the noisy Z basis. theta, phi =0 
%correspond to the usual projective X and Z measurements

%Define non-projective Z Kraus operator measurements, b = 1 implies Bob
%measures in noisy Z basis with angle phi
c1 = cos(phi);
s1 = sin(phi);
Bby(:,:,1,1) = [c1,0;0,s1];
Bby(:,:,2,1) = [s1,0;0,c1];

%Define non-projective X Kraus operator measurements, b = 2 implies Bob
%measures in noisy X basis with angle theta
c2 = cos(theta);
s2 = sin(theta);

Bby(:,:,1,2) = (1/2)*[c2+s2,c2-s2; c2-s2,c2+s2];
Bby(:,:,2,2) = (1/2)*[c2+s2,s2-c2; s2-c2,c2+s2];

%Define POVM measurement operators, M_{b1|y1} =
%(B_{b1|y1})^\dagger(B_{b1|y1})
Mby(:,:,1,1) = (Bby(:,:,1,1)')*(Bby(:,:,1,1));
Mby(:,:,2,1) = (Bby(:,:,2,1)')*(Bby(:,:,2,1));
Mby(:,:,1,2) = (Bby(:,:,1,2)')*(Bby(:,:,1,2));
Mby(:,:,2,2) = (Bby(:,:,2,2)')*(Bby(:,:,2,2));

%ob is the number of measurement outcomes, mb is the number of possible
%measurement choices Bob has
[dB,dB,ob,mb] = size(Mby);
