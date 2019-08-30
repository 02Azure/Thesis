function y = Lamination_ParaSearch(Xi)
global U h_web t_ply Nx maxtweb maxtskin b f

%skin: SSSS Boundary, Long Flexural Anisotropic Plate, 0 +-45 90 layer only
%web: SSSF Boundary, Long Flexural Anisotropic Plate, 0 +-45 90 layer only

%matriks Xi skin dan web
XiD = [1 Xi(1) Xi(2) 0 0;...
       1 -Xi(1) Xi(2) 0 0;...
       0 0 -Xi(2) 1 0;...
       0 0 -Xi(2) 0 1;...
       0 0.5*Xi(3) 0 0 0;...
       0 0.5*Xi(3) 0 0 0];    
       
D = XiD*U; %excluding h^3/12, matriks D web dan skin
% (D(1,1)) -> D(1)
% (D(2,2)) -> D(2)
% (D(1,2)) -> D(3)
% (D(3,3)) -> D(4)
% (D(1,3)) -> D(5)
% (D(2,3)) -> D(6)

%Critical Buckling Load Web, without h
Pcr1w = (12/h_web^2)*(D(4)-((D(6)^2)/D(2)));
Pcr2w = (1/h_web^2) *(12*D(4)-10*((D(5)^2)/(D(1))));

%non-dimensional skin-flange variable for Pcr skin calculation
% alpha = (Ds(2)/Ds(1))^0.25;
beta = (D(3)+2*D(4))/((D(1)*D(2))^0.5);
gamma = D(5)/(((D(1)^3)*D(2))^0.25);
delta = D(6)/(((D(2)^3)*D(1))^0.25);   

%Kx skin
Kx = ((2*(1-(4*delta*gamma)-(3*(delta^4))+(2*(delta^2)*beta))^0.5)+2*(beta-3*(delta^2)));

for n = t_ply*6:t_ply:maxtweb %n = tebal ply web/skin
%Pcr calculation web
   if(Pcr1w <= Pcr2w)
      Pcrw = Pcr1w*((n^3)/12);
   else
      Pcrw = Pcr2w*((n^3)/12);
   end

%Pcr calculation skin + flange
   Pcrs = ((n^3)/12)*Kx*((pi()^2)*((D(1)*D(2))^0.5)/(b)^2);

   if(Pcrs>=Nx && Pcrw>=Nx ||n==maxtweb)
      RFs = 1-(Pcrs/Nx);
      RFw = 1-(Pcrw/Nx);
      break
   end
end

y = (n/t_ply)+((RFs+RFw)/20);
