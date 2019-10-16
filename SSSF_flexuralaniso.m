function y = SSSF_flexuralaniso(Xi)
global U h_web t_ply Nx maxtweb
%SSSF Boundary, Long Flexural Anisotropic Plate, 0 +-45 90 layer only
%matriks Xi
XiD = [1 Xi(1) Xi(2) 0 0;...
       1 -Xi(1) Xi(2) 0 0;...
       0 0 -Xi(2) 1 0;...
       0 0 -Xi(2) 0 1;...
       0 0.5*Xi(3) 0 0 0;...
       0 0.5*Xi(3) 0 0 0;];    
   
D = XiD*U; %excluding h^3/12
% (D(1,1)) -> D(1)
% (D(2,2)) -> D(2)
% (D(1,2)) -> D(3)
% (D(3,3)) -> D(4)
% (D(1,3)) -> D(5)
% (D(2,3)) -> D(6)

Pcr1 = (12/h_web^2)*(D(4)-((D(6)^2)/D(2)));
Pcr2 = (1/h_web^2) *(12*D(4)-10*((D(5)^2)/(D(1))));

for n = t_ply*6:2*t_ply:maxtweb
      if(Pcr1*((n^3)/12)>=Nx && Pcr2*((n^3)/12)>=Nx ||n==maxtweb)
         Nxcr(1) = Pcr1*((n^3)/12);
         Nxcr(2) = Pcr2*((n^3)/12);
         RF = 1-(min(Nxcr)/Nx);
        break
      end
end
y = (n/t_ply)+(RF/10);
