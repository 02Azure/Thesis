function y = SSSS_flexuralaniso(Xi)
global U b t_ply Nx maxtskin
%SSSS Boundary, Long Flexural Anisotropic Plate, 0 +-45 90 layer only
%matriks Xi       
XiD = [1 Xi(1) Xi(2) 0 0;...
       1 -Xi(1) Xi(2) 0 0;...
       0 0 -Xi(2) 1 0;...
       0 0 -Xi(2) 0 1;...
       0 0.5*Xi(3) 0 0 0;...
       0 0.5*Xi(3) 0 0 0;];    
   
D_= XiD*U; %for non dimensional calculation, excluding h^3/12
% (D(1,1)) -> D(1)
% (D(2,2)) -> D(2)
% (D(1,2)) -> D(3)
% (D(3,3)) -> D(4)
% (D(1,3)) -> D(5)
% (D(2,3)) -> D(6)

%non-dimensional variable calculation
% alpha = (D_(2)/D_(1))^0.25;
beta = (D_(3)+2*D_(4))/((D_(1)*D_(2))^0.5);
gamma = D_(5)/(((D_(1)^3)*D_(2))^0.25);
delta = D_(6)/(((D_(2)^3)*D_(1))^0.25);   
Kx = ((2*(1-(4*delta*gamma)-(3*(delta^4))+(2*(delta^2)*beta))^0.5)+2*(beta-3*(delta^2)));

for n = t_ply*6:t_ply:maxtskin
      if(((n^3)/12)*Kx*((pi()^2)*((D_(1)*D_(2))^0.5)/(b)^2)>=Nx||n==maxtskin)
         Nxcr = ((n^3)/12)*Kx*((pi()^2)*((D_(1)*D_(2))^0.5)/(b)^2);
         RF = 1-(Nxcr/Nx);
        break
      end
end
y = (n/t_ply)+(RF/10);