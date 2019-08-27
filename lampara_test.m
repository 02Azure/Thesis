%Menghitung Lamination Parameters
clc
clear
transform = @(theta,C)([cosd(theta)^2 sind(theta)^2 2*sind(theta)*cosd(theta);... %T1
                        sind(theta)^2 cosd(theta)^2 -2*sind(theta)*cosd(theta);...
                       -sind(theta)*cosd(theta) sind(theta)*cosd(theta) cosd(theta)^2-sind(theta)^2]\C*...
                       [cosd(theta)^2 sind(theta)^2 sind(theta)*cosd(theta);... %T2
                        sind(theta)^2 cosd(theta)^2 -sind(theta)*cosd(theta);...
                       -2*sind(theta)*cosd(theta) 2*sind(theta)*cosd(theta) cosd(theta)^2-sind(theta)^2]);

transform1 = @(theta)([cosd(theta)^2 sind(theta)^2 2*sind(theta)*cosd(theta);... %T1
                            sind(theta)^2 cosd(theta)^2 -2*sind(theta)*cosd(theta);...
                            -sind(theta)*cosd(theta) sind(theta)*cosd(theta) cosd(theta)^2-sind(theta)^2]);
                        
transform2 = @(theta)([cosd(theta)^2 sind(theta)^2 sind(theta)*cosd(theta);... %T2
                           sind(theta)^2 cosd(theta)^2 -sind(theta)*cosd(theta);...
                           -2*sind(theta)*cosd(theta) 2*sind(theta)*cosd(theta) cosd(theta)^2-sind(theta)^2]);
                       
mat_name = "example_material"; %nama material
%dari ref
% E11 = 127.5538; %E1(GPa)
% E22 = 11.30747; %E2(GPa)
% v12 = 0.3; %Poisson ratio v12
% G12 = 5.99848; %G12(GPa)

E11 = 135000; %GPa
E22 = 10000; %GPa
G12 = 5000; %GPa
v12 = 0.3; %Poisson's Ratio

v21 = E22*v12/E11;
S = [ 1/E11 -v21/E22 0; -v12/E11 1/E22 0; 0 0 1/G12 ];
C = inv(S);

Q11 = (E11^2)/(E11-(E22*(v12^2)));
Q22 = (E11*E22)/(E11-(E22*(v12^2)));
Q12 = v12*Q22;
Q66 = G12;

inva = [3/8 3/8 1/4 1/2;... 
        1/2 -1/2 0 0;... 
        1/8 1/8 -1/4 -1/2;...
        1/8 1/8 3/4 -1/2;... 
        1/8 1/8 -1/4 1/2 ];
    
Q = [ Q11;Q22;Q12;Q66 ];
U = inva*Q;

%stackseq = 0 45 0 90 -45 90 0 45 0
layernum = 20; %total ply
t(1:layernum) = 0.125; %tebal ply ( mm )
theta = [  45    45    45    45    45    45   -45    45    45    45 45 45 45 -45 45 45 45 45 45 45 ]
% theta = [ 45 -45 45 45 45 45 -45 45 45 -45 45 45 45 45 -45 45 ]
% theta = [ 0 90 90 90 45 45 45 -45 -45 -45 -45 -45 -45 45 45 45 90 90 90 0 ]
% theta = [-45 45 -45 -45 45 45 -45 -45 45 -45];
% theta = [ 0 0 90 90 90 90 90 90 0 0 ]
% theta = [ 0 0 0 90 45 45 45 45 45 45 45 45 -45 -45 -45 -45 -45 -45 -45 -45 -45 -45 -45 -45 -45 -45 -45 -45 45 45 45 45 45 45 45 45 90 0 0 0 ]
%  theta = [ -45 -45 45 -45 -45 -45 -45 45 -45 -45 -45 -45 45 -45 45 -45 -45 45 45 45 45 45 45 -45 -45 45 -45 45 -45 -45 -45 -45 45 -45 -45 -45 -45 45 -45 -45 ]
%   theta =  [0 0 0 90 90 90 45 45 45 45 45 45 -45 -45 -45 -45 -45 -45 -45 -45 -45 -45 -45 -45 45 45 45 45 45 45 90 90 90 0 0 0  ] 
% theta = [45 -45 -45 45 45 -45 45 -45 45 45 -45 -45 -45 45 -45 -45 -45 -45 -45 -45 -45 -45 45 -45 -45 -45 45 45 -45 45 -45 45 45 -45 -45 45 ]
% theta = [ -45 -45 45 45 -45 -45  45 45 -45 45 -45 45 45 -45 45 45 -45 45 45 -45 45 -45 45 45 -45 -45 45 45 -45 -45 ]
%    0.5324    0.7500    0.0312
%4 2 3
% (45/-452/452/-45/45/-45/452/
% -453/45/-454)

% -0.2487   -0.5026    0.1127
% 0     1     5
%pendefinisian nilai z
clear z
z(1:layernum+1) = 0;
z(1) = -0.5*sum(t);

for n = 2:layernum+1
    z(n) = z(n-1)+t(n-1);
end

%Perhitungan matriks ABD
%clear A B D
A(3,3) = 0;
B(3,3) = 0;
D(3,3) = 0;

for n = 1:layernum
    A = A+transform(theta(n),C)*(z(n+1)-z(n));
    B = B+0.5*transform(theta(n),C)*((z(n+1))^2-(z(n)^2));
    D = D+(1/3)*transform(theta(n),C)*((z(n+1))^3-(z(n))^3);
end


%perhitungan lamination parameters
norm_z = 2*z'/sum(t);
Xi_A(1:4) = 0;
Xi_D(1:4) = 0;
for n = 1:layernum
   Xi_A(1) = Xi_A(1)+ norm_z(n+1)*0.5*cosd(2*theta(n))-norm_z(n)*0.5*cosd(2*theta(n));
   Xi_A(2) = Xi_A(2)+ norm_z(n+1)*0.5*cosd(4*theta(n))-norm_z(n)*0.5*cosd(4*theta(n));
   Xi_A(3) = Xi_A(3)+ norm_z(n+1)*0.5*sind(2*theta(n))-norm_z(n)*0.5*sind(2*theta(n));
   Xi_A(4) = Xi_A(4)+ norm_z(n+1)*0.5*sind(4*theta(n))-norm_z(n)*0.5*sind(4*theta(n));
   
   Xi_D(1) = Xi_D(1)+ norm_z(n+1)^3*0.5*cosd(2*theta(n))-norm_z(n)^3*0.5*cosd(2*theta(n));
   Xi_D(2) = Xi_D(2)+ norm_z(n+1)^3*0.5*cosd(4*theta(n))-norm_z(n)^3*0.5*cosd(4*theta(n));
   Xi_D(3) = Xi_D(3)+ norm_z(n+1)^3*0.5*sind(2*theta(n))-norm_z(n)^3*0.5*sind(2*theta(n));
   Xi_D(4) = Xi_D(4)+ norm_z(n+1)^3*0.5*sind(4*theta(n))-norm_z(n)^3*0.5*sind(4*theta(n));
end

inv_pairA = [1 Xi_A(1) Xi_A(2) 0 0;...
             1 -Xi_A(1) Xi_A(2) 0 0;...
             0 0 -Xi_A(2) 1 0;...
             0 0 -Xi_A(2) 0 1;...
             0 0.5*Xi_A(3) Xi_A(4) 0 0;...
             0 0.5*Xi_A(3) -Xi_A(4) 0 0;];
         
inv_pairD = [1 Xi_D(1) Xi_D(2) 0 0;...
             1 -Xi_D(1) Xi_D(2) 0 0;...
             0 0 -Xi_D(2) 1 0;...
             0 0 -Xi_D(2) 0 1;...
             0 0.5*Xi_D(3) Xi_D(4) 0 0;...
             0 0.5*Xi_D(3) -Xi_D(4) 0 0;];         
      
Xi_D
% sum(t)*inv_pairA*U
% disp(A(1,1))
% disp(A(2,2))
% disp(A(1,2))
% disp(A(3,3))
% disp(A(1,3))
% disp(A(2,3))
% 
% ((sum(t)^3)/12)*inv_pairD*U
% disp(D(1,1))
% disp(D(2,2))
% disp(D(1,2))
% disp(D(3,3))
% disp(D(1,3))
% disp(D(2,3))