clc
clear
global Nx h_web t_ply U XiDopt maxtweb b a f

    Nx = 500; %Critical Buckling Compression load Constraint, N/mm

    %Geometri Pelat
    a = 300; %panjang pelat, mm
    b = 170; %lebar pelat skin, mm
    f = 80; %lebar total flange, mm
    h_web = 50; %tinggi web stiffener, mm
    maxtweb = 7; %tebal web maksimum yang diperbolehkan

    %Input Properti Material
    E11 = 135000; %MPa
    E22 = 10000; %MPa
    G12 = 5000; %MPa
    v12 = 0.3; %Poisson's Ratio
    t_ply = 0.125; %tebal per ply, mm

    %Perhitungan Matriks Q dan invariant material
    Q11 = (E11^2)/(E11-(E22*(v12^2)));
    Q22 = (E11*E22)/(E11-(E22*(v12^2)));
    Q12 = v12*Q22;
    Q66 = G12;
    Q = [ Q11;Q22;Q12;Q66 ];
    inva = [3/8 3/8 1/4 1/2;... 
            1/2 -1/2 0 0;... 
            1/8 1/8 -1/4 -1/2;...
            1/8 1/8 3/4 -1/2;... 
            1/8 1/8 -1/4 1/2 ];
    U = inva*Q;

%GA lamination optimum parameter search - Stiffener+Skin----------------------
ObjectiveFunction = @Lamination_ParaSearch;
X0 = [0 -1 0 0 -1 0]; %gene awal - opsional
A = [2,-1,0,0,0,0; %kiri ke kanan: Xi_D skin 123, Xi_D stiffener 123
     0,1,2,0,0,0;
     -2,-1,0,0,0,0;
     0,1,-2,0,0,0;
     0,0,0,2,-1,0;
     0,0,0,0,1,2;
     0,0,0,-2,-1,0;
     0,0,0,0,1,-2];
B = [1;1;1;1;1;1;1;1];
Aeq = []; beq = [];
LB = [ -1 -1 -1 -1 -1 -1 ];
UB = [ 1 1 1 1 1 1];
nvars = 6;

options = optimoptions('ga','PlotFcn', @gaplotbestf,'InitialPopulation',X0);
[x,fval] = ga(ObjectiveFunction,nvars,A,B,Aeq,beq,LB,UB,[],[],options)                                       
                                      
XiDopt(1,:) = x(1:3); %Xi_skin
XiDopt(2,:) = x(4:6); %Xi_web
total_plyskin = round(fval)/2; 
eigen_val = (round(fval)-fval)*10;

clear x A b Aeq beq LB UB

%GA Most fit Stacking sequence search - skin and web-------------------------
options = optimoptions('ga','PlotFcn', @gaplotbestf);
ObjectiveFunction = @seqsearch;
nvars = ceil(total_plyskin/2)+1;
A = [];
B = [];
Aeq = []; beq = [];

if mod(total_plyskin,2) == 1
   LB(1:nvars) = 1; %mid symteric
   UB(1:nvars-1) = 4;
   UB(nvars) = 1;
else
   LB(1:nvars-1) = 1;
   LB(nvars) = 0; %symetric
   UB(1:nvars-1) = 4;
   UB(nvars) = 0;
end

Intcon = 1:nvars;
error = 1;
n = 1;
max = 15;

while error > 0.00001 && n ~= max                    
    [x,fval,output] = ga(ObjectiveFunction,...
                             nvars,A,B,Aeq,beq,LB,UB,[],Intcon,options);
   
    n = n+1;
    if fval <= error 
       clear seq_skin
       seq_skin = x;
    end
    error = fval;
end

%Convert to degree
 for n = 1:length(seq_skin)-1
    if  seq_skin(n) == 1
        seq_skin(n) = 0;
    end
    if  seq_skin(n) == 2
        seq_skin(n) = 45;
    end
    if  seq_skin(n) == 3
        seq_skin(n) = -45;
    end
    if  seq_skin(n) == 4
        seq_skin(n) = 90;
    end
 end
 clear LB UB
 
 clc
 disp('Optimization Finished')
 disp('most optimum stacking sequence for skin:')
 disp(seq_skin(1:end-1))
 
 if(seq_skin(end)==0)
    disp('Symmetric')
    disp('Total number of layers:')
    disp((length(seq_skin)-1)*2)
 else
    disp('Mid-plane Symmetric')
    disp('Total number of layers:')
    disp(((length(seq_skin)-2)*2)+1)
 end
 eig_val
 
