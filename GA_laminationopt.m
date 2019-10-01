clc
clear
global Nx h_web t_ply U XiDopt XiAopt maxtweb b a f flag total_plyweb

    Nx = 500; %Critical Buckling Compression load Constraint, N/mm

    %Geometri Pelat
    a = 300; %panjang pelat, mm
    b = 170; %lebar pelat skin, mm
    f = 40; %lebar setengah  total flange, mm
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
    
%GA lamination optimum parameter search - Web------------------------------
ObjectiveFunction = @SSSF_flexuralaniso;
X0 = [0 -1 0];
A = [2,-1,0;0,1,2;-2,-1,0;0,1,-2];
B = [1;1;1;1];
Aeq = []; beq = [];
LB = [ -1 -1 -1 ];
UB = [ 1 1 1 ];
nvars = 3;
options = optimoptions('ga','PlotFcn', @gaplotbestf,'InitialPopulation',X0);
[x,fval] = ga(ObjectiveFunction,nvars,A,B,Aeq,beq,LB,UB,[],[],options) 
total_plyweb = round(fval);
eigval_web = ((round(fval)-fval)*10)+1;
XiDopt(1,:) = x; %XiD optimal untuk web



%GA lamination optimum parameter search - Skin+flange----------------------
ObjectiveFunction = @SSSS_skinflange;
X0 = [0 -1 0 0 -1 0 0 -1 0 0 -1 0]; %gene awal - opsional
%kiri ke kanan: Xi_D skin 123, Xi_D flange 123, Xi_A skin 123, Xi_A flange 123
A = [2,-1,0,0,0,0,0,0,0,0,0,0; 
     0,1,2,0,0,0,0,0,0,0,0,0;
     -2,-1,0,0,0,0,0,0,0,0,0,0;
     0,1,-2,0,0,0,0,0,0,0,0,0;
     0,0,0,2,-1,0,0,0,0,0,0,0;
     0,0,0,0,1,2,0,0,0,0,0,0;
     0,0,0,-2,-1,0,0,0,0,0,0,0;
     0,0,0,0,1,-2,0,0,0,0,0,0;
     0,0,0,0,0,0,2,-1,0,0,0,0;
     0,0,0,0,0,0,0,1,2,0,0,0;
     0,0,0,0,0,0,-2,-1,0,0,0,0;
     0,0,0,0,0,0, 0,1,-2,0,0,0;
     0,0,0,0,0,0,0,0,0,2,-1,0;
     0,0,0,0,0,0,0,0,0,0,1,2;
     0,0,0,0,0,0,0,0,0,-2,-1,0;
     0,0,0,0,0,0,0,0,0,0,1,-2];
B = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
Aeq = []; beq = [];
LB = [ -1 -1 0 -1 -1 0 -1 -1 -1 -1 -1 -1 ];
UB = [ 1 1 0 1 1 0 1 1 1 1 1 1];
nvars = 12;

options = optimoptions('ga','PlotFcn', @gaplotbestf,'InitialPopulation',X0);
[x,fval] = ga(ObjectiveFunction,nvars,A,B,Aeq,beq,LB,UB,[],[],options)                                       
%                                       
XiDopt(2,:) = x(1:3); %XiD_optimalskin
XiDopt(3,:) = x(4:6); %XiD_optimalflange
XiAopt(1,:) = x(7:9); %XiA_optimalskin
XiAopt(2,:) = x(10:12); %XiA_optimalflange


total_plyskin = round(fval); 
eigval_skin = ((round(fval)-fval)*10)+1;
% 
clear x A b Aeq beq LB UB

% GA Most fit Stacking sequence search - skin-------------------------
options = optimoptions('ga','PlotFcn', @gaplotbestf);
ObjectiveFunction = @seqsearch;
nvars = ceil(total_plyskin/2)+1;
A = [];
B = [];
Aeq = []; beq = [];
flag = 1;

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
 clear LB UB Intcon


%GA Most fit Stacking sequence search - web-------------------------
nvars = ceil(total_plyweb/4)+1;
flag = 2;

if mod(total_plyweb/2,2) == 1
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

while error > 0.00001 && n ~= max                    
    [x,fval,output] = ga(ObjectiveFunction,...
                             nvars,A,B,Aeq,beq,LB,UB,[],Intcon,options);
   
    n = n+1;
    if fval <= error 
       clear seq_web
       seq_web = x;
    end
    error = fval;
end

%Convert to degree
for n = 1:length(seq_web)-1
    if  seq_web(n) == 1
        seq_web(n) = 0;
    end
    if  seq_web(n) == 2
        seq_web(n) = 45;
    end
    if  seq_web(n) == 3
        seq_web(n) = -45;
    end
    if  seq_web(n) == 4
        seq_web(n) = 90;
    end
end
seq_web(n+1:n*2) = seq_web(1:n);
 clc
 disp('Optimization Finished')
 disp('most optimum stacking sequence for skin:')
 disp(seq_skin(1:end-1))
 
 if(seq_skin(end)==0)
    disp('Symmetric')
 else
    disp('Mid-plane Symmetric')
 end
 disp('Total number of layers:')
 disp(total_plyskin)
 eigval_skin
 
disp('most optimum stacking sequence for web:')
 disp(seq_web)
 
 if(mod(total_plyweb,2)==0)
    disp('Symmetric')
 else
    disp('Mid-plane Symmetric')
 end
 disp('Total number of layers:')
 disp(total_plyweb)
 eigval_web
