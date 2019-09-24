clc
clear
syms x y a b K L J D11 D12 D66 D22 Nx f D11sf D12sf D66sf D22sf
%Hanya untuk tebal bervariasi secara konstan pada bagian pelat tertentu
%f panjang setengah flange satu stiffened panel skin
terms = 9; %number of terms
m = terms;
n = m;

% matriks Xi skin, Xi_4 A dan D = 0
Xi = [ 0 -1 -0.027778 0 -1 0 ];
XiD = [1 Xi(1) Xi(2) 0 0;...
       1 -Xi(1) Xi(2) 0 0;...
       0 0 -Xi(2) 1 0;...
       0 0 -Xi(2) 0 1;...
       0 0.5*Xi(3) 0 0 0;...
       0 0.5*Xi(3) 0 0 0];    
   
XiA = [1 Xi(4) Xi(5) 0 0;...
       1 -Xi(4) Xi(5) 0 0;...
       0 0 -Xi(5) 1 0;...
       0 0 -Xi(5) 0 1;...
       0 0.5*Xi(6) 0 0 0;...
       0 0.5*Xi(6) 0 0 0];
   
E11 = 135000; %MPa
E22 = 10000; %MPa
G12 = 5000; %MPa
% test isotropic mat
% E11 = 72000; %MPa
% E22 = 72000; %MPa
% G12 = 26900; %MPa
v12 = 0.3; %Poisson'n Ratio
t_ply = 0.125; %tebal per ply, mm
h = 12*t_ply;
% Perhitungan Matriks Q dan invariant material
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
D = XiD*U*(h^3)/12;
A = XiA*U*h;

% stiffness flange
Xif = [0 -1 -0.444444444  0 -1 0];

XiDf = [1 Xif(1) Xif(2) 0 0;...
       1 -Xif(1) Xif(2) 0 0;...
       0 0 -Xif(2) 1 0;...
       0 0 -Xif(2) 0 1;...
       0 0.5*Xif(3) 0 0 0;...
       0 0.5*Xif(3) 0 0 0];    
   
XiAf = [1 Xif(4) Xif(5) 0 0;...
       1 -Xif(4) Xif(5) 0 0;...
       0 0 -Xif(5) 1 0;...
       0 0 -Xif(5) 0 1;...
       0 0.5*Xif(6) 0 0 0;...
       0 0.5*Xif(6) 0 0 0];
hf = 6*t_ply;
Df = XiDf*U*(hf^3)/12;
Af = XiAf*U*hf;
% Dummy = Dsf;
% Dsf = D;
% D = Dummy;
% (D(1,1)) -> D(1)
% (D(2,2)) -> D(2)
% (D(1,2)) -> D(3)
% (D(3,3)) -> D(4)
% (D(1,3)) -> D(5)
% (D(2,3)) -> D(6)
%mencari kekakuan tumpukan skin+flange
Asf = [ A(1) A(3) A(5); A(3) A(2) A(6); A(5) A(6) A(4)] + [Af(1) Af(3) Af(5); Af(3) Af(2) Af(6); Af(5) Af(6) Af(4)];
Bsf = hf*0.5*[A(1) A(3) A(5); A(3) A(2) A(6); A(5) A(6) A(4)]-h*0.5*[Af(1) Af(3) Af(5); Af(3) Af(2) Af(6); Af(5) Af(6) Af(4)];
Dsf = [ D(1) D(3) D(5); D(3) D(2) D(6); D(5) D(6) D(4)] + [ Df(1) Df(3) Df(5); Df(3) Df(2) Df(6); Df(5) Df(6) Df(4)]...
      + (hf^2)*0.25*[ A(1) A(3) A(5); A(3) A(2) A(6); A(5) A(6) A(4)] + (h^2)*0.25*[Af(1) Af(3) Af(5); Af(3) Af(2) Af(6); Af(5) Af(6) Af(4)];
Dsmea = Dsf-(Bsf*inv(Asf)*Bsf);
clear Dsf
Dsf(1) = Dsmea(1,1);
Dsf(2) = Dsmea(2,2);
Dsf(3) = Dsmea(1,2);
Dsf(4) = Dsmea(3,3);
Dsf(5) = Dsmea(1,3);
Dsf(6) = Dsmea(2,3);
Dsf = Dsf';
% (D(1,1)) -> D(1)
% (D(2,2)) -> D(2)
% (D(1,2)) -> D(3)
% (D(3,3)) -> D(4)
% (D(1,3)) -> D(5)
% (D(2,3)) -> D(6)

panjang = 300;
lebar = 170;
flange = 45; %lebar setengah flange, max setengah lebar skin
load = 100;


a_K = K*pi/a;
b_L = L*pi/b;
b_J = J*pi/b;
fun_sf1 = (D11sf*(a_K^4)+D22sf*(b_J^2)*(b_L^2)+D12sf*((a_K^2)*(b_L^2)+(a_K^2)*(b_J^2)))*((1/(b_J-b_L))*sin(b_J*f-b_L*f)-(1/(b_J+b_L))*sin(b_J*f+b_L*f));
fun_sf2 = 4*D66sf*(a_K^2)*(b_L*b_J)*((1/(b_J-b_L))*sin(b_J*f-b_L*f)+(1/(b_J+b_L))*sin(b_J*f+b_L*f));
fun_sf3 = (D11sf*(a_K^4)+D22sf*(b_J^2)*(b_L^2)+D12sf*((a_K^2)*(b_L^2)+(a_K^2)*(b_J^2)))*((1/(b_J-b_L))*(sin((b_J-b_L)*b)-sin((b_J-b_L)*(b-f)))-((1/(b_J+b_L))*(-sin((b_J+b_L)*(b-f)))));
fun_sf4 = 4*D66sf*(a_K^2)*(b_L*b_J)*((1/(b_J-b_L))*(sin((b_J-b_L)*b)-sin((b_J-b_L)*(b-f)))+((1/(b_J+b_L))*(-sin((b_J+b_L)*(b-f)))));
fun_sfl = (a/8)*(fun_sf1+fun_sf2);
fun_sfr = (a/8)*(fun_sf3+fun_sf4);
fun_s1 = (D11*(a_K^4)+D22*(b_J^2)*(b_L^2)+D12*((a_K^2)*(b_L^2)+(a_K^2)*(b_J^2)))*((1/(b_J-b_L))*(sin((b_J-b_L)*(b-f))-sin(b_J*f-b_L*f))-(1/(b_J+b_L))*(sin((b_J+b_L)*(b-f))-sin(b_J*f+b_L*f)));
fun_s2 = 4*D66*(a_K^2)*(b_L*b_J)*((1/(b_J-b_L))*(sin((b_J-b_L)*(b-f))-sin(b_J*f-b_L*f))+(1/(b_J+b_L))*(sin((b_J+b_L)*(b-f))-sin(b_J*f+b_L*f)));
fun_s = (a/8)*(fun_s1+fun_s2);
funLHS = fun_sfl+fun_sfr+fun_s;
funRHS = -0.5*Nx*(a_K^2)*(a*b/4);

funLHS = subs(funLHS,b,lebar);
funLHS = subs(funLHS,f,flange);
funLHS = subs(funLHS,a,panjang);
funLHS = subs(funLHS,D11,D(1));
funLHS = subs(funLHS,D12,D(3));
funLHS = subs(funLHS,D66,D(4));
funLHS = subs(funLHS,D22,D(2));
funLHS = subs(funLHS,D11sf,Dsf(1));
funLHS = subs(funLHS,D12sf,Dsf(3));
funLHS = subs(funLHS,D66sf,Dsf(4));
funLHS = subs(funLHS,D22sf,Dsf(2));
funRHS = subs(funRHS,b,lebar);
funRHS = subs(funRHS,Nx,load);
funRHS = subs(funRHS,a,panjang);
    baris = 1;
    kolom = 1;
LHS(m*n,m*n) = 0;
RHS(m*n,m*n) = 0; 

tic
for K_ = 1:m
    funnew = subs(funLHS,K,K_);
    for L_ = 1:n
        funnew2 = subs(funnew,L,L_);
        for i_=1:m
            for j_=1:n
               if((i_==K_)&&(j_~=L_))
                  LHS(baris,kolom) = subs(funnew2,J,j_);
               elseif((i_==K_)&&(j_==L_))
                  LHS(baris,kolom) = limit(funnew2,J,j_);
                  RHS(baris,kolom) = subs(funRHS,K,K_);
               else
                  LHS(baris,kolom) = 0;
               end
               kolom = kolom+1;
            end
        end
    kolom = 1;
    baris = baris+1;
    end
end
toc

tic
LHS = double(LHS);
RHS = double(RHS);
lambda = eig(LHS,RHS);
toc2 = toc
max(lambda)
sort(lambda(:),'descend')
