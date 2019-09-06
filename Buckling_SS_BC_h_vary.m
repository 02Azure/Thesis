clc
clear
syms x y a b m n PI D11 D12 D66 D22 Nx f D11sf D12sf D66sf D22sf
%Hanya untuk tebal dan D konstant sepanjang pelat
terms = 9; %number of terms
mmax = terms;
nmax = mmax;

%skin
funD11 = D11*4*(m*pi/a)^3*(b/(n*pi))*cos(pi*f/b);
funD1266 = (2*D12+4*D66)*4*(m*pi/a)*(n*pi/b)*cos(pi*f/b);
funD22 = D22*4*(n*pi/b)^3*(a/(m*pi))*cos(pi*f/b);
funRHS = Nx*4*b*m/(a*n);

%skin-flange
funD11sf = D11sf*4*(m*pi/a)^3*(b/(n*pi))*(1-cos(pi*f/b));
funD1266sf = (2*D12sf+4*D66sf)*4*(m*pi/a)*(n*pi/b)*(1-cos(pi*f/b));
funD22sf = D22sf*4*(n*pi/b)^3*(a/(m*pi))*(1-cos(pi*f/b));

funD11 = funD11+funD11sf;
funD1266 = funD1266+funD1266sf;
funD22 = funD22+funD22sf;

for m_max = 1:mmax
    for n_max = 1:nmax
       funnew = funD11;
       funnew = subs(funnew,m,m_max);
       funnew = subs(funnew,n,n_max);
       fun(m_max,n_max) = funnew;
    end
end
funD11 = fun;

for m_max = 1:mmax
    for n_max = 1:nmax
       funnew = funD1266;
       funnew = subs(funnew,m,m_max);
       funnew = subs(funnew,n,n_max);
       fun(m_max,n_max) = funnew;    
    end
end
funD1266 = fun;

for m_max = 1:mmax
    for n_max = 1:nmax
       funnew = funD22;
       funnew = subs(funnew,m,m_max);
       funnew = subs(funnew,n,n_max);
       fun(m_max,n_max) = funnew;    
    end
end
funD22 = fun;
funLHS = funD11+funD1266+funD22;

for m_max = 1:mmax
    for n_max = 1:nmax
       funnew = funRHS;
       funnew = subs(funnew,m,m_max);
       funnew = subs(funnew,n,n_max);
       fun(m_max,n_max) = funnew;
    end
end

funRHS = fun;

%matriks stiffness skin
Xi = [ 0 -1 ];
XiD = [1 Xi(1) Xi(2) 0 0;...
       1 -Xi(1) Xi(2) 0 0;...
       0 0 -Xi(2) 1 0;...
       0 0 -Xi(2) 0 1;...
       0 0 0 0 0;...
       0 0 0 0 0];   
E11 = 135000; %MPa
E22 = 10000; %MPa
G12 = 5000; %MPa
v12 = 0.3; %Poisson'n Ratio
t_ply = 0.125; %tebal per ply, mm
h = 10*t_ply;
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
D = XiD*U*(h^3)/12;

%stiffness skin-flange
Xisf = Xi;
XiDsf =[1 Xisf(1) Xisf(2) 0 0;...
       1 -Xisf(1) Xisf(2) 0 0;...
       0 0 -Xisf(2) 1 0;...
       0 0 -Xisf(2) 0 1;...
       0 0 0 0 0;...
       0 0 0 0 0]; 
hsf = 2*h;
Dsf = XiDsf*U*(hsf^3)/12
   
panjang = 300;
lebar = 170;
flange = 45; %lebar setengah flange, max setengah lebar skin
load = 1;

% 
% funD1266 = subs(funD1266,D12sf,Dsf(3));
% funD1266 = subs(funD1266,D66sf,Dsf(4));
% funD1266 = subs(funD1266,a,panjang);
% funD1266 = subs(funD1266,b,lebar);
% funD1266 = subs(funD1266,f,flange);
% vpa(funD1266)

funLHS = subs(funLHS,D11,D(1));
funLHS = subs(funLHS,D12,D(3));
funLHS = subs(funLHS,D66,D(4));
funLHS = subs(funLHS,D22,D(2));
funLHS = subs(funLHS,D11sf,Dsf(1));
funLHS = subs(funLHS,D12sf,Dsf(3));
funLHS = subs(funLHS,D66sf,Dsf(4));
funLHS = subs(funLHS,D22sf,Dsf(2));
funLHS = subs(funLHS,a,panjang);
funLHS = subs(funLHS,b,lebar);
funLHS = subs(funLHS,f,flange)
funLHS = vpa(funLHS);

% % 
% 
funRHS = subs(funRHS,a,panjang);
funRHS = subs(funRHS,b,lebar);
funRHS = subs(funRHS,Nx,load);
funRHS = vpa(funRHS);

lambda = funLHS./funRHS
sort(lambda(:))
