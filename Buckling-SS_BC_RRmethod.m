clc
clear
syms x y a b m n PI D11 D12 D66 D22 Nx 
%Hanya untuk tebal dan D konstant sepanjang pelat
terms = 6; %number of terms
mmax = terms;
nmax = mmax;

% funD11 = ((m*PI/a)^4)*sin(m*PI*x/a)*sin(n*PI*y/b);
% funD11 = int(funD11,x,0,a);
% funD11 = int(funD11,y,0,b)
% 
% funD1266 = ((m*PI/a)^2)*((n*PI/b)^2)*sin(m*PI*x/a)*sin(n*PI*y/b);
% funD1266 = int(funD1266,x,0,a);
% funD1266 = int(funD1266,y,0,b)
% 
% funD22 = ((n*PI/b)^4)*sin(m*PI*x/a)*sin(n*PI*y/b);
% funD22 = int(funD22,x,0,a);
% funD22 = int(funD22,y,0,b)
% 
% funRHS = ((m*PI/a)^2)*sin(m*PI*x/a)*sin(n*PI*y/b);
% funRHS = int(funRHS,x,0,a);
% funRHS = int(funRHS,y,0,b)

% funD1266 = subs(funD1266,PI*m,PI*(2*m-1));
% funD1266 = subs(funD1266,PI*n,PI*(2*n-1));
% funD22 = subs(funD22,PI*m,PI*(2*m-1));
% funD22 = subs(funD22,PI*n,PI*(2*n-1));
% funRHS = subs(funRHS,PI*m,PI*(2*m-1));
% funRHS = subs(funRHS,PI*n,PI*(2*n-1));
% 

funD11 = 4*(m*pi/a)^3*(b/(n*pi));
funD1266 = 4*(m*pi/a)*(n*pi/b);
funD22 = 4*(n*pi/b)^3*(a/(m*pi));
funRHS = 4*b*m/(a*n);

for m_max = 1:mmax
    for n_max = 1:nmax
       funnew = funD11;
       funnew = subs(funnew,m,m_max);
       funnew = subs(funnew,n,n_max);
       fun(m_max,n_max) = funnew;
       
    end
end
funD11 = D11*fun;

for m_max = 1:mmax
    for n_max = 1:nmax
       funnew = funD1266;
       funnew = subs(funnew,m,m_max);
       funnew = subs(funnew,n,n_max);
       fun(m_max,n_max) = funnew;    
    end
end
funD1266 = (2*D12+4*D66)*fun;

for m_max = 1:mmax
    for n_max = 1:nmax
       funnew = funD22;
       funnew = subs(funnew,m,m_max);
       funnew = subs(funnew,n,n_max);
       fun(m_max,n_max) = funnew;    
    end
end
funD22 = D22*fun;
funLHS = funD11+funD1266+funD22;

for m_max = 1:mmax
    for n_max = 1:nmax
       funnew = funRHS;
       funnew = subs(funnew,m,m_max);
       funnew = subs(funnew,n,n_max);
       fun(m_max,n_max) = funnew;
    end
end

funRHS = Nx*fun;

%matriks Xi  
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
panjang = 170;
lebar = 30;
load = 1;

funLHS = subs(funLHS,D11,D(1));
funLHS = subs(funLHS,D12,D(3));
funLHS = subs(funLHS,D66,D(4));
funLHS = subs(funLHS,D22,D(2));
funLHS = subs(funLHS,a,panjang);
funLHS = subs(funLHS,b,lebar);
funLHS = vpa(funLHS);

% % 
% 
funRHS = subs(funRHS,a,panjang);
funRHS = subs(funRHS,b,lebar);
funRHS = subs(funRHS,Nx,load);
funRHS = vpa(funRHS);

lambda = funLHS./funRHS;
sort(lambda(:))
