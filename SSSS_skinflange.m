function y = SSSS_skinflange(Xi)
global U a b f t_ply Nx maxtweb total_plyweb

%SSSS Boundary, 0 +-45 90 layer only
%Matriks Xi D 
%skin
XiDs = [1 Xi(1) Xi(2) 0 0;...
       1 -Xi(1) Xi(2) 0 0;...
       0 0 -Xi(2) 1 0;...
       0 0 -Xi(2) 0 1;...
       0 0.5*Xi(3) 0 0 0;...
       0 0.5*Xi(3) 0 0 0]; 
%flange
XiDf = [1 Xi(4) Xi(5) 0 0;...
       1 -Xi(4) Xi(5) 0 0;...
       0 0 -Xi(5) 1 0;...
       0 0 -Xi(5) 0 1;...
       0 0.5*Xi(6) 0 0 0;...
       0 0.5*Xi(6) 0 0 0]; 
   
%Matriks Xi A
%skin
XiAs = [1 Xi(7) Xi(8) 0 0;...
       1 -Xi(7) Xi(8) 0 0;...
       0 0 -Xi(8) 1 0;...
       0 0 -Xi(8) 0 1;...
       0 0.5*Xi(9) 0 0 0;...
       0 0.5*Xi(9) 0 0 0]; 
%flange
XiAf = [1 Xi(10) Xi(11) 0 0;...
       1 -Xi(10) Xi(11) 0 0;...
       0 0 -Xi(11) 1 0;...
       0 0 -Xi(11) 0 1;...
       0 0.5*Xi(12) 0 0 0;...
       0 0.5*Xi(12) 0 0 0]; 

%tebal flange --> 1/2 tebal web
hf = 0.5*total_plyweb*t_ply;

%matriks kekakuan flange
Df = (XiDf*U*hf^3)/12;
Af = XiAf*U*hf; 

%definisi Ds Df, As Af ke matriks kekakuan masing-masing
% (D(1,1)) -> D(1)
% (D(2,2)) -> D(2)
% (D(1,2)) -> D(3)
% (D(3,3)) -> D(4)
% (D(1,3)) -> D(5)
% (D(2,3)) -> D(6)

%matriks kekakuan skin
for h = 12*t_ply:t_ply:maxtweb   
   Ds = (XiDs*U*h^3)/12; 
   As = XiAs*U*h; 
   
%mencari matriks kekakuan tumpukan skin+flange
   Asf = [ As(1) As(3) As(5); As(3) As(2) As(6); As(5) As(6) As(4)] + [Af(1) Af(3) Af(5); Af(3) Af(2) Af(6); Af(5) Af(6) Af(4)];
   Bsf = hf*0.5*[As(1) As(3) As(5); As(3) As(2) As(6); As(5) As(6) As(4)]-h*0.5*[Af(1) Af(3) Af(5); Af(3) Af(2) Af(6); Af(5) Af(6) Af(4)];
   Dsf = [ Ds(1) Ds(3) Ds(5); Ds(3) Ds(2) Ds(6); Ds(5) Ds(6) Ds(4)] + [ Df(1) Df(3) Df(5); Df(3) Df(2) Df(6); Df(5) Df(6) Df(4)]...
         + (hf^2)*0.25*[ As(1) As(3) As(5); As(3) As(2) As(6); As(5) As(6) As(4)] + (h^2)*0.25*[Af(1) Af(3) Af(5); Af(3) Af(2) Af(6); Af(5) Af(6) Af(4)];
%Reduced bending stiffness skin+flange
   Dsmea = Dsf-(Bsf*inv(Asf)*Bsf);
   clear Dsf
   Dsf(1) = Dsmea(1,1);
   Dsf(2) = Dsmea(2,2);
   Dsf(3) = Dsmea(1,2);
   Dsf(4) = Dsmea(3,3);
   Dsf(5) = Dsmea(1,3);
   Dsf(6) = Dsmea(2,3);
   Dsf = Dsf';

   term = 7;
   m = term;
   n = m;
   baris = 1;
   kolom = 1;
   LHS(m*n,m*n) = 0;
   RHS(m*n,m*n) = 0; 
   
   for K = 1:m
      a_K = K*pi/a;
      for L = 1:n
         b_L = L*pi/b;
         for i_=1:m
            for J=1:n
               b_J = J*pi/b;
               if((i_==K)&&(J~=L))
                  fun_sf1 = (Dsf(1)*(a_K^4)+Dsf(2)*(b_J^2)*(b_L^2)+Dsf(3)*((a_K^2)*(b_L^2)+(a_K^2)*(b_J^2)))*((1/(b_J-b_L))*sin(b_J*f-b_L*f)-(1/(b_J+b_L))*sin(b_J*f+b_L*f));
                  fun_sf2 = 4*Dsf(4)*(a_K^2)*(b_L*b_J)*((1/(b_J-b_L))*sin(b_J*f-b_L*f)+(1/(b_J+b_L))*sin(b_J*f+b_L*f));
                  fun_sf3 = (Dsf(1)*(a_K^4)+Dsf(2)*(b_J^2)*(b_L^2)+Dsf(3)*((a_K^2)*(b_L^2)+(a_K^2)*(b_J^2)))*((1/(b_J-b_L))*(sin((b_J-b_L)*b)-sin((b_J-b_L)*(b-f)))-((1/(b_J+b_L))*(-sin((b_J+b_L)*(b-f)))));
                  fun_sf4 = 4*Dsf(4)*(a_K^2)*(b_L*b_J)*((1/(b_J-b_L))*(sin((b_J-b_L)*b)-sin((b_J-b_L)*(b-f)))+((1/(b_J+b_L))*(-sin((b_J+b_L)*(b-f)))));
                  fun_sfl = (a/8)*(fun_sf1+fun_sf2);
                  fun_sfr = (a/8)*(fun_sf3+fun_sf4);
                  fun_s1 = (Ds(1)*(a_K^4)+Ds(2)*(b_J^2)*(b_L^2)+Ds(3)*((a_K^2)*(b_L^2)+(a_K^2)*(b_J^2)))*((1/(b_J-b_L))*(sin((b_J-b_L)*(b-f))-sin(b_J*f-b_L*f))-(1/(b_J+b_L))*(sin((b_J+b_L)*(b-f))-sin(b_J*f+b_L*f)));
                  fun_s2 = 4*Ds(4)*(a_K^2)*(b_L*b_J)*((1/(b_J-b_L))*(sin((b_J-b_L)*(b-f))-sin(b_J*f-b_L*f))+(1/(b_J+b_L))*(sin((b_J+b_L)*(b-f))-sin(b_J*f+b_L*f)));
                  fun_s = (a/8)*(fun_s1+fun_s2);
                  LHS(baris,kolom) = fun_sfl+fun_sfr+fun_s;
               elseif((i_==K)&&(J==L))
                  fun_sf1 = (Dsf(1)*(a_K^4)+Dsf(2)*(b_L^4)+2*Dsf(3)*(a_K^2)*(b_L^2))*(f-((1/(2*b_L))*sin(2*b_L*f)));
                  fun_sf2 = 4*Dsf(4)*(a_K^2)*(b_L^2)*(f+((1/(2*b_L))*sin(2*b_L*f)));
                  fun_sf3 = (Dsf(1)*(a_K^4)+Dsf(2)*(b_L^4)+2*Dsf(3)*(a_K^2)*(b_L^2))*(f+((1/(2*b_L))*sin(2*b_L*(b-f))));
                  fun_sf4 = 4*Dsf(4)*(a_K^2)*(b_L^2)*(f-((1/(2*b_L))*sin(2*b_L*(b-f))));
                  fun_sfl = (a/8)*(fun_sf1+fun_sf2);
                  fun_sfr = (a/8)*(fun_sf3+fun_sf4);
                  fun_s1 = (Ds(1)*(a_K^4)+Ds(2)*(b_L^4)+2*Ds(3)*(a_K^2)*(b_L^2))*((b-2*f)-((1/(2*b_L))*(sin(2*b_L*(b-f))-sin(2*b_L*f))));
                  fun_s2 = 4*Ds(4)*(a_K^2)*(b_L^2)*((b-2*f)+((1/(2*b_L))*(sin(2*b_L*(b-f))-sin(2*b_L*f))));
                  fun_s = (a/8)*(fun_s1+fun_s2);
                  LHS(baris,kolom) = fun_sfl+fun_sfr+fun_s;
                  RHS(baris,kolom) = -0.5*Nx*(a_K^2)*(a*b/4);
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
   
   lambda = eig(LHS,RHS);
   if max(lambda) <= -1
       RF = 1+(max(lambda));
       break
   end
end

y = (h/t_ply)+(RF/10);
%yang beda --> i =k, j=l, fun_sf2