function y = seqsearch(seq)
global t_ply XiAopt XiDopt flag
%flag 1 = untuk skin, flag 2 = untuk web
layernum = 2*(length(seq)-1)-seq(length(seq));
halflayer = ceil(layernum/2);

%Convert to degree
 for n = 1:halflayer
    if seq(n) == 1
        seq(n) = 0;
    end
    if seq(n) == 2
        seq(n) = 45;
    end
    if seq(n) == 3
        seq(n) = -45;
    end
    if seq(n) == 4
        seq(n) = 90;
    end
 end

theta(1:layernum) = 0;

for n = 1:halflayer
    theta(n) = seq(n);
end

m = 1;
for n = layernum:-1:halflayer+1
    theta(n) = seq(m);
    m = m+1;
end

%pendefinisian nilai z, skin
z(1:layernum+1) = 0;
z(1) = -0.5*t_ply*layernum;

for n = 2:layernum+1
    z(n) = z(n-1)+t_ply;
end

%perhitungan lamination parameters
norm_z = 2*z'/(layernum*t_ply);
Xi_A(1:3) = 0;
Xi_D(1:3) = 0;


for n = 1:layernum   
   Xi_A(1) = Xi_A(1)+ norm_z(n+1)*0.5*cosd(2*theta(n))-norm_z(n)*0.5*cosd(2*theta(n));
   Xi_A(2) = Xi_A(2)+ norm_z(n+1)*0.5*cosd(4*theta(n))-norm_z(n)*0.5*cosd(4*theta(n));
   Xi_A(3) = Xi_A(3)+ norm_z(n+1)*0.5*sind(2*theta(n))-norm_z(n)*0.5*sind(2*theta(n));
   Xi_D(1) = Xi_D(1)+ norm_z(n+1)^3*0.5*cosd(2*theta(n))-norm_z(n)^3*0.5*cosd(2*theta(n));
   Xi_D(2) = Xi_D(2)+ norm_z(n+1)^3*0.5*cosd(4*theta(n))-norm_z(n)^3*0.5*cosd(4*theta(n));
   Xi_D(3) = Xi_D(3)+ norm_z(n+1)^3*0.5*sind(2*theta(n))-norm_z(n)^3*0.5*sind(2*theta(n));
end

if flag == 2
    theta(layernum+1:layernum*2) = theta(1:layernum);    
    layernum = layernum*2;
    clear z norm_z
    z(1:layernum+1) = 0;
    z(1) = -0.5*t_ply*layernum;

    for n = 2:layernum+1
        z(n) = z(n-1)+t_ply;
    end

    norm_z = 2*z'/(layernum*t_ply);
    Xi_Dw(1:3) = 0;
    
    for n = 1:layernum   
       Xi_Dw(1) = Xi_Dw(1)+ norm_z(n+1)^3*0.5*cosd(2*theta(n))-norm_z(n)^3*0.5*cosd(2*theta(n));
       Xi_Dw(2) = Xi_Dw(2)+ norm_z(n+1)^3*0.5*cosd(4*theta(n))-norm_z(n)^3*0.5*cosd(4*theta(n));
       Xi_Dw(3) = Xi_Dw(3)+ norm_z(n+1)^3*0.5*sind(2*theta(n))-norm_z(n)^3*0.5*sind(2*theta(n));
    end
    
end

if flag == 2
   y = sum((XiDopt(3,:)-Xi_D).^2)+10*sum((XiDopt(1,:)-Xi_Dw).^2)+sum((XiAopt(2,1:2-Xi_A(1:2)).^2))+100*Xi_A(3)^2; 
else
   y = sum((XiDopt(2,:)-Xi_D).^2)+sum((XiAopt(1,2)-Xi_A(1:2)).^2+100*Xi_A(3)^2);
end

