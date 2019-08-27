function y = seqsearch(seq)
global t_ply XiDopt lamparaflag point
%ex seq = [ 1 1 1 4 4 4 2 3 2 3 2 3 0 ]
%--> 24 ply symmetric (0 0 0 90 90 90 45 -45 45 -45 45 -45);
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

Xi_D(1:3) = 0;

%pendefinisian nilai z, skin
z(1:layernum+1) = 0;
z(1) = -0.5*t_ply*layernum;

for n = 2:layernum+1
    z(n) = z(n-1)+t_ply;
end

%perhitungan lamination parameters, skin
norm_z = 2*z'/(layernum*t_ply);
Xi_D(1:3) = 0;

for n = 1:layernum   
   Xi_D(1) = Xi_D(1)+ norm_z(n+1)^3*0.5*cosd(2*theta(n))-norm_z(n)^3*0.5*cosd(2*theta(n));
   Xi_D(2) = Xi_D(2)+ norm_z(n+1)^3*0.5*cosd(4*theta(n))-norm_z(n)^3*0.5*cosd(4*theta(n));
   Xi_D(3) = Xi_D(3)+ norm_z(n+1)^3*0.5*sind(2*theta(n))-norm_z(n)^3*0.5*sind(2*theta(n));
end

if lamparaflag == 2
    y = sum((XiDopt(point,:)-Xi_D(1:2)).^2);
else
    y = sum((XiDopt(point,:)-Xi_D).^2);
end

%sequence web
%pendefinisian nilai z, web
%perhitungan lamination parameters, web
