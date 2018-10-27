clear
clf

NX = 50; %number of rows
NY = 70; %number of columns

d  = 5; % grid spacing in meters
dx = d; % keep dx=dy for simplicity
dy = d;

LX=NX*dx;
LY=NY*dy;

Z = zeros(NX,NY);
for i = 1:NX*NY
    Z(i) = rand();
end

dt = 50;  %dx^2/(4*D);

%Area exponent A^m, standard m=1
m=1;

%gradient exponent g^n, standard n=1
n=1;

%threshold 
theta_c = 1; 

K = 1e-6; % meters^(1-2m)/yr, standard 10^-6

D = 0.005; % m^2/yr

x = 1:NX;
y = 1:NY;
[X,Y] = meshgrid(y,x); %strange that y goes first !!!

global A ZS drain drainage pool;
pool_check8(Z,NX,NY)

%pne = not-a-pool(NX,NY)
pne = zeros(NX,NY);
pne(pool==0) = 1;

A = zeros(NX,NY);

Znew = Z;
jt = 1;
for it = 1:2000*(600/dt)
    calculate_collection_area2(Z,NX,NY);
    A = A.*(dx*dy);
    
    for i = 2:NX
        iL = mod(i-1-1,NX)+1; % normally i-1 but observe p.b.c.
        iR = mod(i+1-1,NX)+1; % normally i+1 but observe p.b.c.

        for j = 2:NY
            jD = mod(j-1-1,NY)+1; % normally j-1 but observe p.b.c.
            jU = mod(j+1-1,NY)+1; % normally j+1 but observe p.b.c.
  
            if drain(i,j)>0 %this cell is a drain
                s1 = (Z(i +1, j) - Z(i -1, j))/(2*dx);
                s2 = (Z(i, j + 1) - Z(i, j - 1))/(2*dy);
                s3 = (Z(i + 1, j -1) - Z(i - 1, j + 1))/(2*sqrt(dx^2 + dy^2));
                s4 = (Z(i + 1, j + 1) - Z(i -1, j -1))/(2*sqrt(dx^2 + dy^2));
                gradient = (sqrt(s1^2 + s2^2) + sqrt(s3^2 + s4^2))/2;
                
                Gz = D*((Z(i +1, j) - 2*Z(i, j) + Z(i - 1, j))/dx^2 + (Z(i, j + 1) - 2*Z(i, j) + Z(i, j + 1))/dy^2) - K*(A.^m*gradient.^n - theta_c);              

            elseif drainage(i,j)>0 %this cell is a drainage point (it drains a pool)
                
                if (Z(i,j)>=Z(iR,j)) && pool(iR,j)~=drainage(i,j) 
                    gradient = (Z(i,j)-Z(iR,j))/dx; %pool is on my left, I drain to the right, use this gradiant
                elseif (Z(i,j)>=Z(iL,j)) && pool(iL,j)~=drainage(i,j)
                    gradient = (Z(i,j)-Z(iL,j))/dx;
                elseif (Z(i,j)>=Z(i,jU)) && pool(i,jU)~=drainage(i,j)
                    gradient = (Z(i,j)-Z(i,jU))/dy;
                elseif (Z(i,j)>=Z(i,jD)) && pool(i,jD)~=drainage(i,j)
                    gradient = (Z(i,j)-Z(i,jD))/dy;
                else
                    error
%                    gradient = 0.02; % ??? This does happen (maybe when two pools merge)
                end
                Gz = D*((Z(i +1, j) - 2*Z(i, j) + Z(i - 1, j))/dx^2 + (Z(i, j + 1) - 2*Z(i, j) + Z(i, j + 1))/dy^2) - K*(A.^m*gradient.^n - theta_c);         
            else %this cell is a pool, assume it has some mass diffusion but no erosion!
                Gz = 0;
            end
            if (Gz<0) Gz = 0; end

            % diffusion term
            phiz = D*((Z(i +1, j) - 2*Z(i, j) + Z(i - 1, j))/dx^2 + (Z(i, j + 1) - 2*Z(i, j) + Z(i, j + 1))/dy^2);
           
            Znew(i,j) = Z(i,j) + dt*(phiz - Gz) + 0.3*(dt/600);  % 0.3 is some uplift rate
            
            if (Znew(i,j)<0)
                Znew(i,j) = 0;  % ???? yes this does happen at the boundary when kept at zero
            end            
        end
    end
    
    Znew(1,:) = 0;
%     Znew(ZS(1,2),ZS(1,3)) = Z(ZS(1,2),ZS(1,3));
    Z = Znew;
    
    pool_check8(Z,NX,NY)
    pne = zeros(NX,NY);
    pne(pool==0) = 1;
    
    if (mod(it,round(jt))==0)
        if it>jt 
            jt = jt*1.1;
        end
        if pool(:,:) pool=1; end
        fprintf('Iteration %d\n',it);

        subplot(1,2,1);
       
        pcolor(X,Y,Z);
        surf(X,Y,Z);
        axis([0 NY 0 NX 0 140])
        title('Surface Relief');
        xlabel('x');
        ylabel('y');
        view([+34.5 14]);

        subplot(1,2,2);

        hold off
        pcolor(X,Y,Z)
        shading interp
        hold on
        contour(X,Y,Z,30,'k');  
        drawnow;
        
%        saveas(gcf, ['erosion_valley_',sprintf('%03d',f),'.png'])
%        f = f+1;
    end
end % end of main loop over iterations

%final plot
subplot(1,2,1);
hold off
surf(X,Y,Z);
title('Surface Relief')
xlabel('x');
ylabel('y');
view([+34.5 14]);

subplot(1,2,2);
hold off
pcolor(X,Y,Z)
shading interp
hold on
contour(X,Y,Z,30,'k');  

fprintf('Simulation finished.\n\n');
