clear
clf

global A ZS drain drainage pool ocean ZBeachLevel AOcean;

NX = 1.5*50; %1st index  = row index
NY = 1.5*70; %2nd index  = column index in matrix below

d  = 2; % grid spacing in meters
dx = d; % keep dx=dy for simplicity
dy = d;

LX=NX*dx;
LY=NY*dy;

%set initial conditions
Z = zeros(NX,NY);
for i=1:4
   Z = initial_conditions_add_hill(Z,NX,NY,rand()*NX,rand()*NY,rand()*NX,rand()*70);
end
for i=1:4
   Z = initial_conditions_add_hill(Z,NX,NY,rand()*NX,rand()*NY,rand()*NX/5,rand()*50);
end

%set the ocean volume
oceanLevelParameter = 0; % what does this parameter do?
Zmin = min(min(Z));
Zmax = max(max(Z));
ZBeachLevel = Zmin+oceanLevelParameter*(Zmax-Zmin); 
VOcean=0.0;
AOcean=0;
for i = 1:NX
    for j = 1:NY
        if (Z(i,j)<=ZBeachLevel)
            AOcean = AOcean+1;
            VOcean = VOcean+ZBeachLevel-Z(i,j);
        end
    end
end
fprintf('Minimum elevation      %f\n',Zmin);
fprintf('Maximum elevation      %f\n',Zmax);
fprintf('Beach level            %f\n',ZBeachLevel);
fprintf('Ocean volume           %f\n',VOcean);
fprintf('Ocean surface fraction %f %%\n\n',AOcean/(NX*NY)*100);

%time step
dt = 50;  % dt<dx^2/(4*D) not just dx^2/(2*D)

%Area exponent, setting m=1
m=1;

%gradient exponent, setting n=1
n=1;

%threshold 
theta_c = 1;

K = 1e-6; % meters^(1-2m)/yr

D = 0.0005; % m^2/yr
% D = 0.0;

% any uplift?
U = 0.0;

eta1 = D*dt/dx/dx;
eta2 = D*dt/dy/dy;
if (eta1>=0.25 || eta2>=0.25)
    fprintf('Check if both eta values are below 0.25: %f %f\n',eta1,eta2); 
    pause
end

fprintf('erosion_lab8: Stardard simulation n=%.2f m=%.2f and small dt=%.0f\n\n',n,m,dt);

x = 1:NX;
y = 1:NY;
[X,Y] = meshgrid(y,x); %strange the y goes first !!!

pool_check10(Z,NX,NY,VOcean)

subplot(1,2,1);
pcolor(X,Y,Z);
surf(X,Y,Z);
axis([0 NY 0 NX 0 140])
title('Surface Relief');
xlabel('x');
ylabel('y');
view([+34.5 14]);

subplot(1,2,2);
pne = zeros(NX,NY); %pne = not-a-pool(NX,NY)
pne(pool==0) = 1;
pcolor(X,Y,pne);
colormap jet;
title('Pool Locations');
xlabel('x');
ylabel('y');
set(figure(1), 'position', [10, 200, 1200, 400], 'PaperPositionMode', 'auto','color', 'white')
drawnow;

A = zeros(NX,NY);

Znew = Z;
jt = 1;
f = 1;
t = 0; %time in years
for it = 1:2000*(600/dt)
    t=t+dt;
    
    calculate_collection_area2(Z,NX,NY);
    A = A.*(dx*dy);
    
    for i = 1:NX
        iL = mod(i-1-1,NX)+1; % normally i-1 but observe p.b.c.
        iR = mod(i+1-1,NX)+1; % normally i+1 but observe p.b.c.

        for j = 1:NY
            jD = mod(j-1-1,NY)+1; % normally j-1 but observe p.b.c.
            jU = mod(j+1-1,NY)+1; % normally j+1 but observe p.b.c.

            uplift = U;
            if ocean(i,j)
                Gz     = 0;
                phiz   = 0;
                uplift = 0;
            elseif drain(i,j)>0 %this cell is a drain
                s1 = (Z(iR,j)  - Z(iL,j))/(2*dx);
                s2 = (Z(i,jU)  - Z(i,jD))/(2*dy);
                s3 = (Z(iR,jD) - Z(iL,jU))/(2*sqrt(dx^2+dy^2));
                s4 = (Z(iR,jU) - Z(iL,jD))/(2*sqrt(dx^2+dy^2));
                gradient = ( sqrt(s1^2 + s2^2) + sqrt(s3^2 + s4^2) )/2;
                
                Gz = K*( A(i,j)^m*gradient^n - theta_c );                  
                % diffusion term
                phiz = D*((Z(iL,j) - 2*Z(i,j) + Z(iR,j))/(dx*dx) + (Z(i,jU) - 2*Z(i,j) + Z(i,jD))/(dy*dy));

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
                    error('error',i,j); % should not happen
                end
                
                Gz = K*( A(i,j)^m*gradient^n - theta_c);                
                % diffusion term as before
                phiz = D*((Z(iL,j) - 2*Z(i,j) + Z(iR,j))/(dx*dx) + (Z(i,jU) - 2*Z(i,j) + Z(i,jD))/(dy*dy));

            else %this cell is a pool, assume it has the regular mass diffusion but no erosion!
                Gz = 0;
                phiz = D*((Z(iL,j) - 2*Z(i,j) + Z(iR,j))/(dx*dx) + (Z(i,jU) - 2*Z(i,j) + Z(i,jD))/(dy*dy));
            end
            
            if (Gz<0) Gz = 0; end % no erosion if A^m * g^n < theta_c

            Znew(i,j) = Z(i,j) + dt*(phiz - Gz) + uplift*(dt/600);  % add some uplift
            
            if (Znew(i,j)<0)
                Znew(i,j) = 0;  % ???? yes this does happen at the boundary when kept at zero
            end            
        end
    end
    
    Z = Znew;
    pool_check10(Z,NX,NY,VOcean)
    if (mod(it,round(jt))==0)
        if it>jt 
            jt = jt*1.1;
        end
        if pool(:,:) pool=1; end
        fprintf('Iteration=%5d  time= %6d  ocean level=%.4f   ocean surface fraction=%.3f\n',it,t,ZBeachLevel,100*AOcean/(NX*NY));

        subplot(1,2,1);
       
        pcolor(X,Y,Z);
        surf(X,Y,Z);
        axis([0 NY 0 NX Zmin Zmax])
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
    end
    if mod(f, 200) == 0
    
        saveas(gcf, ['heat_', sprintf('%03d', f), '.png'])
    end
    f = f + 1;
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
