%Mandelbrot set

clf
clear

xmin = -2;
xmax = 2;
ymin = -2;
ymax = 2;

%number of points = resolution in x and y direction
nPoints = 5000;

%upper limit for |z|, consider series has diverged if exceeded
zmax = 1e6;
x = linspace(xmin, xmax, nPoints);
y = linspace(ymin, ymax, nPoints);
%the 'i' index walks along the 'x' direction
for i=1:nPoints
    % the 'j' index walks along the 'y' direction
    for j=1:nPoints 

        %%%% map the index 'i' onto the interval [xmin,xmax]
        %%%% i=1       must yield 'xmin'
        %%%% i=nPoints must yield 'xmax'
        %%%% map the index 'j' onto the interval [ymin,ymax]
        %%%% j=1       must yield 'ymin'
        %%%% j=nPoints must yield 'ymax'
   

        z = complex(x(i),y(j)); % this converts (x,y) to a complex number 'c'
        c = complex(-0.23, 0.8);   % iteration starts with complex number = 0
        
        it_max = 100000;            % no more than 50 iterations please
        for n = 1:it_max
            %%%% add the iteration formula for the Mandelbrot set
            z = z^2 + c;
        
        %store the results in the matrix ZZ(nPoints,nPoints)
        if (abs(z)<=zmax) ZZ(j,i)=0;  % make it '0' for blue if |z| is small and series has not diverged
        elseif (abs(z) > zmax) 
            ZZ(j, i) = n; % make it '1' for red  if |z| is large and series has diverged
            break
        end
        end
    end
    fprintf('i= %d of %d points in X direction: x= %f\n',i,nPoints,x(i));
end

imagesc(x,y,ZZ); % plot the image
set(gca,'YDir','normal');

m = max(max(ZZ));
for i = 1:1000
    ZZ = mod(ZZ + 1, m);
    imagesc(x, y, ZZ);
    set (gca, 'YDir', 'normal');
    pause(0.2);
end