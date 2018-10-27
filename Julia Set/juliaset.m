%Julia Set: We fix c, but vary z
%Mandlebrot Set: if instead we fix z = 0, and vary c we can generate the
                 %mandlebrot set

clf
clear

xmin = -1.5;
xmax = 1.5;
ymin = -1.5;
ymax = 1.5;

%number of points, resolution in x and y direction
nPoints = 1000;

%upper limit for |z|, consider series has diverged if exceeded
zmax = 1e6;

%interval initialization
x = linspace(xmin, xmax, nPoints);
y = linspace(ymin, ymax, nPoints);

%matrix to store z = z^2 + c values initialization
ZZ = zeros(nPoints, nPoints);
%the 'i' index walks along the 'x' direction
for i=1:nPoints
    % the 'j' index walks along the 'y' direction
    for j=1:nPoints 
        z = complex(x(i), y(j)); % this converts (x,y) to a complex number 'c'
        c = complex(-0.4, 0.6);  % iteration starts with fixed complex number
                                 % interesting choices: 
                                 % complex(-0.4, 0.6) 
                                 % complex(-0.85, 0.22)
        
        it_max = 100;            % increase iterations for more depth/resolution.
        for n = 1:it_max
            %iteration formula
            z = z^2 + c;
        
        %store the results in the matrix ZZ(nPoints,nPoints)
        if (abs(z)<=zmax) 
            ZZ(j,i)=0;  % make it '0' for blue if |z| is small and series has not diverged
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

%animated plot of blow up points
m = max(max(ZZ));
for i = 1:1000
    ZZ = mod(ZZ + 1, m);
    imagesc(x, y, ZZ);
    pause(0.2);
end