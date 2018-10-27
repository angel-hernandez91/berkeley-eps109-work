function pool_check8(Z,NX,NY)

% input  Z(current elevation profile)
% output ZS, drain, pool, drainage
global ZS drain pool drainage;

%%%The matrices below have a border of zeros to make checking the boundary easier                                
% checked = zeros(NX,NY);   %matrix of checked points
pool     = zeros(NX,NY);    %matrix of pooled areas
drain    = zeros(NX,NY);    %matrix of draining points
drainage = zeros(NX,NY);    %matrix of drainage points(points connecting drains and pools)

ZS = zeros(NX*NY,3);            %Matrix of Z to be sorted.

%copy all NX*NY cell in 1D array ZS 
for j = 1:NY
    for i = 1:NX
        ZS(i+NX*(j-1),1) = Z(i,j);
        ZS(i+NX*(j-1),2) = i;
        ZS(i+NX*(j-1),3) = j;
    end
end

ZS = sortrows(ZS);               %sort the matrix with ascending elevations

%Set the main draining point, lowest points of all
x = ZS(1,2);
y = ZS(1,3);
drain(x,y) = 1;

P = 1;  %Pool number to identify pools and associated drainages (pool index)

%Loop over all cell starting the one above the lowest
for I = 2:NX*NY
    x = ZS(I,2);
    y = ZS(I,3);
    
    xU = mod(x-1-1,NX)+1;       % normally x-1 but observe p.b.c.
    xD = mod(x+1-1,NX)+1;       % normally x+1 but observe p.b.c.
    yL = mod(y-1-1,NY)+1;       % normally y-1 but observe p.b.c.
    yR = mod(y+1-1,NY)+1;       % normally y+1 but observe p.b.c.
    
    pL = pool(x,yL);            %Store values for surrounding pools and current maximum elevation of that pool
    pR = pool(x,yR);
    pU = pool(xU,y);
    pD = pool(xD,y);
    
    if pL & all(all(drainage~=pL));  %%%Establish if pool is present and  
        PL = pL;
    else
        PL = 0;
    end
    if pR & all(all(drainage~=pR));  %%%no drainage point already
        PR = pR;
    else
        PR = 0;
    end
    if pU & all(all(drainage~=pU));  %%%exists for that pool (ie. 
        PU = pU;
    else
        PU = 0;
    end
    if pD & all(all(drainage~=pD));  %%%non-draining pool)
        PD = pD;
    else
        PD = 0;
    end
    
    
    DL = pL & any(any(drainage==pL));  %%%Establish if pool is present and  
    DR = pR & any(any(drainage==pR));  %%%drainage point does already
    DU = pU & any(any(drainage==pU));  %%%exist for that pool (ie. 
    DD = pD & any(any(drainage==pD));  %%%draining pool)
    
    dL = drain(x,yL) | drainage(x,yL); %%%Establish if point contacts any
    dR = drain(x,yR) | drainage(x,yR); %%%draining areas or drainage points
    dU = drain(xU,y) | drainage(xU,y);
    dD = drain(xD,y) | drainage(xD,y);
    
    %%%%If connected to a pool that does not yet have a drainage point, the
    %%%%point must either be a pool or a point draining that pool (ie.
    %%%%drainage)
    if PL|PR|PU|PD
        if dL|dR|dU|dD|DL|DR|DU|DD
            drainage(x,y) = max([PL,PR,PU,PD]);
%             checked(x,y) = 1;
        else
            pool(x,y) = max([PL,PR,PU,PD]);
%             checked(x,y) = 1;
        end
        if PL
            pool(pool==PL) = max([PL,PR,PU,PD]);
        end
        if PR
            pool(pool==PR) = max([PL,PR,PU,PD]);
        end
        if PU
            pool(pool==PU) = max([PL,PR,PU,PD]);
        end
        if PD
            pool(pool==PD) = max([PL,PR,PU,PD]);
        end
          
    %%%%if connected to a pool with a drainage point it drains freely
    %%%%to that pool as the drainage point is at a lower elevation.
    %%%%Also, if the point is touching a drain or drainage point it
    %%%%also drains.
    elseif pL|pR|pU|pD|dL|dR|dU|dD
        drain(x,y) = 1;
%         checked(x,y) = 1;
        
    %%%%If the point doesn't contact a pool, drain, or drainage point,
    %%%%then it must be a new pool.
    else
        pool(x,y) = P;
%         checked(x,y) = 1;
        P = P+1; % increase pool index
    end
end 