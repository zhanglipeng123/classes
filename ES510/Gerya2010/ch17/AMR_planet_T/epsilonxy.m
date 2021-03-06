% Function epsilonxy()
% This function compose EPSILONxy for specific node
% Function return EPSILONxy
function[epsxy,dvx,dvy]=epsilonxy(S,ni,nodx,nody,nodcel,celnod,celvar,celres,bupper,blower,bleft,bright,xsize,ysize)

%  ci1  vx1  ci3
%        |
%  vy1---+---vy2
%        |
%  ci2  vx2  ci4
% 
% Indexes for 4 cells
ci1=nodcel(ni,1); % upper left
ci2=nodcel(ni,2); % lower left
ci3=nodcel(ni,3); % upper right
ci4=nodcel(ni,4); % lower right
% Maximal distance to the left
ciL=ci1; 
cvL=3;
cnL=2;
if(ciL==0 || (ci2>0 && celres(ciL)>celres(ci2)))
    ciL=ci2;
    cvL=2;
    cnL=1;
end
% Maximal distance to the right
ciR=ci3; 
cvR=3;
cnR=4;
if(ciR==0 || (ci4>0 && celres(ciR)>celres(ci4)))
    ciR=ci4;
    cvR=2;
    cnR=3;
end
% Maximal distance to above
ciA=ci1; 
cvA=4;
cnA=3;
if(ciA==0 || (ci3>0 && celres(ciA)>celres(ci3)))
    ciA=ci3;
    cvA=1;
    cnA=1;
end
% Maximal distance to below
ciB=ci2; 
cvB=4;
cnB=4;
if(ciB==0 || (ci4>0 && celres(ciB)>celres(ci4)))
    ciB=ci4;
    cvB=1;
    cnB=2;
end
%
% dvy/dx
dvydx=0;
dvy=0;
% Left boundary
if(nodx(ni)==0 && ciR>0)
    % Free slip: dvy/dx=0
    if(bleft==0)
        % No slip: dvy/dx=2*vy2/dx
        dx=nodx(celnod(ciR,cnR))-nodx(ni);
        dvy=2*S(celvar(ciR,cvR));
        dvydx=dvy/dx;
    end
end
% Right boundary
if(nodx(ni)==xsize && ciL>0)
    % Free slip: dvy/dx=0
    if(bright==0)
        % No slip: dvy/dx=-2*vy1/dx
        dx=nodx(ni)-nodx(celnod(ciL,cnL));
        dvy=-2*S(celvar(ciL,cvL));
        dvydx=dvy/dx;
    end
end
% Internal nodes
if(nodx(ni)>0 && nodx(ni)<xsize && ciL>0 && ciR>0)
    dx=(nodx(celnod(ciR,cnR))-nodx(celnod(ciL,cnL)))/2;
    dvy=(S(celvar(ciR,cvR))-S(celvar(ciL,cvL)));
    dvydx=dvy/dx;
end
%
% dvx/dy
dvxdy=0;
dvx=0;
% Upper boundary
if(nody(ni)==0 && ciB>0)
    % Free slip: dvx/dy=0
    if(bupper==0)
        % No slip: dvx/dy=2*vx2/dy
        dy=nody(celnod(ciB,cnB))-nody(ni);
        dvx=2*S(celvar(ciB,cvB));
        dvxdy=dvx/dy;
    end
end
% Lower boundary
if(nody(ni)==ysize && ciA>0)
    % Free slip: dvx/dy=0
    if(blower==0)
        % No slip: dvx/dy=-2*vx1/dy
        dy=nody(ni)-nody(celnod(ciA,cnA));
        dvx=-2*S(celvar(ciA,cvA));
        dvxdy=dvx/dy;
    end
end
% Internal nodes
if(nody(ni)>0 && nody(ni)<ysize && ciB>0 && ciA>0)
        dy=(nody(celnod(ciB,cnB))-nody(celnod(ciA,cnA)))/2;
        dvx=(S(celvar(ciB,cvB))-S(celvar(ciA,cvA)));
        dvxdy=dvx/dy;
end
% Compute EPSILONxy
epsxy=(dvxdy+dvydx)/2;
dvx=abs(dvx);
dvy=abs(dvy);
end

