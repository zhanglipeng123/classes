% Function Stokes_Continuity_solver_ghost()
% This function formulates and solves  
% Stokes and Continuity equations defined on 2D staggered grid
% with specified resolution (xnum, ynum) and gridsteps (xstp, ystp)
% given distribution of right parts for all equations (RX,RY,RC) on the grid 
% and given variable shear (etas) and normal (etan) viscosity distributions 
% pressure is normalized relative to given value (prnorm) in the first cell
%
% Velocity Boundary condition specified by bleft,bright,btop,bbot 
% (1=free slip -1=no slip) are implemented from ghost nodes 
% directly into Stokes and continuity equations
%
% Function returns solution for velocity and pressure (vx,vy,pr)
% and distribution of residuals (resx,resy,resc)
function[vx,resx,vy,resy,pr,resc]=Stokes_Continuity_solver_channel(prnorm,etas,etan,xnum,ynum,xstp,ystp,RX,RY,RC,bleft,bright,btop,bbottom)
% 
% Staggered Grid for Multigrid
% 
%     vx       vx       vx    
%
% vy  +---vy---+---vy---+   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  +---vy---+---vy---+   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  +---vy---+---vy---+   vy
%
%     vx       vx       vx    
% 
% Lines show basic grid
% Basic (density) nodes are shown with +
% Ghost nodes shown outside the basic grid
% are used for boundary conditions

% Boundary conditions
% Pressure boundary condition
bpres=0;
% Upper boundary
% Free slip, No slip
btopx=btop;
btopy=0;
% Open boundary: dvx/dy=0, dvy/dy=0
if(btop==0)
    btopx=1;
    btopy=1;
    % Channel flow top->bottom
    bpres=1;
end
% Lower boundary
% Free slip, No slip
bbottomx=bbottom;
bbottomy=0;
% Open boundary: dvx/dy=0, dvy/dy=0
if(bbottom==0)
    bbottomx=1;
    bbottomy=1;
end
% Left Boundary
% Free slip, No slip
bleftx=0;
blefty=bleft;
% Open boundary: dvx/dx=0, dvy/dx=0
if(bleft==0)
    bleftx=1;
    blefty=1;
    % Channel flow left->right
    bpres=2;
end
% Right Boundary
% Free slip, No slip
brightx=0;
brighty=bright;
% Open boundary: dvx/dx=0, dvy/dx=0
if(bright==0)
    brightx=1;
    brighty=1;
end



% Poisson-like equations koefficients
xkf=1/xstp^2;
xkf2=2/xstp^2;
ykf=1/ystp^2;
ykf2=2/ystp^2;
xykf=1/(xstp*ystp);

% Koefficient for scaling pressure
pscale=2*etan(1)/(xstp+ystp);

% Horizontal shift index
ynum3=(ynum-1)*3;


% Creating matrix
L=sparse((xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3);
R=zeros((xnum-1)*(ynum-1)*3,1);

% Solving of Stokes and continuity equations on nodes
for i=1:1:ynum-1
    for j=1:1:xnum-1
        % Indexes for P,vx,vy
        ivx=((j-1)*(ynum-1)+(i-1))*3+1;
        ivy=ivx+1;
        ipr=ivx+2;
        
        % x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
        if (j<xnum-1)
            % x-Stokes equation stensil
            %     +-------------------- -+----------------------+   
            %     |                      |                      |
            %     |                      |                      |
            %     |                   vx(i-1,j)                 |    
            %     |                      |                      |
            %     |                      |                      |
            %     +-----vy(i-1,j)---etas(i,j+1)---vy(i-1,j+1)---+
            %     |                      |                      |
            %     |                      |                      |
            % vx(i,j-1)  pr(i,j)      vx(i,j)     P(i,j+1)   vx(i,j+1)    
            %     |     etan(i,j)        |       etan(i,j+1)    |
            %     |                      |                      |
            %     +------vy(i,j)---etas(i+1,j+1)---vy(i,j+1)----+
            %     |                      |                      |
            %     |                      |                      |
            %     |                   vx(i+1,j)                 |    
            %     |                      |                      |
            %     |                      |                      |
            %     +-------------------- -+----------------------+   
            % Right Part
            R(ivx,1)=RX(i+1,j+1);
            % Computing Current x-Stokes coefficients
            % Central Vx node
            L(ivx,ivx)=-xkf2*(etan(i,j+1)+etan(i,j))-ykf*(etas(i+1,j+1)+etas(i,j+1));
            % Left Vx node
            if (j>1)
                ivxleft=ivx-ynum3;
                L(ivx,ivxleft)=xkf2*etan(i,j);
            else
                L(ivx,ivx)=L(ivx,ivx)+bleftx*xkf2*etan(i,j);
            end
            % Right Vx node
            if (j<xnum-2)
                ivxright=ivx+ynum3;
                L(ivx,ivxright)=xkf2*etan(i,j+1);
            else
                L(ivx,ivx)=L(ivx,ivx)+brightx*xkf2*etan(i,j+1);
            end
            % Top Vx node
            if (i>1)
                ivxtop=ivx-3;
                L(ivx,ivxtop)=ykf*etas(i,j+1);
            else
                L(ivx,ivx)=L(ivx,ivx)+btopx*ykf*etas(i,j+1);
            end
            % Bottom Vx node
            if (i<ynum-1)
                ivxbottom=ivx+3;
                L(ivx,ivxbottom)=ykf*etas(i+1,j+1);
            else
                L(ivx,ivx)=L(ivx,ivx)+bbottomx*ykf*etas(i+1,j+1);
            end
            % Top Left Vy node
            if (i>1)
                ivytopleft=ivx-3+1;
                L(ivx,ivytopleft)=xykf*etas(i,j+1);
            else
                ivybottomleft=ivx+1;
                L(ivx,ivybottomleft)=btopy*xykf*etas(i,j+1);
            end
            % Top Right Vy node
            if (i>1)
                ivytopright=ivx-3+1+ynum3;
                L(ivx,ivytopright)=-xykf*etas(i,j+1);
            else
                ivybottomright=ivx+1+ynum3;
                L(ivx,ivybottomright)=-btopy*xykf*etas(i,j+1);
            end
            % Bottom Left Vy node
            if (i<ynum-1)
                ivybottomleft=ivx+1;
                if (i>1)
                    L(ivx,ivybottomleft)=-xykf*etas(i+1,j+1);
                else
                    L(ivx,ivybottomleft)=L(ivx,ivybottomleft)-xykf*etas(i+1,j+1);
                end
            else
                ivytopleft=ivx-3+1;
                L(ivx,ivytopleft)=L(ivx,ivytopleft)-bbottomy*xykf*etas(i+1,j+1);
            end
            % Bottom Right Vy node
            if (i<ynum-1)
                ivybottomright=ivx+1+ynum3;
                if (i>1)
                    L(ivx,ivybottomright)=xykf*etas(i+1,j+1);
                else
                    L(ivx,ivybottomright)=L(ivx,ivybottomright)+xykf*etas(i+1,j+1);
                end
            else
                ivytopright=ivx-3+1+ynum3;
                L(ivx,ivytopright)=L(ivx,ivytopright)+bbottomy*xykf*etas(i+1,j+1);
            end
            % Left P node
            iprleft=ivx+2;
            L(ivx,iprleft)=pscale/xstp;
            % Right P node
            iprright=ivx+2+ynum3;
            L(ivx,iprright)=-pscale/xstp;
            
        % Ghost Vx_parameter=0 used for numbering
        else
            L(ivx,ivx)=1;
            R(ivx,1)=0;
        end

            
            
        % y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
        if (i<ynum-1)
            % y-Stokes equation stensil
            %     +-------------------- -+-------vy(i-1,j)------+----------------------+    
            %     |                      |                      |                      |
            %     |                      |                      |                      |
            %     |                  vx(i,j-1)     P(i,j)    vx(i,j)                   |    
            %     |                      |        etan(i,j)     |                      |
            %     |                      |                      |                      |
            %     +-----vy(i,j-1)---etas(i+1,j)---vy(i,j)--etas(i+1,j+1)---vy(i,j+1)---+
            %     |                      |                      |                      |
            %     |                      |                      |                      |
            %     |                  vx(i+1,j-1)  P(i+1,j)   vx(i+1,j)                 |    
            %     |                      |      etan(i+1,j)     |                      |
            %     |                      |                      |                      |
            %     +----------------------+-------vy(i+1,j)------+----------------------+
            %
            % Right Part
            R(ivy,1)=RY(i+1,j+1);
            % Computing Current y-Stokes coefficients
            % Central Vy node
            L(ivy,ivy)=-ykf2*(etan(i+1,j)+etan(i,j))-xkf*(etas(i+1,j+1)+etas(i+1,j));
            % Top Vy node
            if(i>1)
                ivytop=ivy-3;
                L(ivy,ivytop)=ykf2*etan(i,j);
            else
                L(ivy,ivy)=L(ivy,ivy)+btopy*ykf2*etan(i,j);
            end
            % Bottom Vy node
            if(i<ynum-2)
                ivybottom=ivy+3;
                L(ivy,ivybottom)=ykf2*etan(i+1,j);
            else
                L(ivy,ivy)=L(ivy,ivy)+bbottomy*ykf2*etan(i+1,j);
            end
            % Left Vy node
            if(j>1)
                ivyleft=ivy-ynum3;
                L(ivy,ivyleft)=xkf*etas(i+1,j);
            else
                L(ivy,ivy)=L(ivy,ivy)+blefty*xkf*etas(i+1,j);
            end
            % Right Vy node
            if(j<xnum-1)
                ivyright=ivy+ynum3;
                L(ivy,ivyright)=xkf*etas(i+1,j+1);
            else
                L(ivy,ivy)=L(ivy,ivy)+brighty*xkf*etas(i+1,j+1);
            end
            % Top left Vx node
            if (j>1)
                ivxtopleft=ivy-1-ynum3;
                L(ivy,ivxtopleft)=xykf*etas(i+1,j);
            else
                ivxtopright=ivy-1;
                L(ivy,ivxtopright)=bleftx*xykf*etas(i+1,j);
            end
            % Bottom left Vx node
            if (j>1)
                ivxbottomleft=ivy-1+3-ynum3;
                L(ivy,ivxbottomleft)=-xykf*etas(i+1,j);
            else
                ivxbottomright=ivy-1+3;
                L(ivy,ivxbottomright)=-bleftx*xykf*etas(i+1,j);
            end
            % Top right Vx node
            if (j<xnum-1)
                ivxtopright=ivy-1;
                if(j>1)
                    L(ivy,ivxtopright)=-xykf*etas(i+1,j+1);
                else
                    L(ivy,ivxtopright)=L(ivy,ivxtopright)-xykf*etas(i+1,j+1);
                end
            else
                ivxtopleft=ivy-1-ynum3;
                L(ivy,ivxtopleft)=L(ivy,ivxtopleft)-brightx*xykf*etas(i+1,j+1);
            end
            % Bottom right Vx node
            if (j<xnum-1)
                ivxbottomright=ivy-1+3;
                if(j>1)
                    L(ivy,ivxbottomright)=xykf*etas(i+1,j+1);
                else
                    L(ivy,ivxbottomright)=L(ivy,ivxbottomright)+xykf*etas(i+1,j+1);
                end
            else
                ivxbottomleft=ivy-1+3-ynum3;
                L(ivy,ivxbottomleft)=L(ivy,ivxbottomleft)+brightx*xykf*etas(i+1,j+1);
           end
            % Top P node
            iprtop=ivy+1;
            L(ivy,iprtop)=pscale/ystp;
            % Bottom P node
            iprbottom=ivy+1+3;
            L(ivy,iprbottom)=-pscale/ystp;
            
        % Ghost Vy_parameter=0 used for numbering
        else
            L(ivy,ivy)=1;
            R(ivy,1)=0;
        end

        
        % Continuity equation dvx/dx+dvy/dy=RC
        if ( ((j>1 || i>1) && bpres==0) || (i>1 && i<ynum-1 && bpres==1) || (j>1 && j<xnum-1 && bpres==2) ) 
            % Continuity equation stensil
            %     +-----vy(i-1,j)--------+
            %     |                      |
            %     |                      |
            % vx(i,j-1)  pr(i,j)      vx(i,j)    
            %     |                      |
            %     |                      |
            %     +------vy(i,j)---------+
            %
            % Right Part
            R(ipr,1)=RC(i,j);
            % Computing Current Continuity coefficients
            % Left Vx node
            if (j>1)
                ivxleft=ipr-2-ynum3;
                L(ipr,ivxleft)=-pscale/xstp;
            end
            % Right Vx node
            if (j<xnum-1)
                ivxright=ipr-2;
                L(ipr,ivxright)=pscale/xstp;
            end
            % Top Vy node
            if (i>1)
                ivytop=ipr-1-3;
                L(ipr,ivytop)=-pscale/ystp;
            end
            % Bottom Vy node
            if (i<ynum-1)
                ivybottom=ipr-1;
                L(ipr,ivybottom)=pscale/ystp;
            end
            
        % Pressure definition for the boundary condition regions
        else
            % Pressure definition in one cell
            if (bpres==0)
                L(ipr,ipr)=2*pscale/(xstp+ystp);
                R(ipr,1)=2*prnorm/(xstp+ystp);
            end
            % Pressure definition at the top and bottom
            if (bpres==1)
                L(ipr,ipr)=2*pscale/(xstp+ystp);
                if (i==1)
                    R(ipr,1)=2*prnorm/(xstp+ystp);
                else
                    R(ipr,1)=0;
                end
            end
            % Pressure definition at the left and right
            if (bpres==2)
                L(ipr,ipr)=2*pscale/(xstp+ystp);
                if (j==1)
                    R(ipr,1)=2*prnorm/(xstp+ystp);
                else
                    R(ipr,1)=0;
                end
            end
        end
             
    end            
end


% Solve matrix
S=L\R;

% Reload solution
for i=1:1:ynum-1
    for j=1:1:xnum-1
        % Indexes for P,vx,vy
        ivx=((j-1)*(ynum-1)+(i-1))*3+1;
        ivy=ivx+1;
        ipr=ivx+2;
        % Reload Vx
        if (j<xnum-1)
            vx(i+1,j+1)=S(ivx);
        end
        % Reload Vy
        if (i<ynum-1)
            vy(i+1,j+1)=S(ivy);
        end
        % Reload P
        pr(i,j)=S(ipr)*pscale;
    end
end

% Apply vx boundary conditions
% Left Boundary
vx(:,1)=bleftx*vx(:,2);
% Right Boundary
vx(:,xnum)=brightx*vx(:,xnum-1);
% Upper Boundary
 vx(1,:)=btopx*vx(2,:);
% Lower Boundary
 vx(ynum+1,:)=bbottomx*vx(ynum,:);

% Apply vy boundary conditions
% Upper Boundary
 vy(1,:)=btopy*vy(2,:);
% Lower Boundary
 vy(ynum,:)=bbottomy*vy(ynum-1,:);
% Left Boundary
vy(:,1)=blefty*vy(:,2);
% Right Boundary
vy(:,xnum+1)=brighty*vy(:,xnum);

% Computing residuals
for i=1:1:ynum+1;
    for j=1:1:xnum+1;

        % x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
        if (j<xnum+1)
            % vx-Boundrary conditions 
            if (i==1 | i==ynum+1 | j==1 | j==xnum);
                resx(i,j)=0;
            % x-Stokes equation
            else
                % Computing Current x-Stokes residual
                % dSIGMAxx/dx-dP/dx
                resx(i,j)=RX(i,j)-(xkf2*(etan(i-1,j)*(vx(i,j+1)-vx(i,j))-etan(i-1,j-1)*(vx(i,j)-vx(i,j-1)))-(pr(i-1,j)-pr(i-1,j-1))/xstp);
                % dSIGMAxy/dy
                resx(i,j)=resx(i,j)-(etas(i,j)*((vx(i+1,j)-vx(i,j))*ykf+(vy(i,j+1)-vy(i,j))*xykf)-etas(i-1,j)*((vx(i,j)-vx(i-1,j))*ykf+(vy(i-1,j+1)-vy(i-1,j))*xykf));
            end
        end
            
        % y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
        if (i<ynum+1)
            % vy-Boundrary conditions 
            if (i==1 | i==ynum | j==1 | j==xnum+1);
                resy(i,j)=0;
            %y-Stokes equation
            else
                % Computing current residual
                % dSIGMAyy/dy-dP/dy
                resy(i,j)=RY(i,j)-(ykf2*(etan(i,j-1)*(vy(i+1,j)-vy(i,j))-etan(i-1,j-1)*(vy(i,j)-vy(i-1,j)))-(pr(i,j-1)-pr(i-1,j-1))/ystp);
                % dSIGMAxy/dx
                resy(i,j)=resy(i,j)-(etas(i,j)*((vy(i,j+1)-vy(i,j))*xkf+(vx(i+1,j)-vx(i,j))*xykf)-etas(i,j-1)*((vy(i,j)-vy(i,j-1))*xkf+(vx(i+1,j-1)-vx(i,j-1))*xykf));

            end
        end
            
        % Continuity equation dvx/dx+dvy/dy=RC
        if (i<ynum && j<xnum);
            % Computing current residual
            resc(i,j)=RC(i,j)-((vx(i+1,j+1)-vx(i+1,j))/xstp+(vy(i+1,j+1)-vy(i,j+1))/ystp);
        end
             
    end            
end


