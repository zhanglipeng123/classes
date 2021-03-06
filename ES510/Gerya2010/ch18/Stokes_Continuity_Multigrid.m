% Solution of 2D Stokes and continuity equations
% for a constant viscosity medium
% with multigrid based on V-cycle and external functions 
% Stokes_Continuity_smoother.m,
% Stokes_Continuity_restriction.m,
% Stokes_Continuity_prolongation.m;
% resolution between multigrid levels changes by factor of two

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

% Clearing all variables and arrays
clear all;
% Clearing figures
clf;

% % Model parameters
% Density distribution in the model corresponds to falling block test
% Block density
rhoblock=3100.0;
% Medium density
rhomedium=3000.0;

% Acceleration of Gravity, m/s^2
g=9.81;
% Model viscosity, Pa*s
etas=1e+20;
% Pressure in the upermost, lwftmost (first) cell
prfirst=0;
% Model size, m
xsize=100000.0;
ysize=100000.0;

% Multigrid parameters
% Number of resolution levels
levelnum=4;
% Resolution (number of nodal points) on the coarsest (last) grid
xnumcur=7;
ynumcur=7;

% Iteration Parameters
% Total number of iteration cycles
inum=1000;
% Relaxation coefficients for Gauss-Seidel iterations
% Stokes equations
krelaxs=0.9;
vkoef=1.0;
% Continuity equation
krelaxc=0.3;
pkoef=1.0;

% Number of smoothing iterations
iternum(1)=5;
iternum(2)=5*2;
iternum(3)=5*4;
iternum(4)=5*8;
iternum(5)=5*16;
iternum(6)=5*32;
iternum(7)=5*64;
iternum(8)=5*128;



% Defining resolutions and gridsteps on all levels
% Number 
for k=levelnum:-1:1
    % Numerical resolution for current level
    xnum(k)=xnumcur;
    ynum(k)=ynumcur;
    % Grid steps
    xstp(k)=xsize./(xnum(k)-1);
    ystp(k)=ysize./(ynum(k)-1);
    % Numerical resolution for the next finer level
    % Doubling number of CELLS
    xnumcur=(xnumcur-1)*2+1;
    ynumcur=(ynumcur-1)*2+1;
end

% FINEST (PRINCIPAL) GRID
% Defining density structure rho()
% Defining initial guesses for velocity vx() vy() and pressure pr()
% Computing right part of Stokes (RX, RY) and Continuity (RC) equation
% vx, vy, P
vx1=zeros(ynum(1)+1,xnum(1));
vy1=zeros(ynum(1),xnum(1)+1);
pr1=zeros(ynum(1)-1,xnum(1)-1);
% Right parts of equations
RX1=zeros(ynum(1)+1,xnum(1));
RY1=zeros(ynum(1),xnum(1)+1);
RC1=zeros(ynum(1)-1,xnum(1)-1);
% Grid points cycle
for i=1:1:ynum(1);
    for j=1:1:xnum(1);
        % Relative distances for (i,j) basic node inside the grid
        dx=(j-1)*xstp(1)/xsize;
        dy=(i-1)*ystp(1)/ysize;
        % Density Structure (falling block)
        rho1(i,j)=rhomedium;
        if(dx>=0.4 && dx<=0.6 && dy>=0.4 && dy<=0.6)
            rho1(i,j)=rhoblock;
        end
        % Pressure update in internal nodes (lithostatic pressure)
        if(i<ynum(1) && j<xnum(1))
            if(i==1)
                pr1(i,j)=prfirst;
            else
                pr1(i,j)=pr1(i-1,j)+ystp(1)*g*(rho1(i,j)+rho1(i,j+1))/2;
            end
        end
        % Right part of x-Stokes Equation
        if(j>1 && i>1 && j<xnum(1))
            RX1(i,j)=0;
        end
        % Right part of y-Stokes Equation
        if(j>1 && i>1 && i<ynum(1))
            RY1(i,j)=-g*(rho1(i,j)+rho1(i,j-1))/2;
        end
    end
end


% Main Multigrid cycle
% We put 8 levels (for the case if we want to increase resolution)
% We use different arrays for different levels 
% since grid resolutions is different
% Multigrid schedule = V-cycle 
for niter=1:1:inum;
    
    
    % Smoothing+restriction cycle
    for k=1:1:levelnum;
        % Action depends on the level of resolution
        switch k
            
            % Leve 1 (principal grid)
            case 1
            % Smoothing operation: solving of Stokes and Continuity equations on nodes
            % and computing residuals
            % by calling function Stokes_smoother()
            [vx1,resx1,vy1,resy1,pr1,resc1]=Stokes_Continuity_smoother(prfirst,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX1,vx1,RY1,vy1,RC1,pr1);
            % Restriction operation: 
            % Interpolating residuals to coarcer level (k+1) 
            % to produce right parts for this coarcer level
            if (levelnum>k)
                [RX2,RY2,RC2]=Stokes_Continuity_restriction(k,xnum,ynum,xstp,ystp,resx1,resy1,resc1);
            end
        
            % Level 2
            case 2
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx2=zeros(ynum(k)+1,xnum(k));
            vy2=zeros(ynum(k),xnum(k)+1);
            pr2=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx2,resx2,vy2,resy2,pr2,resc2]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX2,vx2,RY2,vy2,RC2,pr2);
            % Restriction operation: 
            if (levelnum>k)
                [RX3,RY3,RC3]=Stokes_Continuity_restriction(k,xnum,ynum,xstp,ystp,resx2,resy2,resc2);
            end
        
            % Level 3
            case 3
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx3=zeros(ynum(k)+1,xnum(k));
            vy3=zeros(ynum(k),xnum(k)+1);
            pr3=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx3,resx3,vy3,resy3,pr3,resc3]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX3,vx3,RY3,vy3,RC3,pr3);
            % Restriction operation: 
            if (levelnum>k)
                [RX4,RY4,RC4]=Stokes_Continuity_restriction(k,xnum,ynum,xstp,ystp,resx3,resy3,resc3);
            end
            
            % Level 4
            case 4
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx4=zeros(ynum(k)+1,xnum(k));
            vy4=zeros(ynum(k),xnum(k)+1);
            pr4=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx4,resx4,vy4,resy4,pr4,resc4]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX4,vx4,RY4,vy4,RC4,pr4);
            % Restriction operation: 
            if (levelnum>k)
                [RX5,RY5,RC5]=Stokes_Continuity_restriction(k,xnum,ynum,xstp,ystp,resx4,resy4,resc4);
            end
            
             % Level 5
            case 5
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx5=zeros(ynum(k)+1,xnum(k));
            vy5=zeros(ynum(k),xnum(k)+1);
            pr5=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx5,resx5,vy5,resy5,pr5,resc5]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX5,vx5,RY5,vy5,RC5,pr5);
            % Restriction operation: 
            if (levelnum>k)
                [RX6,RY6,RC6]=Stokes_Continuity_restriction(k,xnum,ynum,xstp,ystp,resx5,resy5,resc5);
            end
             
            % Level 6
            case 6
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx6=zeros(ynum(k)+1,xnum(k));
            vy6=zeros(ynum(k),xnum(k)+1);
            pr6=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx6,resx6,vy6,resy6,pr6,resc6]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX6,vx6,RY6,vy6,RC6,pr6);
            % Restriction operation: 
            if (levelnum>k)
                [RX7,RY7,RC7]=Stokes_Continuity_restriction(k,xnum,ynum,xstp,ystp,resx6,resy6,resc6);
            end
             
            % Level 7
            case 7
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx7=zeros(ynum(k)+1,xnum(k));
            vy7=zeros(ynum(k),xnum(k)+1);
            pr7=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx7,resx7,vy7,resy7,pr7,resc7]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX7,vx7,RY7,vy7,RC7,pr7);
            % Restriction operation: 
            if (levelnum>k)
                [RX8,RY8,RC8]=Stokes_Continuity_restriction(k,xnum,ynum,xstp,ystp,resx7,resy7,resc7);
            end
             
            % Level 8
            case 8
            %Making initial approximation for vx,vy, P corrections (zeros)
            vx8=zeros(ynum(k)+1,xnum(k));
            vy8=zeros(ynum(k),xnum(k)+1);
            pr8=zeros(ynum(k)-1,xnum(k)-1);
            % Smoothing operation: 
            [vx8,resx8,vy8,resy8,pr8,resc8]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX8,vx8,RY8,vy8,RC8,pr8);
            % No Restriction operation on this last level 

        end
    end
   
    % Smoothing+prolongation cycle
    for k=levelnum:-1:1;
        % Action depends on the level of resolution
        switch k

            % Leve 1 (principal grid)
            case 1
            % Smoothing operation: solving of Stokes and Continuity equations on nodes
            % and computing residuals
            % by calling function Stokes_smoother()
            [vx1,resx1,vy1,resy1,pr1,resc1]=Stokes_Continuity_smoother(prfirst,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX1,vx1,RY1,vy1,RC1,pr1);
            % No prolongation operation on finest level
        
            % Level 2
            case 2
            % Smoothing operation: 
            [vx2,resx2,vy2,resy2,pr2,resc2]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX2,vx2,RY2,vy2,RC2,pr2);
            % Prolongation operation: 
            % Interpolating dvx,dvy,dpr corrections to finer level (k+1) 
            % and computing solution update for this finer level
            [dvx1,dvy1,dpr1]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx2,vy2,pr2);
            vx1=vx1+dvx1*vkoef;
            vy1=vy1+dvy1*vkoef;
            pr1=pr1+dpr1*pkoef;
        
            % Level 3
            case 3
            % Smoothing operation: 
            [vx3,resx3,vy3,resy3,pr3,resc3]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX3,vx3,RY3,vy3,RC3,pr3);
            % Prolongation operation: 
            [dvx2,dvy2,dpr2]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx3,vy3,pr3);
            vx2=vx2+dvx2*vkoef;
            vy2=vy2+dvy2*vkoef;
            pr2=pr2+dpr2*pkoef;
            
            % Level 4
            case 4
            % Smoothing operation: 
            [vx4,resx4,vy4,resy4,pr4,resc4]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX4,vx4,RY4,vy4,RC4,pr4);
            % Prolongation operation: 
            [dvx3,dvy3,dpr3]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx4,vy4,pr4);
            vx3=vx3+dvx3*vkoef;
            vy3=vy3+dvy3*vkoef;
            pr3=pr3+dpr3*pkoef;
            
             % Level 5
            case 5
            % Smoothing operation: 
            [vx5,resx5,vy5,resy5,pr5,resc5]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX5,vx5,RY5,vy5,RC5,pr5);
            % Prolongation operation: 
            [dvx4,dvy4,dpr4]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx5,vy5,pr5);
            vx4=vx4+dvx4*vkoef;
            vy4=vy4+dvy4*vkoef;
            pr4=pr4+dpr4*pkoef;
            
            
            % Level 6
            case 6
            % Smoothing operation: 
            [vx6,resx6,vy6,resy6,pr6,resc6]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX6,vx6,RY6,vy6,RC6,pr6);
            % Prolongation operation: 
            [dvx5,dvy5,dpr5]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx6,vy6,pr6);
            vx5=vx5+dvx5*vkoef;
            vy5=vy5+dvy5*vkoef;
            pr5=pr5+dpr5*pkoef;
            
            % Level 7
            case 7
            % Smoothing operation: 
            [vx7,resx7,vy7,resy7,pr7,resc7]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX7,vx7,RY7,vy7,RC7,pr7);
            % Prolongation operation: 
            [dvx6,dvy6,dpr6]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx7,vy7,pr7);
            vx6=vx6+dvx6*vkoef;
            vy6=vy6+dvy6*vkoef;
            pr6=pr6+dpr6*pkoef;
            
            % Level 8
            case 8
            % Smoothing operation: 
            [vx8,resx8,vy8,resy8,pr8,resc8]=Stokes_Continuity_smoother(0,etas,iternum(k),krelaxs,krelaxc,xnum(k),ynum(k),xstp(k),ystp(k),RX8,vx8,RY8,vy8,RC8,pr8);
            % Prolongation operation: 
            [dvx7,dvy7,dpr7]=Stokes_Continuity_prolongation(k,xnum,ynum,xstp,ystp,vx8,vy8,pr8);
            vx7=vx7+dvx7*vkoef;
            vy7=vy7+dvy7*vkoef;
            pr7=pr7+dpr7*pkoef;
        end
    end 
 
    % Defining scale for Stokes residuals from y-Stokes equation
    % dSIGMAij/dj-dP/di=-RHO*gi=0  => Stokes scale=abs(RHO*gi)
    stokesscale= rhoblock*g;
    % Defining scale for Continuity residuals from y-Stokes equation
    % dvx/dx+dvy/dy=0 can be transformed to 2ETA(dvx/dx+dvy/dy)/dx=0 
    % which is similar to dSIGMAij/dj and has scale given above
    % therefore continuity scale = scale=abs(RHO*gi/ETA*dx)
    continscale= rhoblock*g/etas*xstp(1);


    % Plotting Residuals for x-Stokes as surface
    figure(1);
    subplot(2,2,1)
    resx0=resx1/stokesscale;
    surf(resx0);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual x-Stokes')
%     title(['x-Stokes residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Plotting Residuals for y-Stokes as surface
    subplot(2,2,2)
    resy0=resy1/stokesscale;
    surf(resy0);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual Y-stokes')
%     title(['y-Stokes residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Plotting Residuals for Continuity as surface
    subplot(2,2,3)
    resc0=resc1/continscale;
    surf(resc0);
    shading interp;
    light;
    lighting phong;
    axis tight;
    zlabel('residual continuity')
%     title(['Continuity residual, Gauss-Seidel cycle = ',num2str(niter)])
%     colorbar;
    
    % Computing mean square residuals
    resx00(niter)=0;
    resy00(niter)=0;
    resc00(niter)=0;
    for i=1:1:ynum(1)
        for j=1:1:xnum(1)
            % x-Stokes
            if (i>1 && j>1 && j<xnum(1))
                resx00(niter)=resx00(niter)+resx0(i,j)^2;
            end
            % y-Stokes
            if (i>1 && j>1 && i<ynum(1))
                resy00(niter)=resy00(niter)+resy0(i,j)^2;
            end
            % Continuity
            if (j<xnum(1) && i<ynum(1))
                resc00(niter)=resc00(niter)+resc0(i,j)^2;
            end
        end
    end
    resx00(niter)=log10((resx00(niter)/((ynum(1)-1)*(xnum(1)-2)))^0.5);
    resy00(niter)=log10((resy00(niter)/((ynum(1)-2)*(xnum(1)-1)))^0.5);
    resc00(niter)=log10((resc00(niter)/((ynum(1)-1)*(xnum(1)-1)))^0.5);
   

    % Plotting Mean residuals
    subplot(2,2,4)
    plot(resc00, 'k');
    hold on
    plot(resx00, 'b');
    hold on
    plot(resy00, 'r');
    xlabel('V-cycles')
    ylabel('log(Residuals)')

    zlabel(['RMS errors, cycle ',num2str(niter)])
%     title(['RMS errors, Multigrid V-cycle = ',num2str(niter)])
%     colorbar;

    pause(1)

end

