% Solving Navier-Stokes and continuity eq.
% in primitive variable formulation
% with variable viscosity 
% for visco-elastic material
% using  FD with staggered grid
% and marker in cell techique

% Clearing memory and figures
clear all; clf

% Define Numerical model
xsize=10; % Horizontal model size, m
ysize=10; % Vertical model size, m
Nx=51; % Horizontal grid resolution
Ny=51; % Vertical grid resolution
Nx1=Nx+1;
Ny1=Ny+1;
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step, m

% Define Gravity
gx=0; % Horizontal gravity acceleration, m/s^2
gy=10; % Vertical gravity acceleration, m/s^2

% Coordinates of different nodal points
% Basic nodes
x=0:dx:xsize; % Horizontal coordinates of basic grid points, m
y=0:dy:ysize; % Vertical coordinates of basic grid points, m
% Vx-Nodes
xvx=0:dx:xsize+dy; % Horizontal coordinates of vx grid points, m
yvx=-dy/2:dy:ysize+dy/2; % Vertical coordinates of vx grid points, m
% Vy-nodes
xvy=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of vy grid points, m
yvy=0:dy:ysize+dy; % Vertical coordinates of vy grid points, m
% P-Nodes
xp=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of P grid points, m
yp=-dy/2:dy:ysize+dy/2; % Vertical coordinates of P grid points, m

% Nodal arrays
% Basic nodes
ETA=zeros(Ny,Nx); % Viscosity, Pa*s
GGG=zeros(Ny,Nx); % Shear modulus, Pa
EXY=zeros(Ny,Nx); % EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % SIGMAxy, 1/s
SXY0=zeros(Ny,Nx); % SIGMA0xy, 1/s
wyx=zeros(Ny,Nx); % Rotation rate, 1/s
% Vx-Nodes
RHOX=zeros(Ny1,Nx1); % Density, kg/m^3
vx=zeros(Ny1,Nx1); % vx-velocity m/s
vx0=zeros(Ny1,Nx1); % Old vx-velocity m/s
% Vy-Nodes
RHOY=zeros(Ny1,Nx1); % Density, kg/m^3
vy=zeros(Ny1,Nx1); % vy-velocity m/s
vy0=zeros(Ny1,Nx1); % Old vy-velocity m/s
% P-nodes
RHO=zeros(Ny1,Nx1); % Density, kg/m^3
BETTA=zeros(Ny1,Nx1); % Compressibility, 1/Pa
ETAP=zeros(Ny1,Nx1); % Viscosity, Pa*s
GGGP=zeros(Ny1,Nx1); % Shear modulus, Pa
EXX=zeros(Ny,Nx); % EPSILONxx, 1/s
SXX=zeros(Ny,Nx); % SIGMA'xx, 1/s
SXX0=zeros(Ny,Nx); % SIGMA0'xx, 1/s
vxp=zeros(Ny1,Nx1); % Vx in pressure nodes, m/s
vyp=zeros(Ny1,Nx1); % Vy in pressure nodes, m/s
pr=zeros(Ny1,Nx1); % Pressure, Pa
pr0=zeros(Ny1,Nx1); % Old pressure, Pa

% Define markers
Nxm=(Nx-1)*4; % Marker grid resolution in horizontal direction
Nym=(Ny-1)*4; % Marker grid resolution in vertical direction
dxm=xsize/Nxm; % Marker grid step in horizontal direction,m
dym=ysize/Nym; % Marker grid step in vertical direction,m
marknum=Nxm*Nym+1; % Number of markers (including traced point)
xm=zeros(1,marknum); % Horizontal coordinates, m
ym=zeros(1,marknum); % Vertical coordinates, m
tm=zeros(1,marknum); % Material type
sxxm=zeros(1,marknum); % SIGMA'xx, Pa
sxym=zeros(1,marknum); % SIGMAxy, Pa
pr0m=zeros(1,marknum); % Marker pressure, Pa
vx0m=zeros(1,marknum); % Marker vx, m/s
vy0m=zeros(1,marknum); % Marker vy, m/s

% Define properties of materials: 
%        Ball   Air
rhom   = [3000   1     ]; % Density, kg/m^3
etam   = [1e+22  1e-3  ]; % Viscosity, Pa s
bettam = [1e-8   1e-8 ]; % Compressibility, 1/Pa
gggm   = [1e+8   1e+8  ]; % Shear Modulus, Pa

% Define marker coordinates, temperature and material type
m=1; % Marker counter
for jm=1:1:Nxm
    for im=1:1:Nym
        % Define marker coordinates
        xm(m)=dxm/2+(jm-1)*dxm+(rand-0.5)*dxm;
        ym(m)=dym/2+(im-1)*dym+(rand-0.5)*dym;
        % Marker properties
        % Air
        tm(m)=2; % Material type
        % Ball
        if(((xm(m)-xsize*0.5)^2+(ym(m)-ysize*0.25)^2)^0.5<ysize*0.2)
            tm(m)=1; % Material type
        end
        % Update marker counter
        m=m+1;
    end
end
% Tracing point = Ball center
m=marknum;
xm(m)=xsize*0.5;
ym(m)=ysize*0.25;
tm(m)=1;


% Define global matrixes 
% Mechanical solution: L(), R()
N=Nx1*Ny1*3; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts

% Mechanical boundary conditions: free slip=-1; No Slip=1
bcleft=-1;
bcright=-1;
bctop=-1;
bcbottom=-1;
% Pressure BC
PCONF=1e+5;

% Timestepping
dtelastic=1; % Maximal computational timestep, s
dtmin=1e-6; % Minimal computational timestep, s
dt=0.1; % Current computational timestep, s
dtkoefv=1.1; % Koefficient to decrese dt for P,T,vx,vy,SIGMA limits 
dtkoefup=1.1; % Koefficient to increase dt
dxymax=0.1; % Max marker movement per time step, grid steps
vpratio=1/3; % Weight of averaged velocity for moving markers
DSmax=3e+5; % Max stress change per time step, Pa
DPmax=3e+5; % Max stress change per time step, Pa
dsubgrids=0; % Subgrid stress diffusion parameter
dsubgridv=0; % Subgrid velocity diffusion parameter
timesum=0; % Time sum, s
etamin=1e-3; % Lower viscosity cut-off, Pa s
etamax=1e+23; % Upper viscosity cut-off, Pa s
visstep=100; % Periodicity of visualization
nsteps=100000; % number of timesteps
imax=1000; % maximal numbers of timestep adjustment iterations
timestep=1;

for timestep=timestep:1:nsteps

% Save old stresses
sxxm00=sxxm; 
sxym00=sxym;    
    
% Interpolate properties from markers to nodes
% Basic nodes
ETASUM=zeros(Ny,Nx);
GGGSUM=zeros(Ny,Nx);
SXYSUM=zeros(Ny,Nx);
WTSUM=zeros(Ny,Nx);
% Vx-nodes
RHOXSUM=zeros(Ny1,Nx1);
VXSUM=zeros(Ny1,Nx1);
WTXSUM=zeros(Ny1,Nx1);
% Vy-nodes
RHOYSUM=zeros(Ny1,Nx1);
VYSUM=zeros(Ny1,Nx1);
WTYSUM=zeros(Ny1,Nx1);
% P-Nodes
ETAPSUM=zeros(Ny1,Nx1);
GGGPSUM=zeros(Ny1,Nx1);
SXXSUM=zeros(Ny1,Nx1);
RHOSUM=zeros(Ny1,Nx1);
BETTASUM=zeros(Ny1,Nx1);
PRSUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);

for m=1:1:marknum
    
    % Interpolation to basic nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    ETASUM(i,j)=ETASUM(i,j)+etam(tm(m))*wtmij;
    GGGSUM(i,j)=GGGSUM(i,j)+1/gggm(tm(m))*wtmij;
    SXYSUM(i,j)=SXYSUM(i,j)+sxym(m)*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    ETASUM(i+1,j)=ETASUM(i+1,j)+etam(tm(m))*wtmi1j;
    GGGSUM(i+1,j)=GGGSUM(i+1,j)+1/gggm(tm(m))*wtmi1j;
    SXYSUM(i+1,j)=SXYSUM(i+1,j)+sxym(m)*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    ETASUM(i,j+1)=ETASUM(i,j+1)+etam(tm(m))*wtmij1;
    GGGSUM(i,j+1)=GGGSUM(i,j+1)+1/gggm(tm(m))*wtmij1;
    SXYSUM(i,j+1)=SXYSUM(i,j+1)+sxym(m)*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    ETASUM(i+1,j+1)=ETASUM(i+1,j+1)+etam(tm(m))*wtmi1j1;
    GGGSUM(i+1,j+1)=GGGSUM(i+1,j+1)+1/gggm(tm(m))*wtmi1j1;
    SXYSUM(i+1,j+1)=SXYSUM(i+1,j+1)+sxym(m)*wtmi1j1;
    WTSUM(i+1,j+1)=WTSUM(i+1,j+1)+wtmi1j1;    
    
    
    % Interpolation to vx-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvx(1))/dx)+1;
    i=fix((ym(m)-yvx(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xvx(j);
    dymi=ym(m)-yvx(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    RHOXSUM(i,j)=RHOXSUM(i,j)+rhom(tm(m))*wtmij;
    VXSUM(i,j)=VXSUM(i,j)+vx0m(m)*rhom(tm(m))*wtmij;
    WTXSUM(i,j)=WTXSUM(i,j)+wtmij;
    % i+1,j Node
    RHOXSUM(i+1,j)=RHOXSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    VXSUM(i+1,j)=VXSUM(i+1,j)+vx0m(m)*rhom(tm(m))*wtmi1j;
    WTXSUM(i+1,j)=WTXSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOXSUM(i,j+1)=RHOXSUM(i,j+1)+rhom(tm(m))*wtmij1;
    VXSUM(i,j+1)=VXSUM(i,j+1)+vx0m(m)*rhom(tm(m))*wtmij1;
    WTXSUM(i,j+1)=WTXSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOXSUM(i+1,j+1)=RHOXSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    VXSUM(i+1,j+1)=VXSUM(i+1,j+1)+vx0m(m)*rhom(tm(m))*wtmi1j1;
    WTXSUM(i+1,j+1)=WTXSUM(i+1,j+1)+wtmi1j1;

    % Interpolation to vy-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvy(1))/dx)+1;
    i=fix((ym(m)-yvy(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xvy(j);
    dymi=ym(m)-yvy(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    RHOYSUM(i,j)=RHOYSUM(i,j)+rhom(tm(m))*wtmij;
    VYSUM(i,j)=VYSUM(i,j)+vy0m(m)*rhom(tm(m))*wtmij;
    WTYSUM(i,j)=WTYSUM(i,j)+wtmij;
    % i+1,j Node
    RHOYSUM(i+1,j)=RHOYSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    VYSUM(i+1,j)=VYSUM(i+1,j)+vy0m(m)*rhom(tm(m))*wtmi1j;
    WTYSUM(i+1,j)=WTYSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOYSUM(i,j+1)=RHOYSUM(i,j+1)+rhom(tm(m))*wtmij1;
    VYSUM(i,j+1)=VYSUM(i,j+1)+vy0m(m)*rhom(tm(m))*wtmij1;
    WTYSUM(i,j+1)=WTYSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOYSUM(i+1,j+1)=RHOYSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    VYSUM(i+1,j+1)=VYSUM(i+1,j+1)+vy0m(m)*rhom(tm(m))*wtmi1j1;
    WTYSUM(i+1,j+1)=WTYSUM(i+1,j+1)+wtmi1j1;

    % Interpolation to P-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    ETAPSUM(i,j)=ETAPSUM(i,j)+etam(tm(m))*wtmij;
    GGGPSUM(i,j)=GGGPSUM(i,j)+1/gggm(tm(m))*wtmij;
    SXXSUM(i,j)=SXXSUM(i,j)+sxxm(m)*wtmij;
    RHOSUM(i,j)=RHOSUM(i,j)+rhom(tm(m))*wtmij;
    BETTASUM(i,j)=BETTASUM(i,j)+bettam(tm(m))*wtmij;
    PRSUM(i,j)=PRSUM(i,j)+pr0m(m)*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    ETAPSUM(i+1,j)=ETAPSUM(i+1,j)+etam(tm(m))*wtmi1j;
    GGGPSUM(i+1,j)=GGGPSUM(i+1,j)+1/gggm(tm(m))*wtmi1j;
    SXXSUM(i+1,j)=SXXSUM(i+1,j)+sxxm(m)*wtmi1j;
    RHOSUM(i+1,j)=RHOSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    BETTASUM(i+1,j)=BETTASUM(i+1,j)+bettam(tm(m))*wtmi1j;
    PRSUM(i+1,j)=PRSUM(i+1,j)+pr0m(m)*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    ETAPSUM(i,j+1)=ETAPSUM(i,j+1)+etam(tm(m))*wtmij1;
    GGGPSUM(i,j+1)=GGGPSUM(i,j+1)+1/gggm(tm(m))*wtmij1;
    SXXSUM(i,j+1)=SXXSUM(i,j+1)+sxxm(m)*wtmij1;
    RHOSUM(i,j+1)=RHOSUM(i,j+1)+rhom(tm(m))*wtmij1;
    BETTASUM(i,j+1)=BETTASUM(i,j+1)+bettam(tm(m))*wtmij1;
    PRSUM(i,j+1)=PRSUM(i,j+1)+pr0m(m)*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    ETAPSUM(i+1,j+1)=ETAPSUM(i+1,j+1)+etam(tm(m))*wtmi1j1;
    GGGPSUM(i+1,j+1)=GGGPSUM(i+1,j+1)+1/gggm(tm(m))*wtmi1j1;
    SXXSUM(i+1,j+1)=SXXSUM(i+1,j+1)+sxxm(m)*wtmi1j1;
    RHOSUM(i+1,j+1)=RHOSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    BETTASUM(i+1,j+1)=BETTASUM(i+1,j+1)+bettam(tm(m))*wtmi1j1;
    PRSUM(i+1,j+1)=PRSUM(i+1,j+1)+pr0m(m)*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;
end
% Compute physical properties
% Basic nodes
YNY=zeros(Ny,Nx);
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            ETA(i,j)=ETASUM(i,j)/WTSUM(i,j);
            GGG(i,j)=1/(GGGSUM(i,j)/WTSUM(i,j));
            SXY0(i,j)=SXYSUM(i,j)/WTSUM(i,j);
        end
    end
end
% Vx-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTXSUM(i,j)>0)
            RHOX(i,j)=RHOXSUM(i,j)/WTXSUM(i,j);
            vx0(i,j)=VXSUM(i,j)/RHOXSUM(i,j);
        end
    end
end
% Vy-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTYSUM(i,j)>0)
            RHOY(i,j)=RHOYSUM(i,j)/WTYSUM(i,j);
            vy0(i,j)=VYSUM(i,j)/RHOYSUM(i,j);
        end
    end
end
% P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM(i,j)>0)
            ETAP(i,j)=ETAPSUM(i,j)/WTPSUM(i,j);
            GGGP(i,j)=1/(GGGPSUM(i,j)/WTPSUM(i,j));
            SXX0(i,j)=SXXSUM(i,j)/WTPSUM(i,j);
            RHO(i,j)=RHOSUM(i,j)/WTPSUM(i,j);
            BETTA(i,j)=BETTASUM(i,j)/WTPSUM(i,j);
            pr0(i,j)=PRSUM(i,j)/WTPSUM(i,j);
        end
    end
end
% Applying vx-velocity boundary conditions for interpolated temperature
% Upper boundary 
vx0(1,2:Nx-1)=-bctop*vx0(2,2:Nx-1); % free slip
% Lower boundary 
vx0(Ny1,2:Nx-1)=-bcbottom*vx0(Ny,2:Nx-1); % free slip
% Left boundary
vx0(:,1)=0; % no slip/constant velocity
% Right boundary
vx0(:,Nx)=0; % no slip/constant velocity
% Applying vy-velocity boundary conditions for interpolated temperature
% Left boundary
vy0(2:Ny-1,1)=-bcleft*vy0(2:Ny-1,2); % free slip
% Right boundary
vy0(2:Ny-1,Nx1)=-bcright*vy0(2:Ny-1,Nx); % free slip
% Upper boundary 
vy0(1,:)=0; % no slip/constant velocity
% Lower boundary 
vy0(Ny,:)=0; % no slip/constant velocity
% Applying Pressure boundary conditions for interpolated pressure
dpr0=pr0(2,2)-PCONF;
pr0=pr0-dpr0;
pr0m=pr0m-dpr0;

% Try to increase computational Timestep
dt=min(dt*dtkoefup,dtelastic);


% Save initial viscosity
ETA00=ETA;
ETAP00=ETAP;

% Start thermomechanical iterations on Nodes
for inum=1:1:imax
% Set lower limit for viscosity in Navier-Stokes equations
etaminNS=max(etamin, dt*min(min(GGG))*1e-12);
% % Limit viscosity in the basic nodes
% ETA=ETA00;
% for i=1:1:Ny
%     for j=1:1:Nx
%         if(ETA(i,j)<etaminNS)
%         ETA(i,j)=etaminNS;
%         end
%     end
% end
% % Limit viscosity in the pressure nodes
% ETAP=ETAP00;
% for i=2:1:Ny
%     for j=2:1:Nx
%         if(ETAP(i,j)<etaminNS)
%         ETAP(i,j)=etaminNS;
%         end
%     end
% end
% Introducing scaled pressure
pscale=etaminNS/dx;

% Composing global matrixes L(), R() for Navier-Stokes and continuity equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*3+1; % Vx
        kvy=kvx+1; % Vy
        kpm=kvx+2; % P
        
        % Vx equation External points
        if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
            % Boundary Condition 
            % Ghost unknowns 1*Vx=0
            if(j==Nx1)
                L(kvx,kvx)=1; % Left part
                R(kvx)=0; % Right part
            end
            % Left Boundary
            if(j==1)
                L(kvx,kvx)=1; % Left part
                R(kvx)=0; % Right part
            end
            % Right Boundary
            if(j==Nx)
                L(kvx,kvx)=1; % Left part
                R(kvx)=0; % Right part
            end
            % Top boundary
            if(i==1 && j>1 && j<Nx)
                L(kvx,kvx)=1; % Left part
                L(kvx,kvx+3)=bctop; % Left part
                R(kvx)=0; % Right part
            end
            % Top boundary
            if(i==Ny1 && j>1 && j<Nx)
                L(kvx,kvx)=1; % Left part
                L(kvx,kvx-3)=bcbottom; % Left part
                R(kvx)=0; % Right part
            end
        else
        % Internal points: x-Navier-Stokes eq.
        %            Vx2
        %             |
        %        Vy1  |  Vy3
        %             |
        %     Vx1-P1-Vx3-P2-Vx5
        %             |
        %        Vy2  |  Vy4
        %             |
        %            Vx4
        %
        % Computational viscosity
        ETA1=ETA(i-1,j)*GGG(i-1,j)*dt/(GGG(i-1,j)*dt+ETA(i-1,j));
        ETA2=ETA(i,j)*GGG(i,j)*dt/(GGG(i,j)*dt+ETA(i,j));
        ETAP1=ETAP(i,j)*GGGP(i,j)*dt/(GGGP(i,j)*dt+ETAP(i,j));
        ETAP2=ETAP(i,j+1)*GGGP(i,j+1)*dt/(GGGP(i,j+1)*dt+ETAP(i,j+1));
        % Old stresses
        SXY1=SXY0(i-1,j)*ETA(i-1,j)/(GGG(i-1,j)*dt+ETA(i-1,j));
        SXY2=SXY0(i,j)*ETA(i,j)/(GGG(i,j)*dt+ETA(i,j));
        SXX1=SXX0(i,j)*ETAP(i,j)/(GGGP(i,j)*dt+ETAP(i,j));
        SXX2=SXX0(i,j+1)*ETAP(i,j+1)/(GGGP(i,j+1)*dt+ETAP(i,j+1));
        % Density gradients
        dRHOdx=(RHOX(i,j+1)-RHOX(i,j-1))/2/dx;
        dRHOdy=(RHOX(i+1,j)-RHOX(i-1,j))/2/dy;
        % Left part
        L(kvx,kvx-Ny1*3)=ETAP1/dx^2; % Vx1
        L(kvx,kvx-3)=ETA1/dy^2; % Vx2
        L(kvx,kvx)=-(ETAP1+ETAP2)/dx^2-...
                      (ETA1+ETA2)/dy^2-...
                      dRHOdx*gx*dt-RHOX(i,j)/dt; % Vx3
        L(kvx,kvx+3)=ETA2/dy^2; % Vx4
        L(kvx,kvx+Ny1*3)=ETAP2/dx^2; % Vx5
        L(kvx,kvy)=ETAP1/dx/dy-ETA2/dx/dy-dRHOdy*gx*dt/4;  % Vy2
        L(kvx,kvy+Ny1*3)=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdy*gx*dt/4;  % Vy4
        L(kvx,kvy-3)=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdy*gx*dt/4;  % Vy1
        L(kvx,kvy+Ny1*3-3)=ETAP2/dx/dy-ETA1/dx/dy-dRHOdy*gx*dt/4;  % Vy3
        L(kvx,kpm)=pscale/dx; % P1
        L(kvx,kpm+Ny1*3)=-pscale/dx; % P2
        % Right part
        R(kvx)=-RHOX(i,j)*(gx+vx0(i,j)/dt)-(SXY2-SXY1)/dy-(SXX2-SXX1)/dx;
        end
        
        % Vy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            % Boundary Condition
            % Ghost unknowns 1*Vx=0
            if(i==Ny1)
                L(kvy,kvy)=1; % Left part
                R(kvy)=0; % Right part
            end
            % Top boundary
            if(i==1)
                L(kvy,kvy)=1; % Left part
                R(kvy)=0; % Right part
            end
            % Bottom boundary
            if(i==Ny)
                L(kvy,kvy)=1; % Left part
                R(kvy)=0; % Right part
            end
            % Left boundary
            if(j==1 && i>1 && i<Ny)
                L(kvy,kvy)=1; % Left part
                L(kvy,kvy+3*Ny1)=bcleft; % Left part
                R(kvy)=0; % Right part
            end
            % Right boundary
            if(j==Nx1 && i>1 && i<Ny)
                L(kvy,kvy)=1; % Left part
                L(kvy,kvy-3*Ny1)=bcright; % Left part
                R(kvy)=0; % Right part
            end
        else
        % Internal points: y-Navier-Stokes eq.
        %            Vy2
        %             |
        %         Vx1 P1 Vx3
        %             |
        %     Vy1----Vy3----Vy5
        %             |
        %         Vx2 P2 Vx4
        %             |
        %            Vy4
        %
        % Computational viscosity
        ETA1=ETA(i,j-1)*GGG(i,j-1)*dt/(GGG(i,j-1)*dt+ETA(i,j-1));
        ETA2=ETA(i,j)*GGG(i,j)*dt/(GGG(i,j)*dt+ETA(i,j));
        ETAP1=ETAP(i,j)*GGGP(i,j)*dt/(GGGP(i,j)*dt+ETAP(i,j));
        ETAP2=ETAP(i+1,j)*GGGP(i+1,j)*dt/(GGGP(i+1,j)*dt+ETAP(i+1,j));
        % Old stresses
        SXY1=SXY0(i,j-1)*ETA(i,j-1)/(GGG(i,j-1)*dt+ETA(i,j-1));
        SXY2=SXY0(i,j)*ETA(i,j)/(GGG(i,j)*dt+ETA(i,j));
        SYY1=-SXX0(i,j)*ETAP(i,j)/(GGGP(i,j)*dt+ETAP(i,j));
        SYY2=-SXX0(i+1,j)*ETAP(i+1,j)/(GGGP(i+1,j)*dt+ETAP(i+1,j));
        % Density gradients
        dRHOdx=(RHOY(i,j+1)-RHOY(i,j-1))/2/dx;
        dRHOdy=(RHOY(i+1,j)-RHOY(i-1,j))/2/dy;
        % Left part
        L(kvy,kvy-Ny1*3)=ETA1/dx^2; % Vy1
        L(kvy,kvy-3)=ETAP1/dy^2; % Vy2
        L(kvy,kvy)=-(ETAP1+ETAP2)/dy^2-...
                      (ETA1+ETA2)/dx^2-...
                      dRHOdy*gy*dt-RHOY(i,j)/dt; % Vy3
        L(kvy,kvy+3)=ETAP2/dy^2; % Vy4
        L(kvy,kvy+Ny1*3)=ETA2/dx^2; % Vy5
        L(kvy,kvx)=ETAP1/dx/dy-ETA2/dx/dy-dRHOdx*gy*dt/4; %Vx3
        L(kvy,kvx+3)=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdx*gy*dt/4; %Vx4
        L(kvy,kvx-Ny1*3)=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdx*gy*dt/4; %Vx1
        L(kvy,kvx+3-Ny1*3)=ETAP2/dx/dy-ETA1/dx/dy-dRHOdx*gy*dt/4; %Vx2
        L(kvy,kpm)=pscale/dy; % P1
        L(kvy,kpm+3)=-pscale/dy; % P2
        
        % Right part
        R(kvy)=-RHOY(i,j)*(gy+vy0(i,j)/dt)-(SXY2-SXY1)/dx-(SYY2-SYY1)/dy;
        end
        
        % P equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1)
            % Boundary Condition
            % 1*P=0
            L(kpm,kpm)=1; % Left part
            R(kpm)=0; % Right part
        else
        % Internal points: compressible continuity eq.
        % dVx/dx+dVy/dy+BETTA*DP/Dt=ALPHA*DT/Dt
        %            Vy1
        %             |
        %        Vx1--P--Vx2
        %             |
        %            Vy2
        %
        % Left part
        L(kpm,kvx-Ny1*3)=-1/dx; % Vx1
        L(kpm,kvx)=1/dx; % Vx2
        L(kpm,kvy-3)=-1/dy; % Vy1
        L(kpm,kvy)=1/dy; % Vy2
        L(kpm,kpm)=pscale*BETTA(i,j)/dt; % P
        % Right part
        R(kpm)=BETTA(i,j)*pr0(i,j)/dt;
        end
        
    end
end

% 4) Solving matrixes, reloading solution
S=L\R; % Obtaining algebraic vector of solutions S()
% Reload solutions S() to vx(), vy(), p()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*3+1; % Vx
        kvy=kvx+1; % Vy
        kpm=kvx+2; % P
        % Reload solution
        vx(i,j)=S(kvx);
        vy(i,j)=S(kvy);
        pr(i,j)=S(kpm)*pscale;
    end
end
% Compute DVX,DVY
DVX=vx-vx0;
DVY=vy-vy0;


% Compute Stress, stress change and strain rate components
% Compute EPSILONxy, SIGMAxy in basic nodes
EXY=zeros(Ny,Nx); % Strain rate EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % Stress SIGMAxy, Pa
DSXY=zeros(Ny,Nx); % Stress change SIGMAxy, Pa
for j=1:1:Nx
    for i=1:1:Ny
        % EXY,SXY, DSXY
        EXY(i,j)=0.5*((vx(i+1,j)-vx(i,j))/dy+...
            (vy(i,j+1)-vy(i,j))/dx);
        SXY(i,j)=2*ETA(i,j)*EXY(i,j)*GGG(i,j)*dt/(GGG(i,j)*dt+ETA(i,j))+...
            SXY0(i,j)*ETA(i,j)/(GGG(i,j)*dt+ETA(i,j));
        DSXY(i,j)=SXY(i,j)-SXY0(i,j);
    end
end
% Compute EPSILONxx, SIGMA'xx in pressure nodes
EXX=zeros(Ny1,Nx1); % Strain rate EPSILONxx, 1/s
EII=zeros(Ny1,Nx1); % Second strain rate invariant, 1/s
SXX=zeros(Ny1,Nx1); % Stress SIGMA'xx, Pa
SII=zeros(Ny1,Nx1); % Second stress invariant, Pa
DSXX=zeros(Ny1,Nx1); % Stress change SIGMA'xx, Pa
DP=zeros(Ny1,Nx1); % Pressure change, Pa
for j=2:1:Nx
    for i=2:1:Ny
        % EXX
        EXX(i,j)=0.5*((vx(i,j)-vx(i,j-1))/dx-(vy(i,j)-vy(i-1,j))/dy);
        % SXX
        SXX(i,j)=2*ETAP(i,j)*EXX(i,j)*GGGP(i,j)*dt/(GGGP(i,j)*dt+ETAP(i,j))+...
            SXX0(i,j)*ETAP(i,j)/(GGGP(i,j)*dt+ETAP(i,j));
        DSXX(i,j)=SXX(i,j)-SXX0(i,j);
        % EII
        EII(i,j)=(EXX(i,j)^2+((EXY(i,j)+EXY(i-1,j)+...
            EXY(i,j-1)+EXY(i-1,j-1))/4)^2)^0.5;
        % SII
        SII(i,j)=(SXX(i,j)^2+((SXY(i,j)+SXY(i-1,j)+...
            SXY(i,j-1)+SXY(i-1,j-1))/4)^2)^0.5;
        % DP
        DP(i,j)=pr(i,j)-pr0(i,j);
    end
end

% Compute shear heating HS in Temperature/Pressure nodes
HS=zeros(Ny1,Nx1); % Shear heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        % Average SXY*EXY
        SXYEXY=(SXY(i,j)^2/ETA(i,j)+SXY(i-1,j)^2/ETA(i-1,j)+...
            SXY(i,j-1)^2/ETA(i,j-1)+SXY(i-1,j-1)^2/ETA(i-1,j-1))/4;
        % HS
        HS(i,j)=SXX(i,j)^2/ETAP(i,j)+SXYEXY;
    end
end



% Nodal adjusment
% External P-nodes: symmetry
pr(:,1)=pr(:,2);
pr(:,Nx1)=pr(:,Nx);
pr(1,:)=pr(2,:);
pr(Ny1,:)=pr(Ny,:);
pr0(:,1)=pr0(:,2);
pr0(:,Nx1)=pr0(:,Nx);
pr0(1,:)=pr0(2,:);
pr0(Ny1,:)=pr0(Ny,:);
SXX(:,1)=SXX(:,2);
SXX(:,Nx1)=SXX(:,Nx);
SXX(1,:)=SXX(2,:);
SXX(Ny1,:)=SXX(Ny,:);
SXX0(:,1)=SXX0(:,2);
SXX0(:,Nx1)=SXX0(:,Nx);
SXX0(1,:)=SXX0(2,:);
SXX0(Ny1,:)=SXX0(Ny,:);

% Apply criteria for decreasing the timestep
yndt=0;
% Define displacement timestep dtm
dtm=dt;
maxvx=max(max(abs(vx)));
maxvy=max(max(abs(vy)));
if(dtm*maxvx>dxymax*dx)
    dtm=dxymax*dx/maxvx/dtkoefv
    yndt=1;
end
if(dtm*maxvy>dxymax*dy)
    dtm=dxymax*dy/maxvy/dtkoefv
    yndt=1;
end
% Apply maximal stress change condition
dts=dt;
maxDScurrent=max(max(max(abs(DSXX))),max(max(abs(DSXY))));
if(maxDScurrent>DSmax)
    dts=dts/maxDScurrent*DSmax/dtkoefv
    yndt=1;
end

% Apply pressure timestepping condition
dtp=dt;
maxDPcurrent=max(max(abs(DP)));
if(maxDPcurrent>DPmax)
    dtp=dtp/maxDPcurrent*DPmax/dtkoefv
    yndt=1;
end

% Check if further dt reduction possible
if(dt<=dtmin)
    yndt=0;
end
    
    
% Reset timestep, restart iteration
ynstop=0;
if(yndt>0)
    dt=max(dtmin,min([dtm, dts, dtp]));
else
    % Exiting iterations
    break
end
% End iterations on Nodes
end
dt


% Reset displacement timestep
dtm=dt;


% Apply subgrid stress, pressure diffusion to markers
if(dsubgrids>0)
SXYSUM=zeros(Ny,Nx);
WTSUM=zeros(Ny,Nx);
SXXSUM=zeros(Ny1,Nx1);
PRSUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
    % SIGMA'xx, P
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node SIGMA'xx, P difference
    dsxxm0=sxxm(m)-(SXX0(i,j)*wtmij+SXX0(i+1,j)*wtmi1j+...
            SXX0(i,j+1)*wtmij1+SXX0(i+1,j+1)*wtmi1j1);
    dprm0=pr0m(m)-(pr0(i,j)*wtmij+pr0(i+1,j)*wtmi1j+...
            pr0(i,j+1)*wtmij1+pr0(i+1,j+1)*wtmi1j1);
    % Relax stress, P difference
    dsxxm1=dsxxm0*exp(-dsubgrids*dtm/(etam(tm(m))/gggm(tm(m))));
    dprm1=dprm0*exp(-dsubgrids*dtm/(etam(tm(m))/gggm(tm(m))));
    % Correct marker stress, P
    ddsxxm=dsxxm1-dsxxm0;
    sxxm(m)=sxxm(m)+ddsxxm;
    ddprm=dprm1-dprm0;
    pr0m(m)=pr0m(m)+ddprm;
    % Update subgrid diffusion on nodes
    % i,j Node
    SXXSUM(i,j)=SXXSUM(i,j)+ddsxxm*wtmij;
    PRSUM(i,j)=PRSUM(i,j)+ddprm*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    SXXSUM(i+1,j)=SXXSUM(i+1,j)+ddsxxm*wtmi1j;
    PRSUM(i+1,j)=PRSUM(i+1,j)+ddprm*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    SXXSUM(i,j+1)=SXXSUM(i,j+1)+ddsxxm*wtmij1;
    PRSUM(i,j+1)=PRSUM(i,j+1)+ddprm*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    SXXSUM(i+1,j+1)=SXXSUM(i+1,j+1)+ddsxxm*wtmi1j1;
    PRSUM(i+1,j+1)=PRSUM(i+1,j+1)+ddprm*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;

    % SIGMAxy
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node SIGMAxy difference
    dsxym0=sxym(m)-(SXY0(i,j)*wtmij+SXY0(i+1,j)*wtmi1j+...
            SXY0(i,j+1)*wtmij1+SXY0(i+1,j+1)*wtmi1j1);
    % Relax stress difference
    dsxym1=dsxym0*exp(-dsubgrids*dtm/(etam(tm(m))/gggm(tm(m))));
    % Correct marker stress
    ddsxym=dsxym1-dsxym0;
    sxym(m)=sxym(m)+ddsxym;
    % Update subgrid diffusion on nodes
    % i,j Node
    SXYSUM(i,j)=SXYSUM(i,j)+ddsxym*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    SXYSUM(i+1,j)=SXYSUM(i+1,j)+ddsxym*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    SXYSUM(i,j+1)=SXYSUM(i,j+1)+ddsxym*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    SXYSUM(i+1,j+1)=SXYSUM(i+1,j+1)+ddsxym*wtmi1j1;
    WTSUM(i+1,j+1)=WTSUM(i+1,j+1)+wtmi1j1;
end
% Compute DSXXsubgrid, DPsubgrid
DSXXsubgrid=zeros(Ny1,Nx1);
DPsubgrid=zeros(Ny1,Nx1);
% P-nodes
for j=2:1:Nx
    for i=2:1:Ny
        if(WTPSUM(i,j)>0)
            DSXXsubgrid(i,j)=SXXSUM(i,j)/WTPSUM(i,j);
            DPsubgrid(i,j)=PRSUM(i,j)/WTPSUM(i,j);
        end
    end
end
% Correct DSXX, P
DSXX=DSXX-DSXXsubgrid;
DP=DP-DPsubgrid;
% Compute DSXYsubgrid
DSXYsubgrid=zeros(Ny,Nx);
% Basic nodes
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            DSXYsubgrid(i,j)=SXYSUM(i,j)/WTSUM(i,j);
        end
    end
end
% Correct DSXY
DSXY=DSXY-DSXYsubgrid;
end

% Interpolate DSXX,P, DSXY to markers
for m=1:1:marknum
    % SIGMA'xx
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<2)
        j=2;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<2)
        i=2;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update marker by SIGMA'xx, P change 
    sxxm(m)=sxxm(m)+(DSXX(i,j)*wtmij+DSXX(i+1,j)*wtmi1j+...
            DSXX(i,j+1)*wtmij1+DSXX(i+1,j+1)*wtmi1j1);
    pr0m(m)=pr0m(m)+(DP(i,j)*wtmij+DP(i+1,j)*wtmi1j+...
            DP(i,j+1)*wtmij1+DP(i+1,j+1)*wtmi1j1);

    % SIGMAxy
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update marker by SIGMA'xx change 
    sxym(m)=sxym(m)+(DSXY(i,j)*wtmij+DSXY(i+1,j)*wtmi1j+...
            DSXY(i,j+1)*wtmij1+DSXY(i+1,j+1)*wtmi1j1);
end


% Apply subgrid velocity diffusion on markers
if(dsubgridv>0)
% Vx-nodes
RHOXSUM=zeros(Ny1,Nx1);
VXSUM=zeros(Ny1,Nx1);
% Vy-nodes
RHOYSUM=zeros(Ny1,Nx1);
VYSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
        
    % Interpolation from/to vx-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvx(1))/dx)+1;
    i=fix((ym(m)-yvx(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xvx(j);
    dymi=ym(m)-yvx(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node vx difference
    dvxm0=vx0m(m)-(vx0(i,j)*wtmij+vx0(i+1,j)*wtmi1j+...
            vx0(i,j+1)*wtmij1+vx0(i+1,j+1)*wtmi1j1);
    % Relax vx difference
    dvxm1=dvxm0*exp(-dsubgridv*etavpm(m)/rhom(tm(m))*dtm*(2/dx^2+2/dy^2));
    % Correct marker vx
    ddvxm=dvxm1-dvxm0;
    vx0m(m)=vx0m(m)+ddvxm;
    % Update properties
    % i,j Node
    RHOXSUM(i,j)=RHOXSUM(i,j)+rhom(tm(m))*wtmij;
    VXSUM(i,j)=VXSUM(i,j)+ddvxm*rhom(tm(m))*wtmij;
    % i+1,j Node
    RHOXSUM(i+1,j)=RHOXSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    VXSUM(i+1,j)=VXSUM(i+1,j)+ddvxm*rhom(tm(m))*wtmi1j;
    % i,j+1 Node
    RHOXSUM(i,j+1)=RHOXSUM(i,j+1)+rhom(tm(m))*wtmij1;
    VXSUM(i,j+1)=VXSUM(i,j+1)+ddvxm*rhom(tm(m))*wtmij1;
    % i+1,j+1 Node
    RHOXSUM(i+1,j+1)=RHOXSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    VXSUM(i+1,j+1)=VXSUM(i+1,j+1)+ddvxm*rhom(tm(m))*wtmi1j1;

    % Interpolation from/to vy-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvy(1))/dx)+1;
    i=fix((ym(m)-yvy(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xvy(j);
    dymi=ym(m)-yvy(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node vx difference
    dvym0=vy0m(m)-(vy0(i,j)*wtmij+vy0(i+1,j)*wtmi1j+...
            vy0(i,j+1)*wtmij1+vy0(i+1,j+1)*wtmi1j1);
    % Relax vx difference
    dvym1=dvym0*exp(-dsubgridv*etavpm(m)/rhom(tm(m))*dtm*(2/dx^2+2/dy^2));
    % Correct marker vx
    ddvym=dvym1-dvym0;
    vy0m(m)=vy0m(m)+ddvym;
    % Update properties
    % i,j Node
    RHOYSUM(i,j)=RHOYSUM(i,j)+rhom(tm(m))*wtmij;
    VYSUM(i,j)=VYSUM(i,j)+ddvym*rhom(tm(m))*wtmij;
    % i+1,j Node
    RHOYSUM(i+1,j)=RHOYSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    VYSUM(i+1,j)=VYSUM(i+1,j)+ddvym*rhom(tm(m))*wtmi1j;
    % i,j+1 Node
    RHOYSUM(i,j+1)=RHOYSUM(i,j+1)+rhom(tm(m))*wtmij1;
    VYSUM(i,j+1)=VYSUM(i,j+1)+ddvym*rhom(tm(m))*wtmij1;
    % i+1,j+1 Node
    RHOYSUM(i+1,j+1)=RHOYSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    VYSUM(i+1,j+1)=VYSUM(i+1,j+1)+ddvym*rhom(tm(m))*wtmi1j1;

end
% Compute DVXsubgrid
DVXsubgrid=zeros(Ny1,Nx1);
% Vx-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(RHOXSUM(i,j)>0)
            DVXsubgrid(i,j)=VXSUM(i,j)/RHOXSUM(i,j);
        end
    end
end
% Correct DVX
DVX=DVX-DVXsubgrid;
% Compute DVYsubgrid
DVYsubgrid=zeros(Ny1,Nx1);
% Vy-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(RHOYSUM(i,j)>0)
            DVYsubgrid(i,j)=VYSUM(i,j)/RHOYSUM(i,j);
        end
    end
end
% Correct DVX
DVY=DVY-DVYsubgrid;
end

% Interpolate DVX, DVY to markers
for m=1:1:marknum
    % Interpolation from vx-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvx(1))/dx)+1;
    i=fix((ym(m)-yvx(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xvx(j);
    dymi=ym(m)-yvx(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update marker vx
    vx0m(m)=vx0m(m)+(DVX(i,j)*wtmij+DVX(i+1,j)*wtmi1j+...
            DVX(i,j+1)*wtmij1+DVX(i+1,j+1)*wtmi1j1);

    % Interpolation from vy-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvy(1))/dx)+1;
    i=fix((ym(m)-yvy(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xvy(j);
    dymi=ym(m)-yvy(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update marker vy
    vy0m(m)=vy0m(m)+(DVY(i,j)*wtmij+DVY(i+1,j)*wtmi1j+...
            DVY(i,j+1)*wtmij1+DVY(i+1,j+1)*wtmi1j1);
end


% Compute velocity in pressure nodes
% vx
for j=2:1:Nx
    for i=2:1:Ny
        vxp(i,j)=(vx(i,j)+vx(i,j-1))/2;
    end
end
% Apply BC
% Top
vxp(1,2:Nx-1)=-bctop*vxp(2,2:Nx-1);    
% Bottom
vxp(Ny1,2:Nx-1)=-bcbottom*vxp(Ny,2:Nx-1);    
% Left
vxp(:,1)=-vxp(:,2);
% Right
vxp(:,Nx1)=-vxp(:,Nx);
% vy
for j=2:1:Nx
    for i=2:1:Ny
        vyp(i,j)=(vy(i,j)+vy(i-1,j))/2;
    end
end    
% Apply BC
% Left
vyp(2:Ny-1,1)=-bcleft*vyp(2:Ny-1,2);    
% Right
vyp(2:Ny-1,Nx1)=-bcright*vyp(2:Ny-1,Nx); % Free slip    
% Top
vyp(1,:)=-vyp(2,:);
% Bottom
vyp(Ny1,:)=-vyp(Ny,:);

% Compute rotation rate wyx=1/2(dVy/dx-dVx/dy) for basic nodes
for i=1:1:Ny
    for j=1:1:Nx
        wyx(i,j)=0.5*((vy(i,j+1)-vy(i,j))/dx-(vx(i+1,j)-vx(i,j))/dy);
    end
end

% Move markers with 4th order Runge-Kutta
vxm=zeros(4,1);
vym=zeros(4,1);
for m=1:1:marknum
    
    % Interpolate local rotation rate
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute vx velocity
    omegam=wyx(i,j)*wtmij+wyx(i+1,j)*wtmi1j+...
    wyx(i,j+1)*wtmij1+wyx(i+1,j+1)*wtmi1j1;
    % Analytical stress rotation using SIGMA'xx=-SIGMA'yy
    THETA=dtm*omegam; % Incremental rotation angle
    sxxmnew=sxxm(m)*cos(THETA)^2-sxxm(m)*sin(THETA)^2-sxym(m)*sin(2*THETA);
    sxymnew=sxxm(m)*sin(2*THETA)+sxym(m)*cos(2*THETA);
    sxxm(m)=sxxmnew; sxym(m)=sxymnew;    
    
    % Save initial marker's coordinates
    xA=xm(m);
    yA=ym(m);
    for rk=1:1:4
        % Interpolate vxp,vyp
        % Define i,j indexes for the upper left node
        j=fix((xm(m)-xp(1))/dx)+1;
        i=fix((ym(m)-yp(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % Compute distances
        dxmj=xm(m)-xp(j);
        dymi=ym(m)-yp(i);
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute vx, vy velocity
        vxm(rk)=vxp(i,j)*wtmij+vxp(i+1,j)*wtmi1j+...
            vxp(i,j+1)*wtmij1+vxp(i+1,j+1)*wtmi1j1;
        vym(rk)=vyp(i,j)*wtmij+vyp(i+1,j)*wtmi1j+...
            vyp(i,j+1)*wtmij1+vyp(i+1,j+1)*wtmi1j1;
        
        % Interpolate vx
        % Define i,j indexes for the upper left node
        j=fix((xm(m)-xvx(1))/dx)+1;
        i=fix((ym(m)-yvx(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % Compute distances
        dxmj=xm(m)-xvx(j);
        dymi=ym(m)-yvx(i);
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute vx velocity
        vxm(rk)=vpratio*vxm(rk)+(1-vpratio)*(vx(i,j)*wtmij+vx(i+1,j)*wtmi1j+...
            vx(i,j+1)*wtmij1+vx(i+1,j+1)*wtmi1j1);
        
        % Interpolate vy
        % Define i,j indexes for the upper left node
        j=fix((xm(m)-xvy(1))/dx)+1;
        i=fix((ym(m)-yvy(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
        end
        % Compute distances
        dxmj=xm(m)-xvy(j);
        dymi=ym(m)-yvy(i);
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute vx velocity
        vym(rk)=vpratio*vym(rk)+(1-vpratio)*(vy(i,j)*wtmij+vy(i+1,j)*wtmi1j+...
            vy(i,j+1)*wtmij1+vy(i+1,j+1)*wtmi1j1);        
        
        % Change coordinates to obtain B,C,D points
        if(rk==1 || rk==2)
            xm(m)=xA+dtm/2*vxm(rk);
            ym(m)=yA+dtm/2*vym(rk);
        elseif(rk==3)
            xm(m)=xA+dtm*vxm(rk);
            ym(m)=yA+dtm*vym(rk);
        end
    end
    % Restore initial coordinates
    xm(m)=xA;
    ym(m)=yA;
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Move markers
    xm(m)=xm(m)+dtm*vxmeff;
    ym(m)=ym(m)+dtm*vymeff;
end  


% Update timesum
timesum=timesum+dtm;

% Tracing point position
xtrace(timestep,1)=xm(marknum); % x-coordinate
ytrace(timestep,1)=ym(marknum); % y-coordinate
ttrace(timestep,1)=timesum; % time-coordinate


% Visualise results
if(timestep==1 || fix(timestep/visstep)*visstep==timestep)

figure(1);colormap('Jet');clf
subplot(3,4,1)
pcolor(x,y,log10(ETA));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['log10ETA, Pa*s timestep=',num2str(timestep)])

subplot(3,4,2)
pcolor(xp,yp,pr)
shading interp;
axis ij image;
colorbar
title('Pressure, Pa')

subplot(3,4,3)
pcolor(xp,yp,vxp)
shading interp;
axis ij image;
colorbar
title(['Vx, m/s time=',num2str(timesum),' s'])
hold on
quiver(xp(3:5:Nx1),yp(3:5:Ny1),vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(3,4,4)
pcolor(xp,yp,vyp)
shading interp;
axis ij image;
colorbar
title('Vy, m/s')
hold on
quiver(xp(3:5:Nx1),yp(3:5:Ny1),vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(3,4,5)
pcolor(xp,yp,HS)
shading interp;
axis ij image;
colorbar
title('Shear heating, W/m^3')

subplot(3,4,6)
pcolor(xp,yp,EXX)
shading interp;
axis ij image;
colorbar
title('EPSILONxx, 1/s')

subplot(3,4,7)
pcolor(x,y,EXY)
shading interp;
axis ij image;
colorbar
title('EPSILONxy, 1/s')

subplot(3,4,8)
pcolor(xp,yp,log10(EII))
shading interp;
axis ij image;
colorbar
title('logEPSILONii, 1/s')

subplot(3,4,9)
pcolor(xp,yp,RHO)
shading interp;
axis ij image;
colorbar
title('Density, kg/m^3')

subplot(3,4,10)
pcolor(xp,yp,SXX)
shading interp;
axis ij image;
colorbar
title('SIGMAxx, Pa')

subplot(3,4,11)
pcolor(x,y,SXY)
shading interp;
axis ij image;
colorbar
title('SIGMAxy, Pa')

subplot(3,4,12)
pcolor(xp,yp,SII)
shading interp;
axis ij image;
colorbar
title('SIGMAii, Pa')
% Save figure
namefig=['Ball1_1_' num2str(timestep)];
print ('-djpeg', '-r150',namefig)


figure(2);colormap('Jet');clf
pcolor(x,y,log10(ETA));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['log10ETA, Pa*s timestep=',num2str(timestep),' time=',num2str(timesum),' s'])
hold on
quiver(xp(3:5:Nx1),yp(3:5:Ny1),vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'w','LineWidth',1.5)
axis([0 xsize 0 ysize]);
% Save figure
namefig=['Ball1_2_' num2str(timestep)];
print ('-djpeg', '-r150',namefig)

figure(3);colormap('Jet');clf
subplot(2,1,1)
plot(ttrace(1:timestep),ytrace(1:timestep));% caxis([17 21])
axis ij
title(['y position vs. time timestep=',num2str(timestep),' time=',num2str(timesum),' s'])
subplot(2,1,2)
plot(ytrace(1:timestep));% caxis([17 21])
axis ij
title(['y position vs. timestep timestep=',num2str(timestep),' time=',num2str(timesum),' s'])
% Save figure
namefig=['Ball1_3_' num2str(timestep)];
print ('-djpeg', '-r150',namefig)



pause(0.1)
end

end
