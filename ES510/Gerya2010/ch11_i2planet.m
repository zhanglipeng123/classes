% Modeling thermomechanical evolution of a self-gravitating planet in 2D
% Solving Stokes, continuity, temperature and Poisson eq.
% in primitive variable formulation
% with variable viscosity and thermal conductivity
% using FD with staggered grid

% Clearing memory and figures
clear all; clf

% Define Numerical model
xsize=10000000; % Horizontal model size, m
ysize=10000000; % Vertical model size, m
Nx=101; % Horizontal grid resolution
Ny=101; % Vertical grid resolution
Nx1=Nx+1;
Ny1=Ny+1;
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step, m

G=6.672e-11; % Gravity constant, N*m^2/kg^2

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
EXY=zeros(Ny,Nx); % EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % SIGMAxy, 1/s
TYP=zeros(Ny,Nx); % material type
% Vx-Nodes
RHOX=zeros(Ny1,Nx1); % Density, kg/m^3
KX=zeros(Ny1,Nx1); % Thermal conductivity, W/m/K
vx=zeros(Ny1,Nx1); % vx-velocity, m/s
gx=zeros(Ny1,Nx1); % gx-gravity, m/s^2
% Vy-Nodes
RHOY=zeros(Ny1,Nx1); % Density, kg/m^3
KY=zeros(Ny1,Nx1); % Thermal conductivity, W/m/K
vy=zeros(Ny1,Nx1); % vy-velocity, m/s
gy=zeros(Ny1,Nx1); % gy-gravity, m/s^2
% P-nodes
RHO=zeros(Ny1,Nx1); % Density, kg/m^3
RHOCP=zeros(Ny1,Nx1); % Volumetric heat capacity, J/m^3/K
ALPHA=zeros(Ny1,Nx1); % Thermal expansion, J/m^3/K
HR=zeros(Ny1,Nx1); % Radioactive heating, W/m^3
HA=zeros(Ny1,Nx1); % Adiabatic heating, W/m^3
HS=zeros(Ny1,Nx1); % Shear heating, W/m^3
ETAP=zeros(Ny1,Nx1); % Viscosity, Pa*s
EXX=zeros(Ny,Nx); % EPSILONxx, 1/s
SXX=zeros(Ny,Nx); % SIGMAxx, 1/s
tk1=zeros(Ny1,Nx1); % Old temperature, K
tk2=zeros(Ny1,Nx1); % New temperature, K
vxp=zeros(Ny1,Nx1); % Vx in pressure nodes, m/s
vyp=zeros(Ny1,Nx1); % Vy in pressure nodes, m/s
pr=zeros(Ny1,Nx1); % Pressure, Pa
PHI=zeros(Ny1,Nx1); % Gravity potential, J/kg

% Define markers
Nxm=(Nx-1)*4; % Marker grid resolution in horizontal direction
Nym=(Ny-1)*4; % Marker grid resolution in vertical direction
dxm=xsize/Nxm; % Marker grid step in horizontal direction,m
dym=ysize/Nym; % Marker grid step in vertical direction,m
marknum=Nxm*Nym; % Number of markers
xm=zeros(1,marknum); % Horizontal coordinates, m
ym=zeros(1,marknum); % Vertical coordinates, m
tm=zeros(1,marknum); % Material type
tkm=zeros(1,marknum); % Marker temperature, K

% Define properties of materials: 
%        space mantle1 mantle2 plume1 plume2
rhom   = [0      3300   3300   3000   3000  ]; % Density, kg/m^3
etam   = [1e+17  1e+21  1e+21  1e+19  1e+19 ]; % Viscosity, Pa s
rhocpm = [3.3e+6 3.3e+6 3.3e+6 3.0e+6 3.0e+6]; % Volumetric heat capacity, kg/m^3
alpham = [0      3e-5   3e-5   2e-5   2e-5  ]; % Thermal expansion, 1/K
km     = [3000   3      3      2      2     ]; % Thermal conductivity, W/m/K
hrm    = [0      2e-8   2e-8   2e-7   2e-7  ]; % Radiogenic heat production, W/m^3

% Define marker coordinates, temperature and material type
m=1; % Marker counter
for jm=1:1:Nxm
    for im=1:1:Nym
        % Define marker coordinates
        xm(m)=dxm/2+(jm-1)*dxm+(rand-0.5)*dxm;
        ym(m)=dym/2+(im-1)*dym+(rand-0.5)*dym;
        % Marker properties
        % Distance from the model centre
        rmark=((xm(m)-xsize/2)^2+(ym(m)-ysize/2)^2)^0.5;
        % Sticky space (to have internal free surface)
%         if(rmark<4800000)
            tm(m)=1; % Material type
            tkm(m)=273; % Temperature
%         else
%             % Rigid shell
%             tm(m)=4; % Material type
%             tkm(m)=273; % Temperature
%         end
        % Square Planet
        tcur=3;
        for mdist=3000000:-300000:1000000
            if(abs(xm(m)-xsize/2)<mdist && abs(ym(m)-ysize/2)<mdist)
                tm(m)=tcur; % Material type
                tkm(m)=1500; % Temperature
            end
            % Change type
            tcur=5-tcur;
        end
%         rmark=((xm(m)-xsize*0.35)^2+(ym(m)-ysize*0.35)^2)^0.5;
        % Low-dinsity core
        tcur=4;
        for mdist=1500000:-300000:300000
%             if(rmark<mdist)
%                 tm(m)=tcur; % Material type
%                 tkm(m)=1800; % Temperature
%             end
            if(abs(xm(m)-xsize/2)<mdist && abs(ym(m)-ysize/2)<mdist)
                tm(m)=tcur; % Material type
                tkm(m)=1800; % Temperature
            end
            % Change type
            tcur=9-tcur;
        end
        % Update marker counter
        m=m+1;
    end
end

% Introducing scaled pressure
pscale=1e+21/dx;

% Define global matrixes 
% Mechanical solution: L(), R()
N=Nx1*Ny1*3; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts
% Thermal solution: LT(), RT()
N=Nx1*Ny1; % Global number of unknowns
LT=sparse(N,N); % Matrix of coefficients (left part)
RT=zeros(N,1); % Vector of right parts
% Gravity solution: LP(), RP()
N=Nx1*Ny1; % Global number of unknowns
LP=sparse(N,N); % Matrix of coefficients (left part)
RP=zeros(N,1); % Vector of right parts

% Mechanical boundary conditions: free slip=-1; No Slip=1
bcleft=-1;
bcright=-1;
bctop=-1;
bcbottom=-1;
% Thermal boundary conditions
tktop=273;
tkbottom=1500;

% Timestepping
dt=1e+10; % initial timestep
dtkoef=1.2; % timestep increment
dxymax=0.5; % max marker movement per time step, grid steps
vpratio=1/3; % Weight of averaged velocity for moving markers
DTmax=20; % max temperature change per time step, K
dsubgridt=1; % subgrid diffusion parameter
timesum=0; % Model time, s
for timestep=1:1:2000
    
% Interpolate properties from markers to nodes
% Basic nodes
ETASUM=zeros(Ny,Nx);
TYPSUM=zeros(Ny,Nx);
WTSUM=zeros(Ny,Nx);
% Vx-nodes
RHOXSUM=zeros(Ny1,Nx1);
KXSUM=zeros(Ny1,Nx1);
WTXSUM=zeros(Ny1,Nx1);
% Vy-nodes
RHOYSUM=zeros(Ny1,Nx1);
KYSUM=zeros(Ny1,Nx1);
WTYSUM=zeros(Ny1,Nx1);
% P-Nodes
ETAPSUM=zeros(Ny1,Nx1);
RHOSUM=zeros(Ny1,Nx1);
RHOCPSUM=zeros(Ny1,Nx1);
ALPHASUM=zeros(Ny1,Nx1);
HRSUM=zeros(Ny1,Nx1);
TKSUM=zeros(Ny1,Nx1);
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
    TYPSUM(i,j)=TYPSUM(i,j)+tm(m)*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    ETASUM(i+1,j)=ETASUM(i+1,j)+etam(tm(m))*wtmi1j;
    TYPSUM(i+1,j)=TYPSUM(i+1,j)+tm(m)*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    ETASUM(i,j+1)=ETASUM(i,j+1)+etam(tm(m))*wtmij1;
    TYPSUM(i,j+1)=TYPSUM(i,j+1)+tm(m)*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    ETASUM(i+1,j+1)=ETASUM(i+1,j+1)+etam(tm(m))*wtmi1j1;
    TYPSUM(i+1,j+1)=TYPSUM(i+1,j+1)+tm(m)*wtmi1j1;
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
    KXSUM(i,j)=KXSUM(i,j)+km(tm(m))*wtmij;
    WTXSUM(i,j)=WTXSUM(i,j)+wtmij;
    % i+1,j Node
    RHOXSUM(i+1,j)=RHOXSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    KXSUM(i+1,j)=KXSUM(i+1,j)+km(tm(m))*wtmi1j;
    WTXSUM(i+1,j)=WTXSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOXSUM(i,j+1)=RHOXSUM(i,j+1)+rhom(tm(m))*wtmij1;
    KXSUM(i,j+1)=KXSUM(i,j+1)+km(tm(m))*wtmij1;
    WTXSUM(i,j+1)=WTXSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOXSUM(i+1,j+1)=RHOXSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    KXSUM(i+1,j+1)=KXSUM(i+1,j+1)+km(tm(m))*wtmi1j1;
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
    KYSUM(i,j)=KYSUM(i,j)+km(tm(m))*wtmij;
    WTYSUM(i,j)=WTYSUM(i,j)+wtmij;
    % i+1,j Node
    RHOYSUM(i+1,j)=RHOYSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    KYSUM(i+1,j)=KYSUM(i+1,j)+km(tm(m))*wtmi1j;
    WTYSUM(i+1,j)=WTYSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOYSUM(i,j+1)=RHOYSUM(i,j+1)+rhom(tm(m))*wtmij1;
    KYSUM(i,j+1)=KYSUM(i,j+1)+km(tm(m))*wtmij1;
    WTYSUM(i,j+1)=WTYSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOYSUM(i+1,j+1)=RHOYSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    KYSUM(i+1,j+1)=KYSUM(i+1,j+1)+km(tm(m))*wtmi1j1;
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
    RHOSUM(i,j)=RHOSUM(i,j)+rhom(tm(m))*wtmij;
    RHOCPSUM(i,j)=RHOCPSUM(i,j)+rhocpm(tm(m))*wtmij;
    ALPHASUM(i,j)=ALPHASUM(i,j)+alpham(tm(m))*wtmij;
    HRSUM(i,j)=HRSUM(i,j)+hrm(tm(m))*wtmij;
    TKSUM(i,j)=TKSUM(i,j)+tkm(m)*rhocpm(tm(m))*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    ETAPSUM(i+1,j)=ETAPSUM(i+1,j)+etam(tm(m))*wtmi1j;
    RHOSUM(i+1,j)=RHOSUM(i+1,j)+rhom(tm(m))*wtmi1j;
    RHOCPSUM(i+1,j)=RHOCPSUM(i+1,j)+rhocpm(tm(m))*wtmi1j;
    ALPHASUM(i+1,j)=ALPHASUM(i+1,j)+alpham(tm(m))*wtmi1j;
    HRSUM(i+1,j)=HRSUM(i+1,j)+hrm(tm(m))*wtmi1j;
    TKSUM(i+1,j)=TKSUM(i+1,j)+tkm(m)*rhocpm(tm(m))*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    ETAPSUM(i,j+1)=ETAPSUM(i,j+1)+etam(tm(m))*wtmij1;
    RHOSUM(i,j+1)=RHOSUM(i,j+1)+rhom(tm(m))*wtmij1;
    RHOCPSUM(i,j+1)=RHOCPSUM(i,j+1)+rhocpm(tm(m))*wtmij1;
    ALPHASUM(i,j+1)=ALPHASUM(i,j+1)+alpham(tm(m))*wtmij1;
    HRSUM(i,j+1)=HRSUM(i,j+1)+hrm(tm(m))*wtmij1;
    TKSUM(i,j+1)=TKSUM(i,j+1)+tkm(m)*rhocpm(tm(m))*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    ETAPSUM(i+1,j+1)=ETAPSUM(i+1,j+1)+etam(tm(m))*wtmi1j1;
    RHOSUM(i+1,j+1)=RHOSUM(i+1,j+1)+rhom(tm(m))*wtmi1j1;
    RHOCPSUM(i+1,j+1)=RHOCPSUM(i+1,j+1)+rhocpm(tm(m))*wtmi1j1;
    ALPHASUM(i+1,j+1)=ALPHASUM(i+1,j+1)+alpham(tm(m))*wtmi1j1;
    HRSUM(i+1,j+1)=HRSUM(i+1,j+1)+hrm(tm(m))*wtmi1j1;
    TKSUM(i+1,j+1)=TKSUM(i+1,j+1)+tkm(m)*rhocpm(tm(m))*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;
end
% Compute physical properties
% Basic nodes
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            ETA(i,j)=ETASUM(i,j)/WTSUM(i,j);
            TYP(i,j)=TYPSUM(i,j)/WTSUM(i,j);
        end
    end
end
% Vx-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTXSUM(i,j)>0)
            RHOX(i,j)=RHOXSUM(i,j)/WTXSUM(i,j);
            KX(i,j)=KXSUM(i,j)/WTXSUM(i,j);
        end
    end
end
% Vy-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTYSUM(i,j)>0)
            RHOY(i,j)=RHOYSUM(i,j)/WTYSUM(i,j);
            KY(i,j)=KYSUM(i,j)/WTYSUM(i,j);
        end
    end
end
% P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM(i,j)>0)
            ETAP(i,j)=ETAPSUM(i,j)/WTPSUM(i,j);
            RHO(i,j)=RHOSUM(i,j)/WTPSUM(i,j);
            RHOCP(i,j)=RHOCPSUM(i,j)/WTPSUM(i,j);
            ALPHA(i,j)=ALPHASUM(i,j)/WTPSUM(i,j);
            HR(i,j)=HRSUM(i,j)/WTPSUM(i,j);
            tk1(i,j)=TKSUM(i,j)/RHOCPSUM(i,j);
        end
    end
end
% Applying thermal boundary conditions for interpolated temperature
% Upper boundary 
tk1(1,2:Nx)=tk1(2,2:Nx); % Insulating boundary
% Lower boundary 
tk1(Ny1,2:Nx)=tk1(Ny,2:Nx); % Insulating boundary
% Left boundary
tk1(:,1)=tk1(:,2); % Insulating boundary
% Right boundary
tk1(:,Nx1)=tk1(:,Nx); % Insulating boundary

% Gravity solution
% Composing global matrixes LT(), RT()
% Going through all points of the 2D grid and
% composing respective equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global index in algebraic space
        gk=(j-1)*Ny1+i;
        % Distance from the model centre
        rnode=((xp(j)-xsize/2)^2+(yp(i)-ysize/2)^2)^0.5;
        % External points
        if(rnode>xsize/2 || i==1 || i==Ny1 || j==1 || j==Nx1)
            % Boundary Condition
            % PHI=0
            LP(gk,gk)=1; % Left part
            RP(gk)=0; % Right part
        else
        % Internal points: Temperature eq.
        % d2PHI/dx^2+d2PHI/dy^2=2/3*4*G*pi*RHO
        %          PHI2
        %           |
        %           |
        %  PHI1----PHI3----PHI5
        %           |
        %           |
        %          PHI4
        %
        % Density gradients
        dRHOdx=(RHO(i,j+1)-RHO(i,j-1))/2/dx;
        dRHOdy=(RHO(i+1,j)-RHO(i-1,j))/2/dy;
        % Left part
        LP(gk,gk-Ny1)=1/dx^2; % PHI1
        LP(gk,gk-1)=1/dy^2; % PHI2
        LP(gk,gk)=-2/dx^2-2/dy^2; % PHI3
        LP(gk,gk+1)=1/dy^2; % PHI4
        LP(gk,gk+Ny1)=1/dx^2; % PHI5
        % Right part
        RP(gk)=2/3*4*G*pi*RHO(i,j);
        end
    end
end

% Solving matrixes
SP=LP\RP; % Obtaining algebraic vector of solutions SP()

% Reload solutions SP() to geometrical array PHI()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Compute global index
        gk=(j-1)*Ny1+i;
        % Reload solution
        PHI(i,j)=SP(gk);
    end
end
% Compute gravity acceleration
% gx
for j=1:1:Nx
    for i=1:1:Ny1
        % gx=-dPHI/dx
        gx(i,j)=-(PHI(i,j+1)-PHI(i,j))/dx;
    end
end
% gy
for j=1:1:Nx1
    for i=1:1:Ny
        % gy=-dPHI/dy
        gy(i,j)=-(PHI(i+1,j)-PHI(i,j))/dy;
    end
end


% Mechanical Solution
% Composing global matrixes L(), R() for Stokes and continuity equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*3+1; % Vx
        kvy=kvx+1; % Vy
        kpm=kvx+2; % P
        
        % Vx equation External points
        if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
            % Boundary Condition
            % 1*Vx=0
            L(kvx,kvx)=1; % Left part
            R(kvx)=0; % Right part
            % Top boundary
            if(i==1 && j>1 && j<Nx)
                L(kvx,kvx+3)=bctop; % Left part
            end
            % Bottom boundary
            if(i==Ny1 && j>1 && j<Nx)
                L(kvx,kvx-3)=bcbottom; % Left part
            end
        else
        % Internal points: x-Stokes eq.
        % ETA*(d2Vx/dx^2+d2Vx/dy^2)-dP/dx=0
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
        % Viscosity points
        ETA1=ETA(i-1,j);
        ETA2=ETA(i,j);
        ETAP1=ETAP(i,j);
        ETAP2=ETAP(i,j+1);
        % Density gradients
        dRHOdx=(RHOX(i,j+1)-RHOX(i,j-1))/2/dx;
        dRHOdy=(RHOX(i+1,j)-RHOX(i-1,j))/2/dy;
        % Left part
        L(kvx,kvx-Ny1*3)=2*ETAP1/dx^2; % Vx1
        L(kvx,kvx-3)=ETA1/dy^2; % Vx2
        L(kvx,kvx)=-2*(ETAP1+ETAP2)/dx^2-...
                      (ETA1+ETA2)/dy^2-...
                      dRHOdx*gx(i,j)*dt; % Vx3
        L(kvx,kvx+3)=ETA2/dy^2; % Vx4
        L(kvx,kvx+Ny1*3)=2*ETAP2/dx^2; % Vx5
        L(kvx,kvy)=-ETA2/dx/dy-dRHOdy*gx(i,j)*dt/4;  % Vy2
        L(kvx,kvy+Ny1*3)=ETA2/dx/dy-dRHOdy*gx(i,j)*dt/4;  % Vy4
        L(kvx,kvy-3)=ETA1/dx/dy-dRHOdy*gx(i,j)*dt/4;  % Vy1
        L(kvx,kvy+Ny1*3-3)=-ETA1/dx/dy-dRHOdy*gx(i,j)*dt/4;  % Vy3
        L(kvx,kpm)=pscale/dx; % P1
        L(kvx,kpm+Ny1*3)=-pscale/dx; % P2
        % Right part
        R(kvx)=-RHOX(i,j)*gx(i,j);
        end
        
        % Vy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            % Boundary Condition
            % 1*Vy=0
            L(kvy,kvy)=1; % Left part
            R(kvy)=0; % Right part
            % Left boundary
            if(j==1 && i>1 && i<Ny)
                L(kvy,kvy+3*Ny1)=bcleft; % Left part
            end
            % Right boundary
            if(j==Nx1 && i>1 && i<Ny)
                L(kvy,kvy-3*Ny1)=bcright; % Left part
            end
        else
        % Internal points: y-Stokes eq.
        % ETA*(d2Vy/dx^2+d2Vy/dy^2)-dP/dy=-RHO*gy
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
        % Viscosity points
        % Viscosity points
        ETA1=ETA(i,j-1);
        ETA2=ETA(i,j);
        ETAP1=ETAP(i,j);
        ETAP2=ETAP(i+1,j);
        % Density gradients
        dRHOdx=(RHOY(i,j+1)-RHOY(i,j-1))/2/dx;
        dRHOdy=(RHOY(i+1,j)-RHOY(i-1,j))/2/dy;
        % Left part
        L(kvy,kvy-Ny1*3)=ETA1/dx^2; % Vy1
        L(kvy,kvy-3)=2*ETAP1/dy^2; % Vy2
        L(kvy,kvy)=-2*(ETAP1+ETAP2)/dy^2-...
                      (ETA1+ETA2)/dx^2-...
                      dRHOdy*gy(i,j)*dt; % Vy3
        L(kvy,kvy+3)=2*ETAP2/dy^2; % Vy4
        L(kvy,kvy+Ny1*3)=ETA2/dx^2; % Vy5
        L(kvy,kvx)=-ETA2/dx/dy-dRHOdx*gy(i,j)*dt/4; %Vx3
        L(kvy,kvx+3)=ETA2/dx/dy-dRHOdx*gy(i,j)*dt/4; %Vx4
        L(kvy,kvx-Ny1*3)=ETA1/dx/dy-dRHOdx*gy(i,j)*dt/4; %Vx1
        L(kvy,kvx+3-Ny1*3)=-ETA1/dx/dy-dRHOdx*gy(i,j)*dt/4; %Vx2
        L(kvy,kpm)=pscale/dy; % P1
        L(kvy,kpm+3)=-pscale/dy; % P2
        
        % Right part
        R(kvy)=-RHOY(i,j)*gy(i,j);
        end
        
        % P equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1 ||...
          (i==2 && j==2))
            % Boundary Condition
            % 1*P=0
            L(kpm,kpm)=1; % Left part
            R(kpm)=0; % Right part
            % Real BC
            if(i==2 && j==2)
                L(kpm,kpm)=1*pscale; %Left part
                R(kpm)=1e+9; % Right part
            end
        else
        % Internal points: continuity eq.
        % dVx/dx+dVy/dy=0
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
        % Right part
        R(kpm)=0;
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

% Compute Stress and strain rate components
% Compute EPSILONxy, SIGMAxy in basic nodes
EXY=zeros(Ny,Nx); % Strain rate EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % Strain rate SIGMAxy, Pa
for j=1:1:Nx
    for i=1:1:Ny
        % EXY,SXY
        EXY(i,j)=0.5*((vx(i+1,j)-vx(i,j))/dy+...
            (vy(i,j+1)-vy(i,j))/dx);
        SXY(i,j)=2*ETA(i,j)*EXY(i,j);
    end
end
% Compute EPSILONxx, SIGMA'xx in pressure nodes
EXX=zeros(Ny1,Nx1); % Strain rate EPSILONxx, 1/s
SXX=zeros(Ny1,Nx1); % Strain rate SIGMAxx, Pa
for j=2:1:Nx
    for i=2:1:Ny
        % EXX
        EXX(i,j)=(vx(i,j)-vx(i,j-1))/dx;
        % SXX
        SXX(i,j)=2*ETAP(i,j)*EXX(i,j);
    end
end

% Compute shear heating HS in pressure nodes
HS=zeros(Ny1,Nx1); % Adiabatic heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        % Average SXY*EXY
        SXYEXY=(SXY(i,j)*EXY(i,j)+SXY(i-1,j)*EXY(i-1,j)+...
            SXY(i,j-1)*EXY(i,j-1)+SXY(i-1,j-1)*EXY(i-1,j-1))/4;
        % HS
        HS(i,j)=2*SXX(i,j)*EXX(i,j)+2*SXYEXY;
    end
end

% Compute adiabatic heating HA in pressure nodes
HA=zeros(Ny1,Nx1); % Shear heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        % Average vy, vx
        VXP=(vx(i,j)+vx(i,j-1))/2;
        VYP=(vy(i,j)+vy(i-1,j))/2;
        % Average gy, gx
        GXP=(gx(i,j)+gx(i,j-1))/2;
        GYP=(gy(i,j)+gy(i-1,j))/2;
        % HA
        HA(i,j)=tk1(i,j)*ALPHA(i,j)*RHO(i,j)*(VXP*GXP+VYP*GYP);
    end
end

% Define timestep
dt=dt*dtkoef;
maxvx=max(max(abs(vx)));
maxvy=max(max(abs(vy)));
if(dt*maxvx>dxymax*dx)
    dt=dxymax*dx/maxvx;
end
if(dt*maxvy>dxymax*dy)
    dt=dxymax*dy/maxvy;
end


% Thermal iterations
for titer=1:1:2
% Composing global matrixes LT(), RT()
% Going through all points of the 2D grid and
% composing respective equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global index in algebraic space
        gk=(j-1)*Ny1+i;
        % External points
        if(i==1 || i==Ny1 || j==1 || j==Nx1)
            % Boundary Condition
            % Top BC: dT/dx=0
            if(i==1 && j>1 && j<Nx1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk+1)=-1; % Left part
                RT(gk)=0; % Right part
            end
            % Bottom BC: dT/dy=0
            if(i==Ny1 && j>1 && j<Nx1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk-1)=-1; % Left part
                RT(gk)=0; % Right part
            end
            % Left BC: dT/dx=0
            if(j==1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk+Ny1)=-1; % Left part
                RT(gk)=0; % Right part
            end
            % Right BC: dT/dx=0
            if(j==Nx1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk-Ny1)=-1; % Left part
                RT(gk)=0; % Right part
            end
        else
        % Internal points: Temperature eq.
        % RHO*CP*dT/dt=-dqx/dx-dqy/dy+Hr+Hs+Ha
        %          Tdt2
        %           |
        %          Ky1
        %           |
        %Tdt1-Kx1-T03,Tdt3-Kx2-Tdt5
        %           |
        %          Ky2
        %           |
        %          Tdt4
        %
        % Left part
        Kx1=KX(i,j-1); 
        Kx2=KX(i,j); 
        Ky1=KY(i-1,j); 
        Ky2=KY(i,j); 
        LT(gk,gk-Ny1)=-Kx1/dx^2; % T1
        LT(gk,gk-1)=-Ky1/dy^2; % FI2
        LT(gk,gk)=RHOCP(i,j)/dt+(Kx1+Kx2)/dx^2+(Ky1+Ky2)/dy^2; % FI3
        LT(gk,gk+1)=-Ky2/dy^2; % FI4
        LT(gk,gk+Ny1)=-Kx2/dx^2; % FI5
        % Right part
        RT(gk)=RHOCP(i,j)/dt*tk1(i,j)+HR(i,j)+HA(i,j)+HS(i,j);
        end
    end
end

% Solving matrixes
ST=LT\RT; % Obtaining algebraic vector of solutions ST()

% Reload solutions ST() to geometrical array Tdt()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Compute global index
        gk=(j-1)*Ny1+i;
        % Reload solution
        tk2(i,j)=ST(gk);
    end
end
% Compute DT
DT=tk2-tk1;
titer
dt
% Apply thermal timestepping condition
maxDTcurrent=max(max(abs(DT)));
if(titer<2 && maxDTcurrent>DTmax)
    dt=dt/maxDTcurrent*DTmax;
else
    break;
end
end
DT0=DT;


% Apply subgrid diffusion on markers
if(dsubgridt>0)
TKSUM=zeros(Ny1,Nx1);
RHOCPSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
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
    % Compute marker-node T difference
    dtm0=tkm(m)-(tk1(i,j)*wtmij+tk1(i+1,j)*wtmi1j+...
            tk1(i,j+1)*wtmij1+tk1(i+1,j+1)*wtmi1j1);
    % Relax temperature difference
    dtmdt=dtm0*exp(-dsubgridt*km(tm(m))*dt/rhocpm(tm(m))*(2/dx^2+2/dy^2));
    % Correct marker temperature
    ddtm=dtmdt-dtm0;
    tkm(m)=tkm(m)+ddtm;
    % Update subgrid diffusion on nodes
    % i,j Node
    TKSUM(i,j)=TKSUM(i,j)+ddtm*rhocpm(tm(m))*wtmij;
    RHOCPSUM(i,j)=RHOCPSUM(i,j)+rhocpm(tm(m))*wtmij;
    % i+1,j Node
    TKSUM(i+1,j)=TKSUM(i+1,j)+ddtm*rhocpm(tm(m))*wtmi1j;
    RHOCPSUM(i+1,j)=RHOCPSUM(i+1,j)+rhocpm(tm(m))*wtmi1j;
    % i,j+1 Node
    TKSUM(i,j+1)=TKSUM(i,j+1)+ddtm*rhocpm(tm(m))*wtmij1;
    RHOCPSUM(i,j+1)=RHOCPSUM(i,j+1)+rhocpm(tm(m))*wtmij1;
    % i+1,j+1 Node
    TKSUM(i+1,j+1)=TKSUM(i+1,j+1)+ddtm*rhocpm(tm(m))*wtmi1j1;
    RHOCPSUM(i+1,j+1)=RHOCPSUM(i+1,j+1)+rhocpm(tm(m))*wtmi1j1;
end
% Compute DTsubgrid
DTsubgrid=zeros(Ny1,Nx1);
% P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(RHOCPSUM(i,j)>0)
            DTsubgrid(i,j)=TKSUM(i,j)/RHOCPSUM(i,j);
        end
    end
end
% Correct DT
DT=DT-DTsubgrid;
end

% Interpolate DT to markers
for m=1:1:marknum
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
    tkm(m)=tkm(m)+DT(i,j)*wtmij+DT(i+1,j)*wtmi1j+...
            DT(i,j+1)*wtmij1+DT(i+1,j+1)*wtmi1j1;
    % Interpolate tk2 at 1st timestep
    if(timestep==1)
        tkm(m)=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;
    end
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


% Move markers with 4th order Runge-Kutta
vxm=zeros(4,1);
vym=zeros(4,1);
for m=1:1:marknum
    % Save initial marker coordinates
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
            xm(m)=xA+dt/2*vxm(rk);
            ym(m)=yA+dt/2*vym(rk);
        elseif(rk==3)
            xm(m)=xA+dt*vxm(rk);
            ym(m)=yA+dt*vym(rk);
        end
    end
    % Restore initial coordinates
    xm(m)=xA;
    ym(m)=yA;
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Move markers
    xm(m)=xm(m)+dt*vxmeff;
    ym(m)=ym(m)+dt*vymeff;
end  


figure(1);colormap('Jet');clf
subplot(3,4,1)
pcolor(x,y,log10(ETA)); %caxis([17 23])
shading flat;
axis ij image;
colorbar
title('colormap of log10ETA')
hold on
quiver(xp(3:5:Nx1),yp(3:5:Ny1),vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'w')

subplot(3,4,2)
pcolor(xp,yp,pr)
shading interp;
axis ij image;
colorbar
title('colormap of Pressure')
hold on
quiver(xp(3:5:Nx1),yp(3:5:Ny1),vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(3,4,3)
pcolor(xp,yp,vxp)
shading interp;
axis ij image;
colorbar
title('colormap of vx')
hold on
quiver(xp(3:5:Nx1),yp(3:5:Ny1),vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(3,4,4)
pcolor(xp,yp,vyp)
shading interp;
axis ij image;
colorbar
title('colormap of vy')
hold on
quiver(xp(3:5:Nx1),yp(3:5:Ny1),vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(3,4,5)
pcolor(xp,yp,HS)
shading interp;
axis ij image;
colorbar
title('colormap of HS')

subplot(3,4,6)
pcolor(xp,yp,HA)
shading interp;
axis ij image;
colorbar
title('colormap of HA')

subplot(3,4,7)
pcolor(xp,yp,RHO)
shading interp;
axis ij image;
colorbar
title('colormap of RHO')

subplot(3,4,8)
pcolor(xvx,yvx,KX)
shading interp;
axis ij image;
colorbar
title('colormap of K')

subplot(3,4,9)
pcolor(xp,yp,tk2)
shading interp;
axis ij image;
colorbar
title('colormap of T')

subplot(3,4,10)
pcolor(xp,yp,PHI)
shading interp;
axis ij image;
colorbar
title('colormap of PHI')

subplot(3,4,11)
pcolor(xvx,yvx,gx)
shading interp;
axis ij image;
colorbar
title('colormap of gx')

subplot(3,4,12)
pcolor(xvy,yvy,gy)
shading interp;
axis ij image;
colorbar
title('colormap of gy')





% figure(2);colormap('Jet');clf
% pcolor(x/1000,y/1000,TYP); %caxis([17 23])
% shading interp;
% axis ij image;
% colorbar
% title(['material  time=',num2str(timesum/(365.25*24*3600*1e+6)),' Myr'])
% % hold on
% % quiver(xp(2:4:Nx1-1)/1000,yp(2:4:Ny1-1)/1000,vxp(2:4:Ny1-1,2:4:Nx1-1),vyp(2:4:Ny1-1,2:4:Nx1-1),'w','LineWidth',1.5)


figure(3);colormap('Jet');clf
xm1=xm;ym1=ym;
hold on
for mtype=1:1:5
mm=0;
for m=1:1:marknum
    if(tm(m)==mtype)
        mm=mm+1;
        xm1(mm)=xm(m);
        ym1(mm)=ym(m);
    end
end
if(mtype==1); plot(xm1(1:mm),ym1(1:mm),'. k');end
if(mtype==2); plot(xm1(1:mm),ym1(1:mm),'. b');end
if(mtype==3); plot(xm1(1:mm),ym1(1:mm),'. g');end
if(mtype==4); plot(xm1(1:mm),ym1(1:mm),'. r');end
if(mtype==5); plot(xm1(1:mm),ym1(1:mm),'. y');end
end
axis ij image
title(['material  time=',num2str(timesum/(365.25*24*3600*1e+6)),' Myr'])
% Save figure
% if(timestep==1 || fix(timestep/10)*10==timestep)
%     namefig=['planet2_' num2str(timestep)];
%     print ('-djpeg', '-r150',namefig)
% end

pause(0.1)

% Advance time
timesum=timesum+dt;
end
