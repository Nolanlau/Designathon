function  MSim(Minf, ax)
%MSim  Run the boundary-layer style MacCormack routine and plot Mach on given axes
%   data = MSim(Minf, ax)
%   Minf : inflow Mach number (scalar)
%   ax   : axes or uiaxes handle where plot will be drawn
%   data : struct with fields x,y,u,v,p,T,delta,Ma (useful for further inspection)

% Physical constants (local)
global mu0 T0 gm Pr R
mu0 = 1.7894e-05;    % Pa s (sea level dynamic viscosity)
T0   = 288.16;       % K
gm   = 1.4;          % gamma
Pr   = 0.71;         % Prandtl
R    = 287;          % J/(kg K)

% Problem geometry and numerics (adopted from your script)
lhori   = 1.0e-5;    % plate length [m]
K       = 0.8;       % Courant number
nx      = 70;
ny      = 70;
maxiter = 1000;

% Inflow primitives (velocity set by Mach)
pinf = 101325;         % Pa
Tinf = 288.16;         % K
% Create inflow primitive state (assumes Primitives accepts u,v,p,T)
inflow = Primitives(0,0,pinf,Tinf);
% set axial velocity from supplied Mach using local sound speed
a_inf = sqrt(gm * R * Tinf);          % a = sqrt(gamma R T)
inflow.u = Minf * a_inf;              % U = M * a

% Reynolds number estimate (calls user function)
try
    [Reinf, ~] = calculateReynoldsNumber(inflow.u, inflow.v, inflow.p, inflow.T, lhori);
catch
    % Fallback: assume ideal gas, compute rho = p/(R T), mu = mu0 (approx)
    rho_inf = pinf / (R * Tinf);
    Reinf = rho_inf * inflow.u * lhori / mu0;
    Reinf = Reinf; % keep scalar
end

% boundary layer scale and vertical domain size
delta = 5 * lhori / sqrt(Reinf);  % as in your script
lvert  = 5 * delta;

% grid (x,y)
[x, y] = meshgrid(linspace(0, lhori, nx), linspace(0, lvert, ny));

% initial condition array of primitives set to inflow
primitives = Primitives(inflow.u * ones(ny, nx), ...
                        inflow.v * ones(ny, nx), ...
                        inflow.p * ones(ny, nx), ...
                        inflow.T * ones(ny, nx));

% Run solver for adiabatic wall (Tw_Tinf = -1 indicates adiabatic in your script)
Tw_Tinf = -1.0;
adiabaticSol = solveMacCormack(primitives, inflow, Tw_Tinf, K, x, y, maxiter);

% Extract fields
[u, v, p, T] = dealPrimitives(adiabaticSol);

% Compute local sound speed array. Prefer getSoundSpeed if present.
if exist('getSoundSpeed', 'file') == 2
    a = getSoundSpeed(T);
else
    a = sqrt(gm * R .* T);
end

% Compute Mach magnitude field
Umag = sqrt(u.^2 + v.^2);
MachField = Umag ./ a;

% Plot to provided axes
contourLevels = 200; % adjust as required
cla(ax);
hold(ax, 'on');

% Use contourf with the axes as first argument
contourf(ax, x./delta, y./delta, MachField, contourLevels, 'LineStyle', 'none');

% Axes and plot configuration
colormap(ax, 'parula');
cb = colorbar(ax);
cb.Label.String = 'Mach Number';
xlabel(ax, 'x ');
ylabel(ax, 'y ', 'Rotation', 0);
box(ax, 'on');
pbaspect(ax, [1 1 1]);
view(ax, 2); % top view equivalent to view(0,90)
set(ax, 'FontName', 'Helvetica');

hold(ax, 'off');



end











































































































%% ========================================================================
%% HELPER FUNCTIONS - All consolidated from separate .m files
%% ========================================================================

function obj = Primitives(uIn,vIn,pIn,TIn)
    % Constructor for Primitives structure
    obj.u = uIn; 
    obj.v = vIn; 
    obj.p = pIn; 
    obj.T = TIn;
end

function r = getDensity(p,T)
    % density of a perfect gas
    global R
    r = p./(R*T);
end

function e = getInternalEnergy(T)
    % internal energy of a calorically perfect gas (constant specific heats)
    global R gm
    cv = R/(gm-1);
    e = cv*T;
end

function Et = getTotalEnergy(r,u,v,T)
    % total energy
    e = getInternalEnergy(T);
    Et = r.*(e + .5*(u.^2 + v.^2));
end

function cv = getCv()
    % ideal gas specific heat at constant volume
    global R gm
    cv = R/(gm-1);
end

function cp = getCp()
    % ideal gas specific heat at constant pressure
    global gm
    cv = getCv();
    cp = gm*cv;
end

function mu = getViscosity(T)
    % Sutherland's Law
    global mu0 T0
    mu = mu0.*(T./T0).^(3/2) .* (T0 + 110)./(T + 110);
end

function lambda = getLambda(T)
    % Stokes's hypothesis
    mu = getViscosity(T);
    lambda = -2/3 * mu;
end

function a = getSoundSpeed(T)
    % ideal gas sound speed
    global gm R
    a = sqrt(gm*R*T);
end

function k = getThermalConductivity(T)
    % Constant Prandtl number
    global Pr
    cp = getCp();
    mu = getViscosity(T);
    k = cp*mu./Pr;
end

function [muOut,lambdaOut,kOut] = getMuLambdaK(T)
    % Get viscosity, lambda, and thermal conductivity efficiently
    global Pr
    muOut = getViscosity(T);
    lambdaOut = -2/3 * muOut;
    cp = getCp();
    kOut = cp*muOut./Pr;
end

function [u,v,p,T] = dealPrimitives(obj)
   % convenience function to deal variables
   u = obj.u; v = obj.v; p = obj.p; T = obj.T;
end

function [ReX,ReY] = calculateReynoldsNumber(u,v,p,T,referenceLength)
    % ReX = [ny,nx double] Reynolds Number for x direction
    % ReY = [ny,nx double] Reynolds Number for y direction
    r = getDensity(p,T);
    mu = getViscosity(T);
    r_mu = r./ mu;
    ReX = referenceLength  .* u .* r_mu;
    ReY = referenceLength  .* v .* r_mu;
end

function dt = calculateTimeStep(u,v,p,T,dx,dy,K)
    global gm Pr
    mu = getViscosity(T);
    r = getDensity(p,T);
    a = getSoundSpeed(T);
    vp = max(4/3*mu, gm*mu./Pr)./r; % Anderson has typos
    dtCFL = 1./( abs(u)./dx + abs(v)./dy + a.*sqrt(1/dx^2 + 1/dy^2) + 2*vp*(1/dx^2 + 1/dy^2));
    dt = K *min(dtCFL(:));
end

function [primitives,massDiffCheck,converged,i,rumTime] = ...
          solveMacCormack(primitives,inflow,Tw_Tinf,K,x,y,maxiter)
% Solve Navier-Stokes equations for supersonic flow over a flat plate using MacCormack method.
%
% INPUTS
% primitives = [Primitives] Initial domain primitives
% inflow     = [Primitives] Primitives at the inflow boundary
% Tw_Tinf    = [double] Wall temperature specification
%                  if Tw_Tinf > 0
%                        Tw/Tinf = Tw_Tinf
%				 else
%                        Adiabatic wall
%                  end
% K          = [double] Courant number (0.5 <= K <= 0.8 recommended)
% x          = [ny,nx double] Grid point x locations (must be uniform spacing)
% x          = [ny,nx double] Grid point y locations (must be uniform spacing)
% maxiter    = [double] maximum number of time steps
% 
% OUTPUTS
% primitives    = [Primitives] Final domain primitives
% massDiffCheck = [double] Percent difference between inflow and outflow mass
% converged     = [logical] True if convergence criteria satisfied
% i             = [double] number of iterations
% rumTime       = [double] solution runtime in seconds
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021						   
tic

% mesh data
[ny,nx] = size(x);
dx = x(1,2)-x(1,1);
dy = y(2,1)-y(1,1);
inX = 2:nx-1; % interior x
inY = 2:ny-1; % interior y

% set boundary conditions
primitives=updateBoundaryConditions(primitives,inflow,Tw_Tinf);

% time march
i = 0;
converged = false;
while ~converged && i < maxiter
    i = i + 1;
    
    % time step size
    [u,v,p,T] = dealPrimitives(primitives);
    dt = calculateTimeStep(u,v,p,T,dx,dy,K);
    
    % solution vector
    U = calculateU(primitives);
    
    % flux vectors
    E = calculateE(primitives,dx,dy,'backward');
    F = calculateF(primitives,dx,dy,'backward');
    
    % forward finite difference predictor
    dUdt_predictor = -ddxf(E,dx) - ddyf(F,dy);
    
    % Predictor step and corrector vector calculations
    U2 = U;
    U2(inY,inX,:) = U(inY,inX,:) + dt*dUdt_predictor(inY,inX,:); % interior points only
    primitives2=decodeSolutionVector(U2);
    E2 = calculateE(primitives2,dx,dy,'forward');
    F2 = calculateF(primitives2,dx,dy,'forward');
    
    % backward finite difference corrector
    dUdt_corrector = -ddxb(E2,dx) - ddyb(F2,dy);
    
    % MacCormack solution step
    dUdt = 0.5*(dUdt_predictor + dUdt_corrector);
    U(inY,inX,:) = U(inY,inX,:) + dt*dUdt(inY,inX,:);% interior points only
    primitives=decodeSolutionVector(U);
    primitives=updateBoundaryConditions(primitives,inflow,Tw_Tinf);
    
    % check density convergence
    [u,v,p,T] = dealPrimitives(primitives);
    rCurrent = getDensity(p,T);
    if i > 1
        deltaR = max(max(abs(rCurrent - rLast)));
        if deltaR < 1e-8
            converged = true;
        end
        fprintf(1,'Iteration:%5d | delta rho: %8e\n',i,deltaR);
    end
    rLast = rCurrent;
    
end
runTime = toc;

% Mass Flow Check
[u,v,p,T] = dealPrimitives(primitives);
r = getDensity(p,T);
massIn = trapz(y(:,1),u(:,1).*r(:,1));
massOut = trapz(y(:,end),u(:,end).*r(:,end));
massDiffCheck = 100*abs(massIn-massOut)./massIn;
fprintf(1,'Mass inflow matches mass outflow within %.3f%%.\n',massDiffCheck);
fprintf(1,'Runtime: %.2f seconds.\n',runTime);
end

function E = calculateE(primitives,dx,dy,direction)
% Calculate E flux array from primitive variables
%
% INPUTS
% primitives = [Primitives] Domain primitives
% dx         = [double] x grid spacing
% dy         = [double] y grid spacing
% direction  = [char] Finite difference direction = 'forward' or 'backward'
%
% OUTPUTS
% E = [ny,ny,4 double] Solution array (See Anderson Eq. 10.4c)
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.
%
% From Anderson page 454:
% To maintain second-order accuracy, the x-derivative terms appearing in E
% are differenced in the opposite direction to that used for dE/dx, while
% the y-derivative terms are approximated using central differences.
% Likewise, the y-derivative terms appearing in F are differenced in the
% opposite direction to that used for dF/dy, while the x-derivative terms
% in F are central differenced. 

% Anthony Ricciardi
% July 2021

% primitives
[u,v,p,T] = dealPrimitives(primitives);
[mu,lambda,k] = getMuLambdaK(T);
r = getDensity(p,T);
Et = getTotalEnergy(r,u,v,T);

% gradients
switch direction
    case 'forward'
        dudx = ddxf(u,dx);
        dvdx = ddxf(v,dx);
        dTdx = ddxf(T,dx);
    case 'backward'
        dudx = ddxb(u,dx);
        dvdx = ddxb(v,dx);
        dTdx = ddxb(T,dx);
    otherwise
        error('direction not allowed')
end
dudy = ddyc(u,dy);
dvdy = ddyc(v,dy);
% dTdy = ddyc(T,dy);

% stresses and heat fluxes
txx = lambda.*(dudx + dvdy)+ 2*mu.*(dudx);
% tyy = lambda.*(dudx + dvdy)+ 2*mu.*(dvdy);
txy = mu.*(dudy + dvdx);
qx = -k.*dTdx;
% qy = -k.*dTdy;

E = zeros(size(r,1),size(r,2),4);
E(:,:,1) = r.*u;
E(:,:,2) = r.*u.^2 + p - txx;
E(:,:,3) = r.*u.*v - txy;
E(:,:,4) = (Et+p).*u- u.*txx - v.*txy + qx;
end

function F = calculateF(primitives,dx,dy,direction)
% Calculate F flux array from primitive variables
%
% INPUTS
% primitives = [Primitives] Domain primitives
% dx         = [double] x grid spacing
% dy         = [double] y grid spacing
% direction  = [char] Finite difference direction = 'forward' or 'backward'
%
% OUTPUTS
% F = [ny,ny,4 double] Solution array (See Anderson Eq. 10.4d [typo in term 3])
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.
%
% From Anderson page 454:
% To maintain second-order accuracy, the x-derivative terms appearing in E
% are differenced in the opposite direction to that used for dE/dx, while
% the y-derivative terms are approximated using central differences.
% Likewise, the y-derivative terms appearing in F are differenced in the
% opposite direction to that used for dF/dy, while the x-derivative terms
% in F are central differenced. 

% Anthony Ricciardi
% July 2021

% primitives
[u,v,p,T] = dealPrimitives(primitives);
[mu,lambda,k] = getMuLambdaK(T);

r = getDensity(p,T);
Et = getTotalEnergy(r,u,v,T);

% gradients
switch direction
    case 'forward'
        dudy = ddyf(u,dy);
        dvdy = ddyf(v,dy);
        dTdy = ddyf(T,dy);
    case 'backward'
        dudy = ddyf(u,dy);
        dvdy = ddyf(v,dy);
        dTdy = ddyf(T,dy);
    otherwise
        error('direction not allowed')
end
dudx = ddxc(u,dx);
dvdx = ddxc(v,dx);
% dTdx = ddxc(T,dx);

% stresses and heat fluxes
% txx = lambda.*(dudx + dvdy)+ 2*mu.*(dudx);
tyy = lambda.*(dudx + dvdy)+ 2*mu.*(dvdy);
txy = mu.*(dudy + dvdx);
% qx = -k.*dTdx;
qy = -k.*dTdy;

F = zeros(size(r,1),size(r,2),4);
F(:,:,1) = r.*v;
F(:,:,2) = r.*u.*v - txy;
F(:,:,3) = r.*v.^2 + p - tyy;
F(:,:,4) = (Et+p).*v- u.*txy - v.*tyy + qy;
end

function U=calculateU(primitives)
% Calculate solution array from primitive variables
%
% INPUTS
% primitives = [Primitives] Domain primitives
%
% OUTPUTS
% U = [ny,ny,4 double] Solution array (See Anderson Eq. 10.4b)
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[u,v,p,T] = dealPrimitives(primitives);
r = getDensity(p,T);
Et = getTotalEnergy(r,u,v,T);

U = zeros(size(u,1),size(u,2),4);
U(:,:,1) = r;
U(:,:,2) = U(:,:,1).*u; % r*u
U(:,:,3) = U(:,:,1).*v; % r*v
U(:,:,4) = Et;
end

function dfdx = ddxb(f,dx)
% Calculates backward-difference derivative wrt x on even-spaced 2D grid
% 
% INPUTS
% f  = [ny,nx,: double] array with function values
% dx = [double] grid spacing
%
% OUTPUTS
% dfdx = [ny,nx,: double] array with function derivatives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[ny,nx,n3] = size(f);
inX = 2:nx;
dfdx = zeros(ny,nx,n3);

% backward difference interior points
dfdx(:,inX,:) = (f(:,inX,:) - f(:,inX-1,:))./dx;

% Set boundary point derivatives equal to adjacent points. Boundary points 
% are essentially unused in the solution, but they are provided for array 
% size consistency.
dfdx(:,1,:) = dfdx(:,2,:);
end

function dfdx = ddxc(f,dx)
% Calculate central-difference derivative wrt x on even-spaced 2D grid 
% 
% INPUTS
% f  = [ny,nx double] array with function values
% dx = [double] grid spacing
%
% OUTPUTS
% dfdx = [ny,nx double] array with function derivatives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[ny,nx] = size(f);
inX = 2:nx-1; % interior x

% initialize
dfdx = zeros(ny,nx);

% central difference interior points  [O(dx^2)]
dfdx(:,inX) = (f(:,inX+1) - f(:,inX-1))./(2*dx);

% % one-sided difference boundary points  [O(dx)]
dfdx(:,1) = (f(:,2) - f(:,1))./(dx);
dfdx(:,nx) = (f(:,nx) - f(:,nx-1))./(dx);

% one-sided difference boundary points [O(dx^2)]
% dfdx(:,1) = (-f(:,3) + 4*f(:,2) - 3*f(:,1))./(2*dx);
% dfdx(:,nx) = (3*f(:,nx) - 4*f(:,nx-1) + f(:,nx-2))./(2*dx);
end

function dfdx = ddxf(f,dx)
% Calculates forward-difference derivative wrt x on even-spaced 2D grid
% 
% INPUTS
% f  = [ny,nx,: double] array with function values
% dx = [double] grid spacing
%
% OUTPUTS
% dfdx = [ny,nx,: double] array with function derivatives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[ny,nx,n3] = size(f);
inX = 1:nx-1;
dfdx = zeros(ny,nx,n3);

% forward difference interior points
dfdx(:,inX,:) = (f(:,inX+1,:) - f(:,inX,:))./dx;

% Set boundary point derivatives equal to adjacent points. Boundary points 
% are essentially unused in the solution, but they are provided for array 
% size consistency.
dfdx(:,nx,:) = dfdx(:,nx-1,:);
end

function dfdy = ddyb(f,dy)
% Calculates backward-difference derivative wrt y on even-spaced 2D grid
% 
% INPUTS
% f  = [ny,nx,: double] array with function values
% dy = [double] grid spacing
%
% OUTPUTS
% dfdy = [ny,nx,: double] array with function derivatives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[ny,nx,n3] = size(f);
inY = 2:ny; 
dfdy = zeros(ny,nx,n3);

% backward difference interior points
dfdy(inY,:,:) = (f(inY,:,:) - f(inY-1,:,:))./dy;

% Set boundary point derivatives equal to adjacent points. Boundary points 
% are essentially unused in the solution, but they are provided for array 
% size consistency.
dfdy(1,:,:) = dfdy(2,:,:);
end

function dfdy = ddyc(f,dy)
% Calculates central-difference derivative wrt y on even-spaced 2D grid 
% 
% INPUTS
% f  = [ny,nx double] array with function values
% dy = [double] grid spacing
%
% OUTPUTS
% dfdy = [ny,nx double] array with function derivatives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[ny,nx] = size(f);
inY = 2:ny-1; % interior y

% initialize
dfdy = zeros(ny,nx);

% central difference interior points  [O(dy^2)]
dfdy(inY,:) = (f(inY+1,:) - f(inY-1,:))./(2*dy);

% % one-sided difference boundary points [O(dy)]
dfdy(1,:) = (f(2,:) - f(1,:))./(dy);
dfdy(ny,:) = (f(ny,:) - f(ny-1,:))./(dy);

% one-sided difference boundary points [O(dy^2)]
% dfdy(1,:) = (-f(3,:) + 4*f(2,:) - 3*f(1,:))./(2*dy);
% dfdy(ny,:) = (3*f(ny,:) - 4*f(ny-1,:) + f(ny-2,:))./(2*dy);
end

function dfdy = ddyf(f,dy)
% Calculates forward-difference derivative wrt y on even-spaced 2D grid
% 
% INPUTS
% f  = [ny,nx,: double] array with function values
% dy = [double] grid spacing
%
% OUTPUTS
% dfdy = [ny,nx,: double] array with function derivatives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[ny,nx,n3] = size(f);
inY = 1:ny-1;
dfdy = zeros(ny,nx,n3);

% forward difference interior points
dfdy(inY,:,:) = (f(inY+1,:,:) - f(inY,:,:))./dy;

% Set boundary point derivatives equal to adjacent points. Boundary points 
% are essentially unused in the solution, but they are provided for array 
% size consistency.
dfdy(ny,:,:) = dfdy(ny-1,:,:);
end

function primitives=decodeSolutionVector(U)
% Calculate primitive variables from solution vector
%
% INPUTS
% U = [ny,ny,4 double] Solution array (See Anderson Eq. 10.4b)
%
% OUTPUTS
% primitives = [Primitives] Domain primitives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

global R gm
r = U(:,:,1);
u = U(:,:,2)./r;
v = U(:,:,3)./r;
Et = U(:,:,4);
e = Et./r - .5*(u.^2 + v.^2);
cv = R/(gm-1);
T = e./cv;
p = r.*R.*T;
primitives = Primitives(u,v,p,T);
end

function primitivesOut=updateBoundaryConditions(primitivesIn,inflow,Tw_Tinf)
% Update solution boundary conditions
%
% INPUTS
% primitivesIn = [Primitives] Domain primitives
% inflow       = [Primitives] Primitives at the inflow boundary
% Tw_Tinf      = [double] Wall temperature specification
%                    if Tw_Tinf > 0
%                          Tw/Tinf = Tw_Tinf
%					 else
%                          Adiabatic wall
%                    end
%
% OUTPUTS
% primitivesOut = [Primitives] Domain primitives with boundary values updated

% Anthony Ricciardi
% July 2021						   

[u,v,p,T] = dealPrimitives(primitivesIn);
[Vinf,~,pinf,Tinf] = dealPrimitives(inflow);

% inflow bounary, x(:,1) == 0
u(:,1) = Vinf;
v(:,1) = 0;
p(:,1) = pinf;
T(:,1) = Tinf;

% upper boundary, y(end,:) == lvert
u(end,:) = Vinf;
v(end,:) = 0;
p(end,:) = pinf;
T(end,:) = Tinf;

% outflow boundary, x(:,end) == lhori
u(:,end) = 2*u(:,end-1) - u(:,end-2);
v(:,end) = 2*v(:,end-1) - v(:,end-2);
p(:,end) = 2*p(:,end-1) - p(:,end-2);
T(:,end) = 2*T(:,end-1) - T(:,end-2);

% plate boundary, y(1,:) == 0
u(1,:) = 0;
v(1,:) = 0;
p(1,:) = 2*p(2,:) - p(3,:);

if Tw_Tinf > 0
    % constant temperature wall
    T(1,:) = Tinf * Tw_Tinf;
else 
    % adiabatic wall
    % 0 = dfdy(1,:) = (-f(3,:) + 4*f(2,:) - 3*f(1,:))./(2*dy);
    % 0 = -f(3,:) + 4*f(2,:) - 3*f(1,:) 
    % 3*f(1,:) = -f(3,:) + 4*f(2,:)
    % f(1,:) = (4*f(2,:)- f(3,:))./3
    
%     T(1,:) = (4*T(2,:)- T(3,:))./3; % [O(dy^2)]
    T(1,:) = T(2,:);  % [O(dy^2)] - Seems to match Anderson more closely
end

% leading edge, x(1,1) == y(1,1) == 0, supersedes others
u(1,1) = 0;
v(1,1) = 0;
p(1,1) = pinf;
T(1,1) = Tinf;

% Export updated primitives
primitivesOut = Primitives(u,v,p,T);
end
