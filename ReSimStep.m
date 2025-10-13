clear all

%
Re=100; %10 to 100
scale=0.8;




xe = 1.0; % Domain length in x
ye = 1.0; % Domain length in y
obsWidth=scale*xe;
obsHeight=scale*ye;

nx = 40; % Number of cells in x
ny = 40; % Number of cells in y

dx = xe/nx; % Cell size x
dy = ye/ny; % Cell size y

x = repmat(dx/2:dx:xe-dx/2,ny,1); % Cell centered x location
y = repmat(dy/2:dy:ye-dy/2,nx,1)'; % Cell centered y location

initialVelocity = 1.0; % Top wall u velocity
rho = 1.0; % Fluid density
mu=rho*initialVelocity*xe/Re;
nu = mu/rho; % Fluid kinematic viscosity

u = zeros(ny+2,nx+2); % Face centered u velocity(to the left face)
v = zeros(ny+2,nx+2); % Face centered v velocity(to the bottom face)
p = zeros(ny+2,nx+2); % Cell centered pressure
div = zeros(ny+2,nx+2); % Cell centered divergence

u_old = zeros(ny+2,nx+2);
v_old = zeros(ny+2,nx+2);

dt = min(0.25*dx*dx/nu,4.0*nu/initialVelocity/initialVelocity); % Time step(linear advection diffusion stability condition)

twfin = 1000 * dt; % Stopping criteria for simulation

% Real grid boundary conditions

% Parallel to wall
ut = 0; % Top wall 
ub = 0.0; % Bottom wall
vl = 0.0; % Left wall
vr = 0.0; % Right wall

% Perpendicular to wall(Not implemented - requires a different BC
% function,where outlet velocity is solved for,and slip or no slip
% conditions are used on the other walls)
%(A possible next direction for flow over a step or two-phase flow)
ul = initialVelocity;
ur = 0.0;
vt = 0.0;
vb = 0.0;


% Main loop

t = 0; % Start time for simulation
realTime = 0;
endTime = 8;
iter = endTime/dt;
while t < iter
    [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny); % Set boundary conditions
  
    [u,v] = obstacleVelocity(u,v,obsWidth,obsHeight,nx,ny,dx,dy);
    [u,v] = intermediateVelocity(u,v,u_old,v_old,rho,mu,nx,ny,dx,dy,dt); % Solve for intermediate velocity condition
    [u,v] = obstacleVelocity(u,v,obsWidth,obsHeight,nx,ny,dx,dy);

    [p,a_p] = pressureSolve(u,v,obsWidth,obsHeight,rho,nx,ny,dx,dy,dt,t); % Pressure iterative solver
    [u,v,u_old,v_old] = accel(p,u,v,rho,nx,ny,dx,dy,dt); % Advance time step
    [u,v] = obstacleVelocity(u,v,obsWidth,obsHeight,nx,ny,dx,dy);
    
    t = t + 1;
    realTime = realTime + dt;
end
% Post processing

U = zeros(ny,nx);
V = zeros(ny,nx);
velMag = zeros(ny,nx);

% Interpolate my cell centered U and V velocites
for i=1:nx
    for j=1:ny
        U(j,i) = 0.5 *(u(j+1,i+1) + u(j+1,i+2));
        V(j,i) = 0.5 *(v(j+1,i+1) + v(j+2,i+1));
        P(j,i) = 0.5 *(p(j+1,i+1) + p(j+2,i+1));
        velMag(j,i) = sqrt(U(j,i)^2 + V(j,i)^2);
    end
end

% Calculate divergence and find max
div = calcDiv(div,u,v,nx,ny,dx,dy);
[m,ii,jj] = maxDiv(div,nx,ny);

%%
% Plotting velocity vectors
figure()
subplot(1,2,1)


contourLevels = 50; % increase this number for smoother contours
xx = linspace(0, xe, nx);
yy = linspace(0, ye, ny);



plot=P;



% Create filled contour plot
contourf(xx, yy, plot, contourLevels, 'LineStyle', 'none'); 
box on
xlabel('x')
ylabel('y')
xlim([0 xe])
ylim([0 max(y(:))]) % ensure y-limits are covered
pbaspect([1 1 1])
% Add colour bar
cb = colorbar;
cb.Label.String = 'u_x'; % label for colour bar


% Ensure full coverage of colour range
colormap('jet')
caxis([min(plot(:)) max(plot(:))]) % ensure colour covers full data range

% Draw rectangle for obstacle
hold on
rectangle('Position', [0 0 obsWidth obsHeight], 'FaceColor', [0 0.5 0.5], 'EdgeColor', 'none')
hold off


% Plotting stream particles
subplot(1,2,2)
[verts,averts] = streamslice(x,y,U,V)
box on
streamline([verts averts])
pbaspect([1 1 1])
xlabel('x')
ylabel('y')
rectangle('Position',[0 0 obsWidth obsHeight],'FaceColor',[0 .5 .5])


set(gcf,'Position',[100 100 1000 350])














% Functions

function [u,v] = obstacleVelocity(u,v,obsWidth,obsHeight,nx,ny,dx,dy)

    obsI = obsWidth/dx;
    obsJ = obsHeight/dy;

    % Set velocity inside the obstacle = 0. Ensure U is divergence free.
    for i=1:obsI+2
        for j=1:obsJ+1
            u(j,i) = 0.0;
        end
    end

    for i=1:obsI+1
        for j=1:obsJ+2
            v(j,i) = 0.0;
        end
    end

    % Apply velocity boundary conditions to the obstacle
    for i=3:obsI+1
        u(obsJ+1,i) = -u(obsJ+2,i);
    end

    for j=3:obsJ+1
        v(j,obsI+1) = -v(j,obsI+2);
    end
end

function  [a_e,a_w,a_n,a_s,u,v] = obstaclePressure(a_e,a_w,a_n,a_s,u,v,obsWidth,obsHeight,nx,ny,dx,dy)

    obsI = obsWidth/dx;
    obsJ = obsHeight/dy;

    % Set pressure coefficients along the wall
    a_s(obsJ+2,2:obsI+1) = 0.0;
    a_w(2:obsJ+1,obsI+2) = 0.0;
end

function dt = stabilityCondition(u,v,nx,ny,dx)
    maxVelMag = 0.0;
    for i=1:nx+2
        for j=1:ny+2
            velMag = sqrt(u(j,i)^2 + v(j,i)^2);
            if(velMag > maxVelMag)
                maxVelMag = velMag;
            end
        end
    end

    dt = 0.3*dx/maxVelMag; % CFL stability condition
end

function [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny)
    u_old(:,2) = ul; % Left wall
    u_old(:,nx+2) = u_old(:,nx+1); % Right wall
    u_old(ny+2,:) = 2*ut + u_old(ny+1,:); % Top wall
    u_old(1,:) = 2*ub - u_old(2,:); % Bottom wall

    v_old(ny+2,:) = vt; % Top wall
    v_old(2,:) = vb; % Bottom wall
    v_old(:,1) = 2.0*vl - v_old(:,2); % Left wall
    v_old(:,nx+2) = 2.0*vr - v_old(:,nx+1); % Right wall;

    u = u_old;
    v = v_old;
end

function [u,v] = intermediateVelocity(u,v,u_old,v_old,rho,mu,nx,ny,dx,dy,dt)
    for i=3:nx+1
        for j=2:ny+1
            % Interpolating velocities
            u_e = 0.5*(u_old(j,i) + u_old(j,i+1));
            u_w = 0.5*(u_old(j,i-1) + u_old(j,i));
            u_n = 0.5*(u_old(j,i) + u_old(j+1,i));
            u_s = 0.5*(u_old(j,i) + u_old(j-1,i));
            
            v_n = 0.5*(v_old(j+1,i-1) + v_old(j+1,i));
            v_s = 0.5*(v_old(j,i-1) + v_old(j,i));
            
            % Solving div(rho*u*u) and div(tau) 
            convection = -(rho*u_e*u_e - rho*u_w*u_w)/dx -(rho*v_n*u_n - rho*v_s*u_s)/dy;
            diffusion = mu*(u_old(j,i-1) - 2.0*u_old(j,i) + u_old(j,i+1))/dx/dx + mu*(u_old(j+1,i) - 2.0*u_old(j,i) + u_old(j-1,i))/dy/dy;
            
            % Calculate intermediate u velocity
            u(j,i) = rho*u_old(j,i) + dt*(diffusion + convection);
            u(j,i) = u(j,i)/rho;
        end
    end

    for i=2:nx+1
        for j=3:ny+1
            % Interpolating velocities
            v_e = 0.5*(v_old(j,i) + v_old(j,i+1));
            v_w = 0.5*(v_old(j,i) + v_old(j,i-1));
            v_n = 0.5*(v_old(j,i) + v_old(j+1,i));
            v_s = 0.5*(v_old(j,i) + v_old(j-1,i));

            u_e = 0.5*(u_old(j-1,i+1) + u_old(j,i+1));
            u_w = 0.5*(u_old(j-1,i) + u_old(j,i));

            % Solving div(rho*u*u) and div(tau)
            convection = -(rho*v_e*u_e - rho*v_w*u_w)/dx -(rho*v_n*v_n - rho*v_s*v_s)/dy;
            diffusion = mu*(v_old(j,i-1) - 2*v_old(j,i) + v_old(j,i+1))/dx/dx + mu*(v_old(j+1,i) - 2*v_old(j,i) + v_old(j-1,i))/dy/dy;
            
            % Calculate intermediate v velocity
            v(j,i) = rho*v_old(j,i) + dt*(diffusion + convection);
            v(j,i) = v(j,i)/rho;
        end
    end
end

function [p,a_p] = pressureSolve(u,v,obsWidth,obsHeight,rho,nx,ny,dx,dy,dt,iter)
    [p,a_p] = pressureOptimized(u,v,rho,nx,ny,dx,dy,dt,iter); % only works for square matrix atm
    % [p,a_p] = sor(u,v,obsWidth,obsHeight,rho,nx,ny,dx,dy,dt,iter);
end

function [p,a_p] = pressureOptimized(u,v,rho,nx,ny,dx,dy,dt,iter)
    p = zeros(ny+2,nx+2);

    for i=2:nx+1
        for j=2:ny+1
            rhs(j,i) =((u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy)/dt;
        end
    end
    
    % Block matrices
    e = ones(nx,1);
    Ax = spdiags([e -e e],-1:1,nx,nx)/rho/dx/dx;

    e = ones(ny,1);
    Ay = spdiags([e -e e],-1:1,ny,ny)/rho/dy/dy;

    % Adjusting block matrices to account for boundary conditions
    e = ones(nx,1);
    e(1) = 0;
    e(nx) = 0;
    
    Ax2 = spdiags([0*e -e 0*e ],-1:1,nx,nx)/rho/dx/dx;
    Iy2 = speye(ny);
    Iy2(2:ny,2:ny) = 0;
    Iy3 = speye(ny);
    Iy3(1:ny-1,1:ny-1) = 0;

    e = 2*ones(nx,1);
    e(1) = 1;
    e(nx) = 1;
    Ax3 = spdiags([0*e -e 0*e ],-1:1,nx,nx)/rho/dx/dx;
    Iy4 = speye(ny);
    Iy4(1,1) = 0;
    Iy4(ny,ny) = 0;

    Ix=speye(nx);
    Iy=speye(ny);
    
    % Building Ax = b matrix: A2,x2,b2 linear matrices
    A2 = kron(Iy,Ax) + kron(Ay,Ix) + kron(Iy2,Ax2) + kron(Iy3,Ax2) + kron(Iy4,Ax3);

    b = zeros(ny,nx);
    for i=2:nx+1
        for j=2:ny+1
            b(j-1,i-1) = rhs(j,i);
        end
    end
    b2 = reshape(b,[ny*nx,1]);
    
    [x2] = bicgstab(A2,b2);
    x = reshape(x2,[ny,nx]);

    for i=2:nx+1
        for j=2:ny+1
            p(j,i) = x(j-1,i-1);
        end
    end
    a_p = 1;
end

function [p,a_p] = sor(u,v,obsWidth,obsHeight,rho,nx,ny,dx,dy,dt,iter)
    p = zeros(ny+2,nx+2);
    rhs = zeros(ny+2,nx+2);
    a_p = zeros(ny+2,nx+2);
    a_e = ones(ny+2,nx+2)/rho/dx/dx;
    a_w = ones(ny+2,nx+2)/rho/dx/dx;
    a_n = ones(ny+2,nx+2)/rho/dy/dy;
    a_s = ones(ny+2,nx+2)/rho/dy/dy;

    % a_e(:,nx+1) = 0.0;
    % a_w(:,2) = 0.0;
    a_n(ny+1,:) = 0.0;
    a_s(2,:) = 0.0;
    
    obsI = obsWidth/dx;
    obsJ = obsHeight/dy;
    [a_e,a_w,a_n,a_s,u,v] = obstaclePressure(a_e,a_w,a_n,a_s,u,v,obsWidth,obsHeight,nx,ny,dx,dy);

    a_p = -(a_e + a_w + a_n + a_s);

    maxError = 10;
    TOL = 1e-7;
    counter = 0;
    beta = 1.872;

    while(abs(maxError) > TOL)
        counter = counter + 1;
        maxError = 0;
        
        % Solve for pressure field
        for i=2:nx+1
            for j=2:ny+1

                if(j < obsJ+1 && i < obsI+1)
                    continue;
                end

                rhs(j,i) =(u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy;
                rhs(j,i) = rhs(j,i)/dt -(a_w(j,i)*p(j,i-1) + a_e(j,i)*p(j,i+1) + a_n(j,i)*p(j+1,i) + a_s(j,i)*p(j-1,i));
                p(j,i) = beta*rhs(j,i)/a_p(j,i) +(1-beta)*p(j,i);
            end
        end
        
        % Determine residuals
        for i=2:nx+1
            for j=2:ny+1

                if(j < obsJ+1 && i < obsI+1)
                    continue;
                end

                rhs(j,i) =(u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy;
                error = a_w(j,i)*p(j,i-1) + a_e(j,i)*p(j,i+1) + a_n(j,i)*p(j+1,i) + a_s(j,i)*p(j-1,i) + a_p(j,i)*p(j,i) - rhs(j,i)/dt;
                if(abs(error) > maxError)
                    maxError = abs(error);
                end
            end
        end
    end
    
   
end

function [u,v,u_old,v_old] = accel(p,u,v,rho,nx,ny,dx,dy,dt)
    for i=3:nx+1
        for j=2:ny+1
            % Update u velocity
            u(j,i) = rho*u(j,i) - dt*(p(j,i) - p(j,i-1))/dx;

            u(j,i) = u(j,i) / rho;
        end
    end

    for i=2:nx+1
        for j=3:ny+1
            % Update v velocity
            v(j,i) = rho*v(j,i) - dt*(p(j,i) - p(j-1,i))/dy;

            v(j,i) = v(j,i) / rho;
        end
    end
    u_old = u;
    v_old = v;
end
%% Utility functions

function div = calcDiv(div,u,v,nx,ny,dx,dy)
    for i=2:nx-1
        for j=2:ny-1
            div(j,i) =(u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy;
        end
    end
end

function [max,ii,jj] = maxDiv(div,nx,ny)
    max = 0.0;
    for i=2:nx-1
        for j=2:ny-1
            if(div(j,i) > max)
                max = div(j,i);
                ii = i;
                jj = j;
            end
        end
    end
end