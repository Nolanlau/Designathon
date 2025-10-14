function out = ReSimStep(mode, varargin)

switch lower(mode)
    case 'compute'
        % ---------- INPUTS ----------
        Re        = varargin{1};       % Reynolds number
        stepW_rel = varargin{2};       % step width   (0..1 of domain)
        stepH_rel = varargin{3};       % step height  (0..1 of domain)

        % ---------- DOMAIN / GRID ----------
        xe = 1.0; ye = 1.0;            % domain size
        nx = 40; ny = 40;              % cells
        dx = xe/nx; dy = ye/ny;

        % obstacle (backward-facing step) dimensions
        obsWidth  = max(0,min(xe, stepW_rel*xe));
        obsHeight = max(0,min(ye, stepH_rel*ye));

        % cell-centered coordinates (for plots)
        xcc = repmat(dx/2:dx:xe-dx/2, ny, 1);
        ycc = repmat(dy/2:dy:ye-dy/2, nx, 1)';

        % ---------- PHYSICALS ----------
        initialVelocity = 1.0; % top-wall speed (sets scale)
        rho = 1.0;
        mu  = rho*initialVelocity*xe/Re;
        nu  = mu/rho;

        % ---------- FIELDS (staggered) ----------
        u = zeros(ny+2,nx+2); v = zeros(ny+2,nx+2); p = zeros(ny+2,nx+2);
        u_old = u; v_old = v;

        % ---------- TIME STEP ----------
        dt = min(0.25*dx*dx/nu, 4.0*nu/(initialVelocity^2));
        endTime = 4;
        iters = max(1, ceil(endTime/dt));
        t = 0;

        % ---------- WALL BC VALUES ----------
        ut = 0; ub = 0; vl = 0; vr = 0;
        ul = initialVelocity; ur = 0; vt = 0; vb = 0;

        % ---------- MAIN LOOP ----------
        for k = 1:iters
            [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny);
            [u,v] = obstacleVelocity(u,v,obsWidth,obsHeight,nx,ny,dx,dy);
            [u,v] = intermediateVelocity(u,v,u_old,v_old,rho,mu,nx,ny,dx,dy,dt);
            [u,v] = obstacleVelocity(u,v,obsWidth,obsHeight,nx,ny,dx,dy);
            [p,~] = pressureSolve(u,v,obsWidth,obsHeight,rho,nx,ny,dx,dy,dt,k);
            [u,v,u_old,v_old] = accel(p,u,v,rho,nx,ny,dx,dy,dt);
            [u,v] = obstacleVelocity(u,v,obsWidth,obsHeight,nx,ny,dx,dy);
            t = t + dt; 
        end

        % ---------- INTERPOLATE TO CELL CENTERS ----------
        U = zeros(ny,nx); V = zeros(ny,nx); P = zeros(ny,nx); velMag = zeros(ny,nx);
        for i = 1:nx
            for j = 1:ny
                U(j,i) = 0.5*(u(j+1,i+1) + u(j+1,i+2));
                V(j,i) = 0.5*(v(j+1,i+1) + v(j+2,i+1));
                P(j,i) = 0.5*(p(j+1,i+1) + p(j+2,i+1));
                velMag(j,i) = hypot(U(j,i), V(j,i));
            end
        end

        % divergence (diagnostic)
        div = calcDiv(zeros(ny+2,nx+2), u, v, nx, ny, dx, dy); 

        % ---------- OUTPUT ----------
        data.x = xcc; data.y = ycc;
        data.U = U; data.V = V; data.P = P; data.Vel = velMag;
        data.xe = xe; data.ye = ye; data.nx = nx; data.ny = ny;
        data.dx = dx; data.dy = dy;
        data.obsWidth = obsWidth; data.obsHeight = obsHeight;
        data.stepW_rel = stepW_rel; data.stepH_rel = stepH_rel;
        data.Re = Re; data.mu = mu; data.nu = nu;
        out = data;

    case 'plot'
        % ---------- INPUTS ----------
        data     = varargin{1};
        fieldStr = lower(string(varargin{2}));
        axHeat   = varargin{3};  % UIAxes_3
        axStream = varargin{4};  % UIAxes2_2

        % ---------- SELECT FIELD ----------
        switch fieldStr
            case {'pressure','p'}
                Z = data.P; cbarLabel = 'Pressure';
            case {'horizontal velocity','u','u_x','ux'}
                Z = data.U; cbarLabel = 'Horizontal Velocity';
            case {'vertical velocity','v','u_y','uy'}
                Z = data.V; cbarLabel = 'Vertical Velocity';
            case {'absolute velocity','total velocity','|u|','magnitude'}
                Z = data.Vel; cbarLabel = 'Absolute Velocity';
        end

        % ---------- HEAT MAP ----------
        cla(axHeat); hold(axHeat, 'on');
        contourLevels = 50;
        xx = linspace(0, data.xe, data.nx);
        yy = linspace(0, data.ye, data.ny);
        contourf(axHeat, xx, yy, Z, contourLevels, 'LineStyle', 'none');
        rectangle(axHeat, 'Position', [0 0 data.obsWidth data.obsHeight], ...
            'FaceColor', [0 0.5 0.5], 'EdgeColor', 'none');
        box(axHeat, 'on');
        xlabel(axHeat, 'x'); ylabel(axHeat, 'y','Rotation',0);
        xlim(axHeat, [0 data.xe]); ylim(axHeat, [0 data.ye]);
        pbaspect(axHeat, [1 1 1]);
        colormap(axHeat, 'parula');
        caxis(axHeat, [min(Z(:)) max(Z(:))]);
        cb = colorbar(axHeat);
        cb.Label.String = cbarLabel;
        hold(axHeat, 'off');

        % ---------- STREAMLINES ----------
        cla(axStream); hold(axStream, 'on');
        % Use streamslice for responsive UI (no blocking seed selection)
        streamslice(axStream, data.x, data.y, data.U, data.V);
        rectangle(axStream, 'Position', [0 0 data.obsWidth data.obsHeight], ...
            'FaceColor', [0 0.5 0.5], 'EdgeColor', 'none');
        box(axStream, 'on');
        xlabel(axStream, 'x'); ylabel(axStream, 'y','Rotation',0);
        xlim(axStream, [0 data.xe]); ylim(axStream, [0 data.ye]);
        pbaspect(axStream, [1 1 1]);
        hold(axStream, 'off');

        out = []; % no data return in plot mode

    otherwise
        error('Mode must be "compute" or "plot".');
end

% ====================== SUBFUNCTIONS ======================

function [u,v] = obstacleVelocity(u,v,obsWidth,obsHeight,nx,ny,dx,dy)
    obsI = max(0, min(nx, floor(obsWidth/dx)));
    obsJ = max(0, min(ny, floor(obsHeight/dy)));

    % zero velocity in obstacle region (staggered correction)
    for i = 1:obsI+2
        for j = 1:obsJ+1
            u(j,i) = 0.0;
        end
    end
    for i = 1:obsI+1
        for j = 1:obsJ+2
            v(j,i) = 0.0;
        end
    end

    % mirror conditions on the obstacle "walls"
    for i = 3:obsI+1
        u(obsJ+1,i) = -u(obsJ+2,i);
    end
    for j = 3:obsJ+1
        v(j,obsI+1) = -v(j,obsI+2);
    end
end

function  [a_e,a_w,a_n,a_s,u,v] = obstaclePressure(a_e,a_w,a_n,a_s,u,v,obsWidth,obsHeight,nx,ny,dx,dy) %#ok<INUSD>
    obsI = max(0, min(nx, floor(obsWidth/dx)));
    obsJ = max(0, min(ny, floor(obsHeight/dy)));
    a_s(obsJ+2,2:obsI+1) = 0.0;
    a_w(2:obsJ+1,obsI+2) = 0.0;
end

function [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny)
    u_old(:,2)     = ul;                  % left inflow (no-penetration via v BC)
    u_old(:,nx+2)  = u_old(:,nx+1);       % right Neumann
    u_old(ny+2,:)  = 2*ut + u_old(ny+1,:);% top
    u_old(1,:)     = 2*ub - u_old(2,:);   % bottom

    v_old(ny+2,:)  = vt;                  % top
    v_old(2,:)     = vb;                  % bottom
    v_old(:,1)     = 2.0*vl - v_old(:,2); % left
    v_old(:,nx+2)  = 2.0*vr - v_old(:,nx+1); % right

    u = u_old; v = v_old;
end

function [u,v] = intermediateVelocity(u,v,u_old,v_old,rho,mu,nx,ny,dx,dy,dt)
    % u-momentum (u is face-centered on vertical cell faces)
    for i = 3:nx+1
        for j = 2:ny+1
            u_e = 0.5*(u_old(j,i) + u_old(j,i+1));
            u_w = 0.5*(u_old(j,i-1) + u_old(j,i));
            u_n = 0.5*(u_old(j,i) + u_old(j+1,i));
            u_s = 0.5*(u_old(j,i) + u_old(j-1,i));
            v_n = 0.5*(v_old(j+1,i-1) + v_old(j+1,i));
            v_s = 0.5*(v_old(j,i-1) + v_old(j,i));
            convection = -(rho*u_e*u_e - rho*u_w*u_w)/dx -(rho*v_n*u_n - rho*v_s*u_s)/dy;
            diffusion  =  mu*(u_old(j,i-1) - 2.0*u_old(j,i) + u_old(j,i+1))/dx/dx ...
                        +  mu*(u_old(j+1,i) - 2.0*u_old(j,i) + u_old(j-1,i))/dy/dy;
            u(j,i) = (rho*u_old(j,i) + dt*(diffusion + convection))/rho;
        end
    end

    % v-momentum (v is face-centered on horizontal cell faces)
    for i = 2:nx+1
        for j = 3:ny+1
            v_e = 0.5*(v_old(j,i) + v_old(j,i+1));
            v_w = 0.5*(v_old(j,i) + v_old(j,i-1));
            v_n = 0.5*(v_old(j,i) + v_old(j+1,i));
            v_s = 0.5*(v_old(j,i) + v_old(j-1,i));
            u_e = 0.5*(u_old(j-1,i+1) + u_old(j,i+1));
            u_w = 0.5*(u_old(j-1,i) + u_old(j,i));
            convection = -(rho*v_e*u_e - rho*v_w*u_w)/dx -(rho*v_n*v_n - rho*v_s*v_s)/dy;
            diffusion  =  mu*(v_old(j,i-1) - 2*v_old(j,i) + v_old(j,i+1))/dx/dx ...
                        +  mu*(v_old(j+1,i) - 2*v_old(j,i) + v_old(j-1,i))/dy/dy; 
            diffusion  =  mu*(v_old(j,i-1) - 2*v_old(j,i) + v_old(j,i+1))/dx/dx ...
                        +  mu*(v_old(j+1,i) - 2*v_old(j,i) + v_old(j-1,i))/dy/dy;
            v(j,i) = (rho*v_old(j,i) + dt*(diffusion + convection))/rho;
        end
    end
end

function [p,a_p] = pressureSolve(u,v,obsWidth,obsHeight,rho,nx,ny,dx,dy,dt,iter) %#ok<INUSD>
    [p,a_p] = pressureOptimized(u,v,rho,nx,ny,dx,dy,dt,iter);
    % alt: [p,a_p] = sor(u,v,obsWidth,obsHeight,rho,nx,ny,dx,dy,dt,iter);
end

function [p,a_p] = pressureOptimized(u,v,rho,nx,ny,dx,dy,dt,iter) %#ok<INUSD>
    p = zeros(ny+2,nx+2);
    rhs = zeros(ny+2,nx+2);
    for i=2:nx+1
        for j=2:ny+1
            rhs(j,i) = ((u(j,i+1)-u(j,i))/dx + (v(j+1,i)-v(j,i))/dy)/dt;
        end
    end

    % 2D Poisson via Kronecker (Dirichlet-like handling on outer bands)
    e  = ones(nx,1);
    Ax = spdiags([e -2*e e], -1:1, nx, nx)/(rho*dx*dx);
    e  = ones(ny,1);
    Ay = spdiags([e -2*e e], -1:1, ny, ny)/(rho*dy*dy);

    Ix = speye(nx); Iy = speye(ny);
    A2 = kron(Iy,Ax) + kron(Ay,Ix);

    b = zeros(ny,nx);
    for i=2:nx+1
        for j=2:ny+1
            b(j-1,i-1) = rhs(j,i);
        end
    end
    b2 = reshape(b, [ny*nx, 1]);

    x2 = bicgstab(A2, b2);
    x  = reshape(x2,[ny,nx]);

    for i=2:nx+1
        for j=2:ny+1
            p(j,i) = x(j-1,i-1);
        end
    end
    a_p = 1;
end

function [p,a_p] = sor(u,v,obsWidth,obsHeight,rho,nx,ny,dx,dy,dt,iter) %#ok<INUSD>
    p = zeros(ny+2,nx+2); rhs = zeros(ny+2,nx+2);
    a_e = ones(ny+2,nx+2)/(rho*dx*dx);
    a_w = ones(ny+2,nx+2)/(rho*dx*dx);
    a_n = ones(ny+2,nx+2)/(rho*dy*dy);
    a_s = ones(ny+2,nx+2)/(rho*dy*dy);

    a_n(ny+1,:) = 0.0; a_s(2,:) = 0.0;
    [a_e,a_w,a_n,a_s,~,~] = obstaclePressure(a_e,a_w,a_n,a_s,u,v,obsWidth,obsHeight,nx,ny,dx,dy);

    a_p = -(a_e + a_w + a_n + a_s);
    TOL = 1e-7; beta = 1.872; maxError = 10;

    while abs(maxError) > TOL
        maxError = 0;
        for i=2:nx+1
            for j=2:ny+1
                rhs(j,i) = ( (u(j,i+1)-u(j,i))/dx + (v(j+1,i)-v(j,i))/dy )/dt;
                res = rhs(j,i) - (a_w(j,i)*p(j,i-1) + a_e(j,i)*p(j,i+1) + ...
                                  a_n(j,i)*p(j+1,i) + a_s(j,i)*p(j-1,i) + a_p(j,i)*p(j,i));
                p(j,i) = p(j,i) + beta*res/a_p(j,i);
                maxError = max(maxError, abs(res));
            end
        end
    end
end

function [u,v,u_old,v_old] = accel(p,u,v,rho,nx,ny,dx,dy,dt)
    for i=3:nx+1
        for j=2:ny+1
            u(j,i) = (rho*u(j,i) - dt*(p(j,i) - p(j,i-1))/dx)/rho;
        end
    end
    for i=2:nx+1
        for j=3:ny+1
            v(j,i) = (rho*v(j,i) - dt*(p(j,i) - p(j-1,i))/dy)/rho;
        end
    end
    u_old = u; v_old = v;
end

function div = calcDiv(div,u,v,nx,ny,dx,dy)
    for i=2:nx-1
        for j=2:ny-1
            div(j,i) = (u(j,i+1)-u(j,i))/dx + (v(j+1,i)-v(j,i))/dy;
        end
    end
end
end
