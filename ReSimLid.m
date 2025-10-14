function out = ReSimLid(mode, arg2, arg3, arg4, arg5, varargin)
% NewReSim_external  Two mode helper for App Designer
% Usage:
%   data = NewReSim_external('compute', Re, nx, ny, xe, ye, Utop, endTime)
%       Runs the solver once and returns a struct 'data' containing grids
%       and fields required for plotting.
%
%   NewReSim_external('plot', data, plotChoice, axContour, axStream)
%       Plots from the supplied data struct without recomputing.
%
% Notes:
%   Viscosity chosen as mu = rho * Utop * L / Re.
%   The function expects helper routines on path:
%     boundaryConditions, intermediateVelocity, pressureSolve, accel,
%     calcDiv, maxDiv

if nargin < 2
    error('NewReSim_external:Insufficient arguments.');
end

switch lower(string(mode))
    case 'compute'
        Re = arg2;
        % Optional parameters with conservative defaults
        if nargin >= 3 && ~isempty(arg3), nx = arg3; else nx = 32; end
        if nargin >= 4 && ~isempty(arg4), ny = arg4; else ny = 32; end
        if nargin >= 5 && ~isempty(arg5), xe = arg5; else xe = 1.0; end
        if numel(varargin) >= 1 && ~isempty(varargin{1}), ye = varargin{1}; else ye = 1.0; end
        if numel(varargin) >= 2 && ~isempty(varargin{2}), Utop = varargin{2}; else Utop = 1.0; end
        if numel(varargin) >= 3 && ~isempty(varargin{3}), endTime = varargin{3}; else endTime = 8.0; end

        % Derived quantities
        dx = xe / nx; dy = ye / ny;
        xc = linspace(dx/2, xe-dx/2, nx);
        yc = linspace(dy/2, ye-dy/2, ny);
        [Xc, Yc] = meshgrid(xc, yc);

        % Allocate fields with ghost cells
        u = zeros(ny+2, nx+2);
        v = zeros(ny+2, nx+2);
        p = zeros(ny+2, nx+2);
        div = zeros(ny+2, nx+2);
        u_old = zeros(ny+2, nx+2);
        v_old = zeros(ny+2, nx+2);

        rho = 1.0;
        mu = rho * Utop * xe / Re; % mu = rho U L / Re
        nu = mu / rho;

        % Time step heuristic
        dt = min(0.25 * dx * dx / nu, 4.0 * nu / (Utop^2));
        iter = ceil(endTime / dt);

        % Boundary values
        ut = Utop; ub = 0.0; vl = 0.0; vr = 0.0;
        ul = 0.0; ur = 0.0; vt = 0.0; vb = 0.0;

        % Time integration loop
        for it = 1:iter
            % The calls below must exist on the MATLAB path with the shown signatures
            [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny);
            [u,v] = intermediateVelocity(u,v,u_old,v_old,rho,mu,nx,ny,dx,dy,dt);
            [p,a_p] = pressureSolve(u,v,rho,nx,ny,dx,dy,dt,it);
            [u,v,u_old,v_old] = accel(p,u,v,rho,nx,ny,dx,dy,dt);
        end

        % Interpolate to cell centres
        U = zeros(ny, nx); V = zeros(ny, nx); P = zeros(ny, nx); velMag = zeros(ny, nx);
        for i = 1:nx
            for j = 1:ny
                U(j,i) = 0.5 * ( u(j+1,i+1) + u(j+1,i+2) );
                V(j,i) = 0.5 * ( v(j+1,i+1) + v(j+2,i+1) );
                P(j,i) = 0.5 * ( p(j+1,i+1) + p(j+2,i+1) );
                velMag(j,i) = sqrt( U(j,i)^2 + V(j,i)^2 );
            end
        end

        % Divergence diagnostic (try, but do not fail the whole compute if absent)
        try
            div = calcDiv(div,u,v,nx,ny,dx,dy);
            [mval, ii, jj] = maxDiv(div,nx,ny);
        catch
            mval = NaN; ii = NaN; jj = NaN;
        end

        % Pack and return
        out = struct();
        out.Re = Re;
        out.params = struct('nx',nx,'ny',ny,'xe',xe,'ye',ye,'dx',dx,'dy',dy,'Utop',Utop,'dt',dt,'iter',iter);
        out.grid = struct('xc',xc,'yc',yc,'Xc',Xc,'Yc',Yc);
        out.fields = struct('U',U,'V',V,'P',P,'velMag',velMag,'pFace',p,'uFace',u,'vFace',v);
        out.divergence = div;
        out.divMax = struct('value',mval,'i',ii,'j',jj);
        return

    case 'plot'
        % arg2 must be data struct
        data = arg2;
        if nargin < 5
            error('NewReSim_external:Plot requires data, plotChoice, axContour and axStream.');
        end
        plotChoice = string(arg3);
        axContour = arg4;
        axStream = arg5;

        % Validate data quickly
        if ~isstruct(data) || ~isfield(data,'fields')
            error('NewReSim_external:Invalid data struct provided to plot mode.');
        end

        % Extract
        xc = data.grid.xc;
        yc = data.grid.yc;
        Xc = data.grid.Xc;
        Yc = data.grid.Yc;
        U = data.fields.U;
        V = data.fields.V;
        P = data.fields.P;
        velMag = data.fields.velMag;
        xe = data.params.xe;
        ye = data.params.ye;

        % Select plot data
        switch lower(plotChoice)
            case 'pressure'
                plotData = P; cbLabel = 'Pressure';
            case 'horizontal velocity'
                plotData = U; cbLabel = 'Horizontal velocity';
            case 'vertical velocity'
                plotData = V; cbLabel = 'Vertical velocity';
            case 'total velocity'
                plotData = velMag; cbLabel = 'Velocity magnitude';
            otherwise
                plotData = U; cbLabel = 'Horizontal velocity';
        end
  
        % Contour axis
        if isempty(axContour) || ~isvalid(axContour)
            axContour = gca;
        end
        cla(axContour);
        hold(axContour,'off');
        contourLevels = 50;
        contourf(axContour, Xc, Yc, plotData, contourLevels, 'LineStyle', 'none');
        box(axContour,'on');
        xlabel(axContour,'x');
        % Ensure ylabel not rotated
        yl = ylabel(axContour,'y');
        try
            set(yl,'Rotation',0);
        catch
            % older MATLAB versions may ignore; continue
        end
        xlim(axContour,[0 xe]);
        ylim(axContour,[0 ye]);
        pbaspect(axContour,[1 1 1]);
        cb = colorbar(axContour);
        cb.Label.String = cbLabel;
        try
            set(cb,'Rotation',0);
        catch
            % older MATLAB versions may ignore; continue
        end


        colormap(axContour,'jet');
        if all(isfinite(plotData(:)))
            caxis(axContour,[min(plotData(:)), max(plotData(:))]);
        end

        % Streamline axis
        if isempty(axStream) || ~isvalid(axStream)
            axStream = gca;
        end
        cla(axStream);
        hold(axStream,'on');
        box(axStream,'on');
 
        for k = 1:1
           
            
                streamslice(axStream, xc,yc,U,V);
          
        end

        xlabel(axStream,'x');
        yl2 = ylabel(axStream,'y');
        try
            set(yl2,'Rotation',0);
        catch
        end
        % Deliberately omit title
        xlim(axStream,[0 xe]);
        ylim(axStream,[0 ye]);
        pbaspect(axStream,[1 1 1]);
        hold(axStream,'off');

        % Plot mode returns nothing
        out = [];
        return

    otherwise
        error('NewReSim_external:Unknown mode. Use ''compute'' or ''plot''.');
end
end











function dt = stabilityCondition(u,v,nx,ny,dx)
    maxVelMag = 0.0;
    for i=1:nx+2
        for j=1:ny+2
            velMag = sqrt(u(j,i)^2 + v(j,i)^2);
            if  velMag > maxVelMag
                maxVelMag = velMag;
            end
        end
    end

    dt = 0.3*dx/maxVelMag; % CFL stability condition
end

function [u,v,u_old,v_old] = boundaryConditions(u,v,u_old,v_old,ut,ub,ul,ur,vt,vb,vl,vr,nx,ny)
    u_old(:,2) = ul; % Left wall
    u_old(:,nx+2) = ur; % Right wall
    u_old(ny+2,:) = 2*ut - u_old(ny+1,:); % Top wall
    u_old(1,:) = 2*ub - u_old(2,:); % Bottom wall

    v_old(ny+2,:) = vt; % Top wall
    v_old(2,:) = vb; % Bottom wall
    v_old(:,1) = 2.0*vl - v_old(:,2); % Left wall
    v_old(:,nx+2) = 2.0*vr - v_old(:,nx+1); % Right wall

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

function [p,a_p] = pressureSolve(u,v,rho,nx,ny,dx,dy,dt,iter)
    [p,a_p] = sor(u,v,rho,nx,ny,dx,dy,dt,iter);
end

function [p,a_p] = sor(u,v,rho,nx,ny,dx,dy,dt,iter)
    p = zeros(ny+2,nx+2);
    rhs = zeros(ny+2,nx+2);

    % Setting pressure coefficients
    a_p = zeros(ny+2,nx+2);
    a_e = ones(ny+2,nx+2)/rho/dx/dx;
    a_w = ones(ny+2,nx+2)/rho/dx/dx;
    a_n = ones(ny+2,nx+2)/rho/dy/dy;
    a_s = ones(ny+2,nx+2)/rho/dy/dy;

    a_e(:,nx+1) = 0.0;
    a_w(:,2) = 0.0;
    a_n(ny+1,:) = 0.0;
    a_s(2,:) = 0.0;

    a_p = -(a_e + a_w + a_n + a_s);

    maxError = Inf;
    TOL = 1e-9;
    counter = 0;
    beta = 1.872;

    while(abs(maxError) > TOL)
        counter = counter + 1;
        maxError = 0;
        
        % Solve for pressure field
        for i=2:nx+1
            for j=2:ny+1
                rhs(j,i) =(u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy;
                rhs(j,i) = rhs(j,i)/dt -(a_w(j,i)*p(j,i-1) + a_e(j,i)*p(j,i+1) + a_n(j,i)*p(j+1,i) + a_s(j,i)*p(j-1,i));
                p(j,i) = beta*rhs(j,i)/a_p(j,i) +(1-beta)*p(j,i);
            end
        end
        
        % Determine residuals
        for i=2:nx+1
            for j=2:ny+1
                rhs(j,i) =(u(j,i+1) - u(j,i))/dx +(v(j+1,i) - v(j,i))/dy;
                error = a_w(j,i)*p(j,i-1) + a_e(j,i)*p(j,i+1) + a_n(j,i)*p(j+1,i) + a_s(j,i)*p(j-1,i) + a_p(j,i)*p(j,i) - rhs(j,i)/dt;
                if abs(error) > maxError
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

function [div] = calcDiv(div,u,v,nx,ny,dx,dy)
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