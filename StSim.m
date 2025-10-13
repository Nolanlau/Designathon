function StSim(ax, St)



if ~(isa(ax,'matlab.ui.control.UIAxes') || isa(ax,'matlab.graphics.axis.Axes')) || ~isvalid(ax)
    error('StSim:BadAxes','First input must be a valid UIAxes/Axes handle.');
end

%% --- User parameters ---

animation_delay = 0.05/5;
remove_particles_on_exit = true;

%% --- Domain and base flow ---
Lx = 2.0; Ly = 1.0;
nx = 201; ny = 101;
x = linspace(-Lx/2, Lx/2, nx);
y = linspace(-Ly/2, Ly/2, ny);
[XX, YY] = meshgrid(x, y);

U = 1.0;
Gamma = 3.0;
xc = 0.2; yc = 0.0;
sigma = 0.15;
L_char = Lx;

% Streamfunction and velocity (Psi unused later)
% Psi = U*YY + Gamma*exp(-((XX-xc).^2 + (YY-yc).^2)/(2*sigma^2));
u = U + Gamma.*( -(YY-yc)./sigma^2 ).*exp(-((XX-xc).^2 + (YY-yc).^2)/(2*sigma^2));
v =        Gamma.*(  (XX-xc)./sigma^2 ).*exp(-((XX-xc).^2 + (YY-yc).^2)/(2*sigma^2));

Fx = griddedInterpolant({x,y}, u.', 'linear', 'nearest'); % note transpose
Fy = griddedInterpolant({x,y}, v.', 'linear', 'nearest');

tau_p = St * L_char / U;

%% --- Particles ---
Np  = 24;
p0x = -Lx/2 + 0.02;
p0y = linspace(-0.8*Ly/2, 0.8*Ly/2, Np)';

xp = p0x * ones(Np,1);
yp = p0y;
vxp = zeros(Np,1);
vyp = zeros(Np,1);
for i = 1:Np
    vxp(i) = Fx(xp(i), yp(i));
    vyp(i) = Fy(xp(i), yp(i));
end
active = true(Np,1);

%% --- Time stepping ---
T_char = L_char / U;
t_end  = 4 * T_char;
dt     = 0.002 * T_char;
tvec   = 0:dt:t_end;
nsteps = numel(tvec);

frame_dt = 1e-2 * T_char;
steps_per_frame = max(1, round(frame_dt/dt));
trail_length_seconds = 0.5 * T_char;
trail_steps = max(1, round(trail_length_seconds / dt));

%% --- Prepare target axes (no new figure) ---
cla(ax);
hold(ax,'on'); box(ax,'on');
axis(ax,'equal');
xlim(ax,[-Lx/2 Lx/2]); ylim(ax,[-Ly/2 Ly/2]);
xlabel(ax,'x'); ylabel(ax,'y');

density = 2;
hs = streamslice(ax, XX, YY, u, v, density, 'cubic');

for h = hs(:).'
    try
        set(h, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.9);
    catch
        try set(h, 'FaceColor', [0 0 0]); catch, end
    end
end

% Particle graphics
particle_lines   = gobjects(Np,1);
particle_markers = gobjects(Np,1);
for i = 1:Np
    particle_lines(i)   = plot(ax, NaN, NaN, '-', 'LineWidth', 1.5, 'Color', 'r');
    particle_markers(i) = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', 'w', ...
                               'MarkerEdgeColor', 'k', 'MarkerSize', 6);
end
drawnow;

%% --- Dynamics helpers ---
accel = @(xq,yq,vxq,vyq) deal( (Fx(xq,yq) - vxq)./tau_p, (Fy(xq,yq) - vyq)./tau_p );

trail_x   = NaN(Np, trail_steps);
trail_y   = NaN(Np, trail_steps);
trail_idx = 1;

%% --- Main loop ---
for n = 1:nsteps-1
    X=xp; Y=yp; VX=vxp; VY=vyp;

    [ax1, ay1] = accel(X,Y,VX,VY);
    k1x = VX;                k1y = VY;                k1vx = ax1;                k1vy = ay1;

    X2 = X + 0.5*dt*k1x;     Y2 = Y + 0.5*dt*k1y;
    VX2= VX+ 0.5*dt*k1vx;    VY2= VY+ 0.5*dt*k1vy;
    [ax2, ay2] = accel(X2,Y2,VX2,VY2);
    k2x = VX2;               k2y = VY2;               k2vx = ax2;                k2vy = ay2;

    X3 = X + 0.5*dt*k2x;     Y3 = Y + 0.5*dt*k2y;
    VX3= VX+ 0.5*dt*k2vx;    VY3= VY+ 0.5*dt*k2vy;
    [ax3, ay3] = accel(X3,Y3,VX3,VY3);
    k3x = VX3;               k3y = VY3;               k3vx = ax3;                k3vy = ay3;

    X4 = X + dt*k3x;         Y4 = Y + dt*k3y;
    VX4= VX+ dt*k3vx;        VY4= VY+ dt*k3vy;
    [ax4, ay4] = accel(X4,Y4,VX4,VY4);
    k4x = VX4;               k4y = VY4;               k4vx = ax4;                k4vy = ay4;

    xp  = X  + dt*(k1x + 2*k2x + 2*k3x + k4x)/6;
    yp  = Y  + dt*(k1y + 2*k2y + 2*k3y + k4y)/6;
    vxp = VX + dt*(k1vx+ 2*k2vx+ 2*k3vx+ k4vx)/6;
    vyp = VY + dt*(k1vy+ 2*k2vy+ 2*k3vy+ k4vy)/6;

    if remove_particles_on_exit
        exited = xp < -Lx/2 | xp > Lx/2 | yp < -Ly/2 | yp > Ly/2;
        active(exited) = false;
    end

    if any(active)
        active_idx = find(active);
        trail_x(active_idx, trail_idx) = xp(active_idx);
        trail_y(active_idx, trail_idx) = yp(active_idx);
    end
    if any(~active)
        inactive_idx = find(~active);
        trail_x(inactive_idx, trail_idx) = NaN;
        trail_y(inactive_idx, trail_idx) = NaN;
    end

    trail_idx = trail_idx + 1;
    if trail_idx > trail_steps, trail_idx = 1; end

    if mod(n, steps_per_frame) == 0
    % Skip updates if axes were deleted (user changed tab, closed figure, etc.)
    if ~isvalid(ax)
        return
    end

    idxs = mod((trail_idx:trail_idx+trail_steps-1)-1, trail_steps) + 1;
    for i = 1:Np
        % If line object is invalid (deleted by cla or otherwise), recreate it
        if ~isvalid(particle_lines(i))
            particle_lines(i) = plot(ax, NaN, NaN, '-', 'LineWidth', 1.5, 'Color', 'r');
        end
        if ~isvalid(particle_markers(i))
            particle_markers(i) = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', 'w', ...
                                       'MarkerEdgeColor', 'k', 'MarkerSize', 6);
        end

        tx = trail_x(i, idxs); ty = trail_y(i, idxs);
        valid = ~isnan(tx);
        if any(valid)
            set(particle_lines(i), 'XData', tx(valid), 'YData', ty(valid));
        else
            set(particle_lines(i), 'XData', NaN, 'YData', NaN);
        end
        if active(i)
            set(particle_markers(i), 'XData', xp(i), 'YData', yp(i));
        else
            set(particle_markers(i), 'XData', NaN, 'YData', NaN);
        end
    end
    drawnow limitrate nocallbacks;
    pause(animation_delay);
end

        drawnow limitrate;
        pause(animation_delay);
    end
end

