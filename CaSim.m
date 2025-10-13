% phasefield_capillary_single.m
% Single-case demonstration of Capillary number influence on a 2D droplet.
% Phase-field (Cahn-Hilliard) with frozen prescribed velocity field.
% Pseudo-spectral in space, semi-implicit in time. Streamslice overlay.

clear; close all; clc;

% ca= 10 to 10000

Ca=9000;


factor=Ca/10000;
%% Domain and discretisation
L = 1.0;            % domain length
Nx = 128/2; Ny = 128; % grid resolution (prefer powers of two)
x = (0:Nx-1)/Nx*L;  y = (0:Ny-1)/Ny*L;
[XX,YY] = meshgrid(x,y);

% Spectral wavenumbers for fft2 (MATLAB ordering)
kx = (2*pi/L)*[0:(Nx/2-1) (-Nx/2):-1];
ky = (2*pi/L)*[0:(Ny/2-1) (-Ny/2):-1];
[KX,KY] = meshgrid(kx,ky);
K2 = KX.^2 + KY.^2;
K4 = K2.^2;
K2(1,1) = 1; K4(1,1) = 1; % prevent NaN, harmless since multiplied by zero where needed

%% Phase-field and physical parameters
epsilon = 0.02;   % interface thickness
M = 1e-3;         % mobility
mu = 1.0;         % dynamic viscosity (reference) 
U0 = 1*factor;        % velocity amplitude (change to vary Ca)
R = 0.18;         % initial droplet radius
center = [0.5,0.5];
sigma=1e-4;

Ca = mu*U0/sigma
% Map sigma to lambda via approximate relation sigma ~ (2*sqrt(2)/3)*lambda/epsilon
lambda = (3*epsilon*sigma)/(2*sqrt(2));

fprintf('Capillary number Ca = %.3g (mu=%g, U0=%g, sigma=%g)\n', Ca, mu, U0, sigma);

%% Time stepping
dt = 2e-4;        % time step; reduce if instability occurs
nsteps = 2400;    % total steps
plotInterval = 30; % plotting frequency

% Initial condition: smooth circle via tanh
rr = sqrt( (XX-center(1)).^2 + (YY-center(2)).^2 );
phi = tanh( (R - rr)/(sqrt(2)*epsilon) ); % +1 inside, -1 outside

% Precompute denominator for semi-implicit update
denom = 1 + dt*M*lambda*epsilon^2 .* K4;

% Prepare figure
figure('Units','normalized','Position',[0.12 0.12 0.7 0.7]);


for n = 1:nsteps %https://uk.mathworks.com/help/matlab/ref/annotation.html

    % Prescribed, divergence-free, periodic velocity field:
    % horizontal shear-like profile: u_x = U0 * sin(2*pi*y/L)
    u_x = U0 * sin(2*pi*YY/L);
    u_y = zeros(size(u_x));

    % Chemical potential nonlinear part f(phi) = lambda*(phi^3 - phi)
    fphi = lambda*(phi.^3 - phi);

    % Spectral transforms
    phi_hat = fft2(phi);
    fhat   = fft2(fphi);

    % Compute convective term via spectral derivatives
    dphi_dx = real(ifft2(1i*KX .* phi_hat));
    dphi_dy = real(ifft2(1i*KY .* phi_hat));
    conv = u_x .* dphi_dx + u_y .* dphi_dy;
    conv_hat = fft2(conv);

    % RHS and semi-implicit update
    rhs_hat = phi_hat + dt*( - conv_hat - M * (K2 .* fhat) );
    phi_hat_new = rhs_hat ./ denom;
    phi = real(ifft2(phi_hat_new));

    % mild clipping to [-1,1]
    phi(phi>1) = 1; phi(phi<-1) = -1;

% Visualization
   if mod(n,plotInterval) == 0 || n==1 || n==nsteps
        % clf;
        
        % --- 1. Define a custom colormap and figure color ---
        droplet_color = [0.3010, 0.7450, 0.9330]; % A light blue color
        % First color (white) is for the background, second (blue) is for the droplet
        custom_colormap = [1, 1, 1; droplet_color]; 
        colormap(custom_colormap);
        set(gcf, 'color', 'w'); % Also set the figure background to white
   
        % --- 2. Plot the filled droplet area ---
        % This fills the background with the first colormap color (white) and the
        % droplet with the second colormap color (blue).
        contourf(XX, YY, phi, 1, 'LineStyle', 'none');
        
        hold on; % <-- ADD THIS LINE to layer the next plot on top
        
        % --- 3. Plot the droplet interface line on top of everything ---
        % This draws a black line at the phi=0 contour.
        contour(XX, YY, phi, [0 0], 'k', 'LineWidth', 1.2); % <-- UNCOMMENT THIS LINE
        % --- 4. Draw three horizontal lines around y = 0.75 with arrows pointing left ---
        y_base = 0.85;               % Central y position
        dy = 0.02;                   % Vertical spacing between lines
        x_start = min(XX(:));
        x_end   = max(XX(:));
        x_mid   = 0.5 * (x_start + x_end);
        % --- 5. Draw three horizontal lines around y = 0.1 with arrows pointing right ---
        y_base2 = 0.15;               % Central y position for lower set
        dy2 = 0.02;                  % Vertical spacing
        y_lines2 = [y_base2 - dy2, y_base2, y_base2 + dy2];

        for i = 1:length(y_lines2)
            % Draw horizontal line
            plot([x_start x_end], [y_lines2(i) y_lines2(i)], 'k', 'LineWidth', 1);

            % Draw arrow pointing right at the mid-point
            quiver(x_mid, y_lines2(i), 0.05*(x_end - x_start), 0, 0, ...
                'k', 'MaxHeadSize', 2, 'LineWidth', 1, 'AutoScale', 'off');
        end

        % Define y positions for the three lines
        y_lines = [y_base - dy, y_base, y_base + dy];

        for i = 1:length(y_lines)
            % Draw horizontal line
            plot([x_start x_end], [y_lines(i) y_lines(i)], 'k', 'LineWidth', 1);

            % Draw arrow pointing left at the mid-point
            quiver(x_mid, y_lines(i), -0.05*(x_end - x_start), 0, 0, ...
                'k', 'MaxHeadSize', 2, 'LineWidth', 1, 'AutoScale', 'off');
        end

        % --- Standard Plot Formatting ---
        axis equal tight;
        xlabel('x'); ylabel('y');
       
        title(sprintf('Phase field at t = %.3f   Ca = %.3g', n*dt, Ca));
        hold off;
        
        drawnow;
    end

end


