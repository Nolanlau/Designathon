function CaSim(app, CaInput)
        % -------- Inputs & targets --------
        Ca = max(CaInput, 1e-8);
        ax = app;

        % -------- Domain & grid --------
        L  = 1.0;
        Nx = 64;                   
        Ny = 128;
        x = (0:Nx-1)/Nx*L;  y = (0:Ny-1)/Ny*L;
        [XX,YY] = meshgrid(x,y);

        % Spectral wavenumbers
        kx = (2*pi/L)*[0:(Nx/2-1), (-Nx/2):-1];
        ky = (2*pi/L)*[0:(Ny/2-1), (-Ny/2):-1];
        [KX,KY] = meshgrid(kx,ky);
        K2 = KX.^2 + KY.^2;  K4 = K2.^2;
        K2(1,1) = 1; K4(1,1) = 1;   

        % -------- Physics params --------
        epsilon = 0.02;             
        M       = 1e-3;              
        mu      = 1.0;               
        sigma   = 1e-4;             
        U0      = Ca*sigma/mu;       
        R       = 0.18;              
        center  = [0.5, 0.5];

        % sigma ≈ (2√2/3) * lambda / epsilon  → lambda
        lambda = (3*epsilon*sigma)/(2*sqrt(2));

        % -------- Time stepping --------
        dt           = 2e-4;
        nsteps       = 2400;         
        plotInterval = 30;           

        % -------- Initial condition --------
        rr  = sqrt( (XX-center(1)).^2 + (YY-center(2)).^2 );
        phi = tanh( (R - rr)/(sqrt(2)*epsilon) );   % +1 inside, -1 outside

        % Semi-implicit denominator
        denom = 1 + dt*M*lambda*epsilon^2 .* K4;

        % -------- Axes aesthetics --------
        cla(ax); hold(ax,'on');
        ax.XLim = [0 1]; ax.YLim = [0 1];
        axis(ax,'equal','tight');
        xlabel(ax,'x'); ylabel(ax,'y');
        %title(ax, sprintf('Phase field   Ca = %.3g', Ca),'FontWeight','bold');
        ax.Color = 'w'; grid(ax,'off'); ax.Layer = 'top'; box(ax,'on');

        droplet_color   = [0.3010, 0.7450, 0.9330];
        custom_colormap = [1,1,1; droplet_color];   
        colormap(ax, custom_colormap);

        % -------- Time loop --------
        for n = 1:nsteps
            % divergence-free shear-like velocity
            u_x = U0 * sin(2*pi*YY/L);
            u_y = zeros(size(u_x));

            % chemical potential nonlinear part
            fphi = lambda*(phi.^3 - phi);

            % spectral transforms
            phi_hat = fft2(phi);
            fhat    = fft2(fphi);

            % convection term u·∇phi
            dphi_dx = real(ifft2(1i*KX .* phi_hat));
            dphi_dy = real(ifft2(1i*KY .* phi_hat));
            conv    = u_x .* dphi_dx + u_y .* dphi_dy;
            conv_hat = fft2(conv);

            % semi-implicit update
            rhs_hat     = phi_hat + dt*( -conv_hat - M*(K2 .* fhat) );
            phi_hat_new = rhs_hat ./ denom;
            phi         = real(ifft2(phi_hat_new));

            % clip
            phi(phi>1) = 1;  phi(phi<-1) = -1;

           
            if mod(n,plotInterval)==0 || n==1 || n==nsteps
                cla(ax);
                contourf(ax, XX, YY, phi, 1, 'LineStyle','none'); hold(ax,'on');
                contour(ax,  XX, YY, phi, [0 0], 'k', 'LineWidth', 1.2);

               
                y_base = 0.85; dy = 0.02;
                y_lines = [y_base-dy, y_base, y_base+dy];
                x_start = 0; x_end = 1; x_mid = 0.5;
                for i = 1:numel(y_lines)
                    plot(ax, [x_start x_end], [y_lines(i) y_lines(i)], 'k', 'LineWidth',1);
                    quiver(ax, x_mid, y_lines(i), -0.05*(x_end-x_start), 0, 0, ...
                        'k','MaxHeadSize',2,'LineWidth',1,'AutoScale','off');
                end

               
                y_base2 = 0.15; dy2 = 0.02;
                y_lines2 = [y_base2-dy2, y_base2, y_base2+dy2];
                for i = 1:numel(y_lines2)
                    plot(ax, [x_start x_end], [y_lines2(i) y_lines2(i)], 'k', 'LineWidth',1);
                    quiver(ax, x_mid, y_lines2(i),  0.05*(x_end-x_start), 0, 0, ...
                        'k','MaxHeadSize',2,'LineWidth',1,'AutoScale','off');
                end

                title(ax, sprintf('Time = %.3f', n*dt), 'FontWeight','bold');
                drawnow limitrate;
            end
        end