function resultsTable = MoleCrabLegPairRFT_v3(phi_d)
    %% Robot Parameters
    l1  = 12/10;   % Length of Link 1 (motor-driven arm)
    l2  = 75.13/10;% Length of Link 2
    h_m = 34.5/10; % height of motor center from hinge point
    w_b = 33.66/10; 
    w_t = 9.8/10;  
    nEl = 200;  
    Lmax = l2 - (abs(h_m) - l1); 

    %% Motor angle 
    revs = 1;
    nFrames = revs*1000;

    % Convert user-specified phi_d (degrees) to radians
    phi = deg2rad(phi_d);

    theta_m = linspace(-pi/2 - phi, revs*(3*pi/2 - phi), nFrames);


    % New center
    x0_new = h_m * sin(phi);
    y0_new = h_m * cos(phi);

    %% Link 1 endpoints
    x1 = l1.*cos(theta_m) + x0_new;
    y1 = l1.*sin(theta_m) + y0_new;

    %% Slope to fixed hinge point
    theta_2 = atan2((y1 - 0), (x1 - 0));

    %% Link 2 endpoints
    x2 = x1 - l2.*cos(theta_2);
    y2 = y1 - l2.*sin(theta_2);

    %% RFT Coeffs (generic for plastic pellets)
    A_coeff = [0.206, 0.169];
    B_coeff = [0.212, 0.358, 0.055];
    C_coeff = [-0.124, 0.253, 0.007];
    D_coeff = 0.088;
    zeta    = 0.111; 

    %% Arrays to store net forces
    F_xAll = zeros(1, nFrames);
    F_zAll = zeros(1, nFrames);

    % For storing base width at each frame (optional)
    w_b_all = zeros(1, nFrames);
    w_x_all = zeros(nEl, nFrames);

    % Initialize for the first frame
    ext = sqrt((x2(1) - 0)^2 + (y2(1) - 0)^2);
    y2_init = linspace(0, -abs(ext*sin(theta_2(1))), nEl);
    x2_init = linspace(0, -abs(ext*cos(theta_2(1))), nEl);

    %% 1) Create Figure for Animation
    figAnim = figure('Name','Mole Crab animation','NumberTitle','off',...
           'Position',[200, 0, 1400, 1000]);
    set(gcf,'WindowState','maximized');
    axis equal;
    hold on; grid on;
    xlim([-(l2 - h_m + 1.5), (h_m + l1)]);
    ylim([-(l2 + l1 - h_m), (l1 + h_m)]);
    xlabel('X','FontSize',12);
    ylabel('Y','FontSize',12);
    title("Mole crab's Leg movement ",'FontSize',14);

    link1Plot = plot([x0_new x1(1)], [y0_new y1(1)],...
                     '-o','LineWidth',1,'Color',[0 0.6 0.8],...
                     'MarkerFaceColor','b','MarkerSize',4);

    link2Plot = plot([x1(1) x2(1)], [y1(1) y2(1)],...
                     '-o','LineWidth',1,'Color',[0.8 0.4 0],...
                     'MarkerFaceColor','r','MarkerSize',4);

    plot(0, 0, 'ks','MarkerSize',4,'MarkerFaceColor','k');

    quiverF_net = quiver(nan, nan, nan, nan, 0, ...
        'Color','m','LineWidth',0.5, ...
        'MaxHeadSize',0.01, 'AutoScale','on');

    scaleX = 50;
    scaleZ = 50;

    ax = gca;
    set(ax, 'Position', [0.05, 0.05, 0.9, 0.9]);

    %% 2) Initialize Video Writer with MP4
    v = VideoWriter('MoleCrabAnimation.mp4','MPEG-4'); 
    v.FrameRate = 30;  % adjust as needed
    open(v);

    %% 3) Animation Loop
    for ii = 1:nFrames
        
        % (A) leg extension
        ext = sqrt((x2(ii) - 0)^2 + (y2(ii) - 0)^2);

        % (B) discretize
        leg_elements   = linspace(0, ext, nEl);
        element_length = leg_elements(2) - leg_elements(1);

        % (C) base width interpolation
        w_b_c = w_t + ((w_b - w_t)*ext)/Lmax;
        w_b_all(ii) = w_b_c;
        w_x = ((w_t - w_b_c)/ext).*leg_elements + w_b_c;
        w_x_all(:,ii) = w_x';

        dA = w_x .* element_length;

        % (D) pen depth
        theta_2_c = theta_2(ii);
        pen_depth = leg_elements .* sin(theta_2_c);
        pen_depth(pen_depth < 0) = 0;

        % positions for each element
        frac = leg_elements / ext;
        xEl = frac * x2(ii);
        yEl = frac * y2(ii);

        % displacement
        del_x2 = xEl - x2_init;
        del_y2 = yEl - y2_init;

        x2_init = xEl;
        y2_init = yEl;

        gamma = -(pi - mean(atan2(del_y2, del_x2)));
        beta  = -atan((yEl(2)-yEl(1)) / (xEl(2)-xEl(1)));

        b = beta / pi;
        g = gamma / (2*pi);

        alpha_z = ...
                A_coeff(1)*cos(2*pi*(0*b + 0*g)) + ...
                A_coeff(2)*cos(2*pi*(1*b + 0*g)) + ...
                B_coeff(1)*sin(2*pi*(1*b + 1*g)) + ...
                B_coeff(2)*sin(2*pi*(0*b + 1*g)) + ...
                B_coeff(3)*sin(2*pi*(-1*b + 1*g));

        alpha_x = ...
            C_coeff(1)*cos(2*pi*(1*b + 1*g)) + ...
            C_coeff(2)*cos(2*pi*(0*b + 1*g)) + ...
            C_coeff(3)*cos(2*pi*(-1*b + 1*g)) + ...
            D_coeff*sin(2*pi*(1*b + 0*g));

        alpha_z = alpha_z * zeta;
        alpha_x = alpha_x * zeta;

        % forces
        f_z_elements = 2 * alpha_z .* pen_depth .* dA;
        f_x_elements = 2 * alpha_x .* pen_depth .* dA;

        Fx_net = sum(f_x_elements);
        Fz_net = sum(f_z_elements);

        % store net forces
        F_xAll(ii) =  Fx_net; 
        F_zAll(ii) =  -Fz_net;

        % quiver
        fx_quiver = scaleX .* f_x_elements;
        fz_quiver = scaleZ .* f_z_elements;

        % update link1, link2
        set(link1Plot, 'XData',[x0_new, x1(ii)], 'YData',[y0_new, y1(ii)]);
        set(link2Plot, 'XData',[x1(ii), x2(ii)], 'YData',[y1(ii), y2(ii)]);

        set(quiverF_net, 'XData', xEl, 'YData', yEl, ...
                         'UData', fx_quiver, 'VData', fz_quiver, ...
                         'AutoScale','off','AutoScaleFactor',1);

        drawnow;

        % (E) capture frame for video
        frame = getframe(figAnim);
        writeVideo(v, frame);

        % pause(0.001);
    end

    % close the video file
    close(v);

    %% 4) Final Plot of Fz vs. Fx
    theta_p = theta_m + pi/2;

    figFinal = figure('Name','Forces vs Theta','NumberTitle','off');
    subplot(2,1,1);
    plot(rad2deg(theta_p), F_zAll, 'LineWidth',1.5);
    xlabel('\theta_p (deg)');
    ylabel('F_z');
    title('F_z vs \theta_p');
    grid on;

    subplot(2,1,2); 
    plot(rad2deg(theta_p), F_xAll, 'LineWidth',1.5);
    xlabel('\theta_p (deg)');
    ylabel('F_x');
    title('F_x vs \theta_p');
    grid on;

    %% 5) Save the final figure as PNG
    saveas(figFinal, 'FinalPlot.png');

    %% 6) Torque along y axis
    T_y = F_zAll .* (11.68/1000);

    % Build resultsTable
    thetaDeg = rad2deg(theta_m(:));
    thetaDeg_p = rad2deg(theta_p(:));
    Fz = F_zAll(:);
    Fx = F_xAll(:);
    Ty = T_y(:);

    resultsTable = table(thetaDeg, thetaDeg_p, Fz, Fx, Ty, ...
        'VariableNames', {'thetaDeg', 'thetaDeg_p', 'Fz','Fx','Ty'});
end
