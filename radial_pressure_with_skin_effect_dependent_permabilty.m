clear
clc

%%% Initial and boundary conditions:

P_wf = 24e6; %[Pa], Dirichlet BC, stops the loop before 2 phases (PVT data for the reservoir fluid)
P_initial = 34.475e6; %[Pa], Dirichlet BC
P_bottomhole = 0.7 * P_initial; %[Pa], Dirichlet BC, regulated from surface
days = 300;

%%% Rock and fluids:

phi = 0.12;
S_oi = 0.87; 
S_wi = 0.13; 
B_oi = 1.25;
C_w = 4.351e-10; %[1/Pa]
C_r = 7.252e-10; %[1/Pa]
C_t = C_w + C_r; %[1/Pa] changes in respect to P neglectable due to scale of reservoir
mi = 5e-3; %[Pa*s], isothermic case
k = 10^-13; %[m^2], sandstone
s = 2; % skin effect ∈ <-5, 10>, can be defined as a vector to show how diffrent skin factors affect the pressure drop

%%% Spatial properties:

h_well_reservoir = 80; %[m]
A = 1200^2; % [m^2] documented reservoir area
r_well = 0.127; %[m]
r_e = sqrt((A)/pi); %[m]
r_skin = 0.6; %[m], pressure buildup test estimate
Volume_phi = 1.445e7; %[m^3]
Volume_phi_oil = Volume_phi * S_oi; %[m^3]

%%% Hawkins' model for change in localized permability due to skin effect:

k_skin = ((s/(k*log(r_skin/r_well))+1))^-1;
display(k_skin, 'Permability at skin affected radius [m^2]');
display(k,'Base permability [m^2]');

%%% Grid

Nr = 100;
Nt = 300000;
t_max = days*24*3600; %large timestep allowed due to nature of pressure drop + for stability of the explicit FDM scheme
dr = (r_e - r_well)/(Nr-1);%
dt = t_max / Nt;
r = linspace(r_well, r_e, Nr)';
P = P_initial * ones(Nr, 1);

%%% Coefficients:

D_base = k/(phi*mi*C_t);
D_skin_effect = k_skin/(phi*mi*C_t);
D = D_base*ones(Nr,1);
display(D_base);
display(D_skin_effect);

%%% Stability condition for the explicit scheme:

D_max = max(D);
dt_stable = (dr^2) / (2 * D_max);

if dt > dt_stable
    warning(['Timestep dt = ', num2str(dt), ...
        ' exceeds stability limit for explicit scheme: dt < ', num2str(dt_stable), ...
        '. Results may be unstable.']);
end

%%% Discretized matrix coefficients:

A_mat = zeros(Nr, Nr);
 for i = 2:Nr-1 % backwards spatial step
 ri = r(i);
 if r(i)<=r_skin
     D(i) = D_skin_effect;
 else
     D(i) = D_base;
 end
 alpha1 = D(i) * dt / (2*dr^2) - D(i) * dt / (4*dr*ri);
 alpha2 = 1 + D(i) * dt / dr^2;
 alpha3 = D(i) * dt / (2*dr^2) + D(i) * dt / (4*dr*ri);
 A_mat(i, i-1) = -alpha1;
 A_mat(i, i)   = alpha2;
 A_mat(i, i+1) = -alpha3;
 end
 A_mat(1,:) = 0; A_mat(1,1) = 1;
 A_mat(end,end-1) = -1; A_mat(end,end) = 1;
 
%%% Solution for time dependent pressure:
 
P_all = zeros(Nr, Nt);
 P_all(:,1) = P;
 
for n = 2:Nt % timestep - n;
    RHS = P;
    RHS(1) = P_wf;
    RHS(end) = 0;
    P = A_mat \ RHS;
    P_all(:,n) = P;
   
end

%%% Gif plot:
 
figure('Color','w');
filename = 'pressure_diffusion.gif';

for n = 1:5000:Nt
    plot(r, P_all(:,n)/1e6, 'b-', 'LineWidth', 2);
    xlabel('r [m]', 'FontWeight','bold');
    ylabel('P [MPa]', 'FontWeight','bold');
    title(sprintf('Radial Pressure Profile  t = %.1f days', n*dt/86400), 'FontWeight','bold');
    ylim([P_wf/1e6 - 1, P_initial/1e6 + 1]);
    grid on;
    set(gca, 'FontSize', 11);

    drawnow;

    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    if n == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.01);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
    end
end

%%% "Heatmap" plot:

time = (0:Nt-1) * dt / 86400;  % time vector in [days]
skip = 5000;
P_plot = P_all(:,1:skip:end)/1e6;      % downsample in time
time_plot = (0:skip:Nt-1) * dt / 86400; % [days]
imagesc(time_plot, r, P_plot);
set(gca, 'YDir', 'normal'); % Flip Y axis
xlabel('Time [days]');
ylabel('Radius [m]');
title('Pressure Evolution in Space and Time');
colorbar;

%%% Pressure drop at diffrent radii plot:
r_indices = [round(Nr/4), round(Nr/2), round(3*Nr/4), Nr];
colors = lines(length(r_indices));

figure;
hold on;
for i = 1:length(r_indices)
    plot(time, P_all(r_indices(i), :) / 1e6, 'LineWidth', 1.5, 'Color', colors(i,:));
end
xlabel('Time [days]');
ylabel('Pressure [MPa]');
title('Pressure Evolution at Different Radii');
legend(arrayfun(@(ri) sprintf('r = %.1f m', r(ri)), r_indices, 'UniformOutput', false));
grid on;

%%% Cumulative oil produced over time with Bo(P)
c_o = 1e-6; % [1/Pa], oil compressibility
Np = zeros(1, Nt); % Preallocate
time = (0:Nt-1) * dt / 86400;  % [days]

for n = 1:Nt
    P_avg = mean(P_all(:, n));
    Bo_dynamic = B_oi * (1 + c_o * (P_initial - P_avg)); % linear approx for Bo(P)
    Np(n) = (phi * Volume_phi * (P_initial - P_avg)) / (C_t * Bo_dynamic * P_initial);
end

%%% Plot
figure;
plot(time, Np / 1e6, 'k-', 'LineWidth', 2);
xlabel('Time [days]');
ylabel('Cumulative Oil Produced [10^6 m^3]');
title('Cumulative Oil Produced vs. Time (with Bo(P) dependency)');
grid on;

