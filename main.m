% Clear workspace
close all; clear; clc;

% Define initial state vector for Earth with position and velocity
rE = [-2.5527e7, 1.4486e8, -6.5836e3]; % Position in km
vE = [-2.9823e1, -5.2819, 3.1683e-4]; % Velocity in km/s
VE = norm(vE); % Magnitude of Earth's velocity
mu = 1.327e11; % Gravitational parameter of the Sun in km^3/s^2
t0_E = 0; % Start time for Earth's orbit simulation
tf_E = 211 * 24 * 60 * 60; % End time (launch time) in seconds
tspan_E = t0_E : 24 * 60 * 60 : tf_E; % Time span for Earth's orbit integration
initial_conditions_E = [rE, vE]; % Initial conditions vector for Earth

% Define initial state vector for Mars with position and velocity
rM = [-1.9717e8, -1.329e8, 2.0529e6]; % Position in km
vM = [1.4449e1, -1.8018e1, -7.3209e-1]; % Velocity in km/s
VM = norm(vM); % Magnitude of Mars' velocity
t0_M = 0; % Start time for Mars' orbit simulation
tf_M = 414 * 24 * 60 * 60; % End time (arrival time) in seconds
tspan_M = t0_M : 24 * 60 * 60 : tf_M; % Time span for Mars' orbit integration
initial_conditions_M = [rM, vM]; % Initial conditions vector for Mars

% Hohmann Transfer orbit initial conditions
r1 = [91930000, -120900000, 5318]; % Position of Earth at launch in km
r2 = [-1583000, 234800000, 4960000]; % Position of Mars at arrival in km
t0_sc = 211 * 24 * 60 * 60; % Start time for spacecraft's transfer orbit
tf_sc = 414 * 24 * 60 * 60; % End time for spacecraft's transfer orbit
tspan_sc = t0_sc : 24 * 60 * 60 : tf_sc; % Time span for spacecraft's transfer orbit
v_earth = [23.2246, 17.9177, -8.8019e-4]; % Earth's velocity at launch in km/s
v_inf_E = 3.016984537409886; % Excess velocity at departure in km/s

% Calculate the direction for the transfer orbit injection
v_norm = cross(r1, r2); % Cross product of position vectors to find normal to plane of motion
dir_vnorm = v_norm ./ norm(v_norm); % Normalizing the direction vector
plane = -v_earth * 1 ./ (dir_vnorm * norm(v_earth)); % Calculating the plane of the transfer
dir_plane = plane ./ norm(plane); % Normalizing the direction of the plane
vel_transit = v_earth + (v_inf_E * dir_plane); % Calculating the transfer orbit velocity
initial_conditions_sc = [r1, vel_transit]; % Initial conditions for the spacecraft

% Numerical integration for Earth, Mars, and spacecraft trajectories
[t_E, X_E] = ode45(@(t_E, X_E) state_vector(X_E, mu), tspan_E, initial_conditions_E);
[t_sc, X_sc] = ode45(@(t_sc, X_sc) state_vector(X_sc, mu), tspan_sc, initial_conditions_sc);
[t_M, X_M] = ode45(@(t_M, X_M) state_vector(X_M, mu), tspan_M, initial_conditions_M);

% Trajectory correction loop to minimize distance to Mars at arrival
d = 1e7; % Initial distance to Mars far greater than target threshold
d_list = d; % Initialize list to track distances
ite = 0; % Iteration counter
while d > 5.5094e+03 && ite < 10000 % Loop until distance is below threshold or max iterations reached
    initial_conditions_sc = [r1, vel_transit]; % Update initial conditions for spacecraft
    [t_sc, X_sc] = ode45(@(t_sc, X_sc) state_vector(X_sc, mu), tspan_sc, initial_conditions_sc);
    dis = [X_M(end, 1) - X_sc(end, 1), X_M(end, 2) - X_sc(end, 2), X_M(end, 3) - X_sc(end, 3)]; % Displacement vector between Mars and spacecraft at arrival
    d = norm(dis); % Compute the norm of the displacement vector to get distance
    vel_transit = vel_transit + dis * 10^-9; % Adjust velocity based on distance (simple correction factor)
    d_list = [d_list, d]; % Append new distance to list
    ite = ite + 1; % Increment iteration counter
end

% Visualization: Combined plot of Earth, Mars, and spacecraft trajectories
filename = 'trajectory.gif';
colordef black;
h = figure;
hold on;
grid on;
plot3(X_E(:, 1), X_E(:, 2), X_E(:, 3), '-b'); % Plot Earth's trajectory in blue
plot3(X_M(:, 1), X_M(:, 2), X_M(:, 3), '-r'); % Plot Mars' trajectory in red
plot3(X_sc(:, 1), X_sc(:, 2), X_sc(:, 3), '-g'); % Plot spacecraft's transfer orbit in dashed green
% Mark initial and final positions for Earth, Mars, and the Sun
plot3(rE(1), rE(2), rE(3), '*b', 'MarkerSize', 6, 'MarkerFaceColor', '[0 0.4470 0.7410]');
plot3(rM(1), rM(2), rM(3), '*r', 'MarkerSize', 6, 'MarkerFaceColor', '[0.6350 0.0780 0.1840]');
plot3(r1(1), r1(2), r1(3), 'ob', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
plot3(r2(1), r2(2), r2(3), 'or', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
plot3(0, 0, 0, 'oy', 'MarkerSize', 18, 'MarkerFaceColor', 'y'); % Plot the Sun as a yellow marker
formattedText = {'\fontsize{12}\color{white}\bfTrajectory to Mars 2020'; '\fontsize{8}\color{white}\rmEclipJ2000 heliocentric frame'}; 
title(formattedText)
ylabel('Y coordinate (km)','FontSize',10)
xlabel('X coordinate (km)','FontSize',10)
zlabel('Z coordinate (km)','FontSize',10)
view(3); % Set the view to 3D
% Create a legend for the plot
lgd = legend('Earth orbit', 'Mars orbit', 'Transfer orbit', 'r_{Earth_{initial}}', 'r_{Mars_{initial}}', 'r_{Earth_{final}}', 'r_{Mars_{final}}', 'r_{Sun}', 'AutoUpdate', 'off');
lgd.FontSize = 10;
set(gcf, 'Units', 'Normalized', 'Outerposition', [0 0 1 1]); % Maximize figure window
set(gca, 'CameraPosition', [2226640057.566024, -3119521064.437298, 69463120.07632765]);

% Create a GIF from the plots
for i = 1:10:length(X_E(:, 1))
    plot3(X_E(i, 1), X_E(i, 2), X_E(i, 3), 'ob', 'MarkerSize', 4);
    plot3(X_M(i, 1), X_M(i, 2), X_M(i, 3), 'or', 'MarkerSize', 4);
    drawnow; % Update plot
    frame = getframe(h); % Capture plot as a frame
    im = frame2im(frame); % Convert frame to image
    [imind, cm] = rgb2ind(im, 256); % Convert image to indexed image
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'LoopCount', Inf); % Write first frame to file
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append'); % Append subsequent frames
    end
end

% Plot extra frames for Mars and spacecraft positions
for i = length(X_E(:, 1)):10:length(X_M(:, 1))
    plot3(X_M(i, 1), X_M(i, 2), X_M(i, 3), 'or', 'MarkerSize', 4); % Mars position update
    % Calculate spacecraft index offset due to different array lengths
    sc_index = i - length(X_E(:, 1)) + 1;
    if sc_index <= length(X_sc(:, 1))
        plot3(X_sc(sc_index, 1), X_sc(sc_index, 2), X_sc(sc_index, 3), '^g', 'MarkerSize', 4); % Spacecraft position update
    end
    drawnow; % Update plot
    frame = getframe(h); % Capture the plot as a frame
    im = frame2im(frame); % Convert the frame to image
    [imind, cm] = rgb2ind(im, 256); % Convert the image to indexed image
    imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append'); % Append the frame to the GIF
end
hold off; % Release plot hold

% Calculate the inclination of the spacecraft's arrival trajectory with respect to Mars' orbital plane
h = cross(X_M(end, 1:3), X_M(end, 4:6)); % Angular momentum vector of Mars
k_hat = [0.4461321045940511, -0.0558363658975008, 0.8932236257109472]'; % Unit vector in the z-direction
inclination = atan2d(norm(cross(h, k_hat)), dot(h, k_hat)); % Calculate inclination angle

% The following sections are commented out and can be enabled for additional analysis or visualization

% % Vector diagram for departure
% colordef black; % Set background color to black
% figure; % Create a new figure
% hold on; % Hold on to plot multiple data sets
% grid on; % Enable grid for better visualization
% plot3(r1(1), r1(2), r1(3), 'ob', 'MarkerSize', 7, 'MarkerFaceColor', 'b'); % Plot Earth's position
% % Plot vectors for Earth's velocity, transfer orbit initial velocity, and departure velocity
% quiver3(r1(1), r1(2), r1(3), v_inf_E * dir_plane(1), v_inf_E * dir_plane(2), v_inf_E * dir_plane(3), '-c');
% quiver3(r1(1), r1(2), r1(3), vel_transit(1), vel_transit(2), vel_transit(3), '-g');
% quiver3(r1(1), r1(2), r1(3), v_earth(1), v_earth(2), v_earth(3), '-b');
% hold off; % Release plot hold

% % Vector diagram for arrival
% colordef black; % Set background color to black
% figure; % Create a new figure
% hold on; % Hold on to plot multiple data sets
% grid on; % Enable grid for better visualization
% plot3(r2(1), r2(2), r2(3), 'or', 'MarkerSize', 7, 'MarkerFaceColor', 'r'); % Plot Mars' position
% % Plot vectors for spacecraft's final velocity, Mars' velocity, and relative arrival velocity
% quiver3(r2(1), r2(2), r2(3), X_sc(end, 4) - X_M(end, 4), X_sc(end, 5) - X_M(end, 5), X_sc(end, 6) - X_M(end, 6), '-m');
% quiver3(r2(1), r2(2), r2(3), X_sc(end, 4), X_sc(end, 5), X_sc(end, 6), '-g');
% quiver3(r2(1), r2(2), r2(3), X_M(end, 4), X_M(end, 5), X_M(end, 6), '-r');
% hold off; % Release plot hold
