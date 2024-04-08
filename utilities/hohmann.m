% Clear workspace
close all; clear; clc;

% Define initial positions of Earth and Mars in their orbits at specific events
rE = [91930000, -120900000, 5318]; % Earth's position vector at launch (in km)
rM = [-1583000, 234800000, 4960000]; % Mars' position vector at arrival (in km)

% Calculate the norm of position vectors to determine orbital radii
r1 = norm(rE); % Orbital radius of Earth at launch (in km)
r2 = norm(rM); % Orbital radius of Mars at arrival (in km)

% Gravitational parameters (in km^3/s^2)
mu = 1.327e11; % Gravitational parameter of the Sun
mu_E = 3.986e5; % Gravitational parameter of Earth
mu_M = 4.2828e4; % Gravitational parameter of Mars

% Radii of planets (in km)
rp_E = 6378; % Radius of Earth
rc_E = rp_E; % Circular orbit assumption, radius equals planet's radius
rp_M = 3396; % Radius of Mars
rc_M = rp_M; % Circular orbit assumption, radius equals planet's radius

% Calculate parameters for Hohmann transfer orbit
[dv, tof, dva, dvp] = hohmann_trans(r1, r2, mu);

% Calculate velocities for escape from Earth
[vp_E, vc_E] = periapse_velocity(dvp, rp_E, rc_E, mu_E);
dv1 = vp_E - vc_E; % Delta-V required for escaping Earth

% Calculate velocities for capture at Mars
[vp_M, vc_M] = periapse_velocity(dva, rp_M, rc_M, mu_M);
dv2 = vp_M - vc_M; % Delta-V required for capture at Mars

% Total Delta-V for the mission
dvtotal = dv1 + dv2;

% Function to calculate parameters of a Hohmann transfer orbit
function [dv, tof, dva, dvp] = hohmann_trans(r1, r2, mu)
    % Calculate tangential velocities at departure and arrival
    vt1 = sqrt(((-2 * mu) / (r1 + r2)) + (2 * mu / r1));
    vt2 = sqrt(((-2 * mu) / (r1 + r2)) + (2 * mu / r2));
    
    % Calculate circular orbit velocities at departure and arrival
    vc1 = sqrt(mu / r1);
    vc2 = sqrt(mu / r2);
    
    % Delta-V calculations for the transfer
    dv = abs(vt1 - vc1) + abs(vt2 - vc2);
    dvp = abs(vt1 - vc1); % Delta-V at departure (Earth)
    dva = abs(vt2 - vc2); % Delta-V at arrival (Mars)
    
    % Time of flight for the Hohmann transfer (in seconds)
    tof = pi * sqrt(((r1 + r2)^3) / (8 * mu));
end

% Function to calculate periapse velocities
function [vp, vc] = periapse_velocity(vinf, rp, rc, mu)
    % Calculate periapse velocity
    vp = sqrt(vinf^2 + (2 * (mu / rp)));
    
    % Calculate circular orbit velocity
    vc = sqrt(mu / rc);
end
