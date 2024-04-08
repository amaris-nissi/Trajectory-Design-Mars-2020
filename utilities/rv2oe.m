% Clear workspace
close all; clear; clc;

% Define the initial state vector (position and velocity) of the object
rv = [-1.9717e8, -1.329e8, 2.0529e6, 1.4449e1, -1.8018e1, -7.3209e-1];
% Define the gravitational parameter of the primary body (Sun) in km^3/s^2
mu = 1.327e11;

% Call the function to convert the state vector to orbital elements
[oe] = rv2oel(rv, mu);

function [oe] = rv2oel(rv, mu)
    % Extract position and velocity vectors from the input state vector
    r_ijk = rv(1:3);
    v_ijk = rv(4:6);

    % Calculate the magnitudes of position and velocity vectors
    r = norm(r_ijk);
    v = norm(v_ijk);

    % Calculate specific angular momentum vector and its magnitude
    hvec = cross(r_ijk, v_ijk);
    h = norm(hvec);

    % Define the unit vector in the z-direction
    zhat = [0, 0, 1];

    % Calculate the node vector (nhat) - the ascending node line
    nhat = cross(zhat, hvec) / norm(cross(zhat, hvec));

    % Calculate the eccentricity vector and its magnitude
    evec = (cross(v_ijk, hvec) / mu) - (r_ijk / r);
    e = norm(evec);

    % Calculate the right ascension of the ascending node (Omega, Ω)
    om = atan2(nhat(2), nhat(1));

    % Calculate the inclination of the orbit (i)
    i = acos(dot(zhat, hvec) / h);

    % Calculate the argument of periapsis (omega, ω)
    w = acos(dot(nhat, evec) / e);
    % Adjust the argument of periapsis based on the direction of evec
    if evec(3) < 0
        w = -w;
    end

    % Calculate the true anomaly (ν)
    nu = acos(dot(r_ijk, evec) / (r * e));
    % Adjust the true anomaly based on the direction of motion
    if dot(r_ijk, v_ijk) < 0
        nu = -nu;
    end

    % Calculate the specific orbital energy
    energy = (v^2 / 2) - (mu / r);

    % Calculate the semi-major axis (a)
    a = -mu / (2 * energy);

    % Compile the computed orbital elements into the output vector
    oe = [a, e, i, w, om, nu];
end
