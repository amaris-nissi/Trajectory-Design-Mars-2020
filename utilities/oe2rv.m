% Clear workspace
close all; clear; clc;

function [rv] = orb2rv(oe, mu)
% Converts orbital elements to position and velocity vectors in the inertial frame.
%
% Inputs:
%   oe - Array of orbital elements [a, e, i, Omega, omega, theta]
%        a     - Semi-major axis
%        e     - Eccentricity
%        i     - Inclination (degrees)
%        Omega - Right ascension of the ascending node (RAAN) (degrees)
%        omega - Argument of periapsis (degrees)
%        theta - True anomaly (degrees)
%   mu - Gravitational parameter of the central body
%
% Outputs:
%   rv - Position and velocity vectors in the inertial frame

% Calculate position vector in the perifocal coordinate system
r = (oe(1) * (1 - oe(2).^2)) / (1 + oe(2) * cos(deg2rad(oe(6))));
phat = [1; 0; 0]; % Unit vector along perifocal coordinate system's p-axis
qhat = [0; 1; 0]; % Unit vector along perifocal coordinate system's q-axis

% Position vector in PQW frame
r_pqw = (r * cos(deg2rad(oe(6))) * phat) + (r * sin(deg2rad(oe(6))) * qhat);

% Calculate velocity vector in the perifocal coordinate system
p = oe(1) * (1 - oe(2)^2); % Semi-latus rectum
v_pqw = (-sin(deg2rad(oe(6))) * sqrt(mu/p) * phat) + (sqrt(mu/p) * (oe(2) + cos(deg2rad(oe(6)))) * qhat);

% Rotation matrix for Omega, RAAN
R3w = [cos(deg2rad(oe(4))), -sin(deg2rad(oe(4))), 0;
       sin(deg2rad(oe(4))), cos(deg2rad(oe(4))), 0;
       0, 0, 1];

% Rotation matrix for i, inclination
R1i = [1, 0, 0;
       0, cos(deg2rad(oe(3))), -sin(deg2rad(oe(3)));
       0, sin(deg2rad(oe(3))), cos(deg2rad(oe(3)))];

% Rotation matrix for omega, argument of periapsis
R3om = [cos(deg2rad(oe(5))), -sin(deg2rad(oe(5))), 0;
        sin(deg2rad(oe(5))), cos(deg2rad(oe(5))), 0;
        0, 0, 1];

% Composite rotation matrix to convert from PQW to IJK (inertial) frame
R = R3w * R1i * R3om;

% Position and velocity in the inertial frame
r_ijk = R * r_pqw;
v_ijk = R * v_pqw;

% Output combined position and velocity vector
rv = [r_ijk; v_ijk];

end
