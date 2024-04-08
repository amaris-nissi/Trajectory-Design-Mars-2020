function r_dot = state_vector(r, mu)
    % Calculates the derivative of the state vector in a two-body orbital dynamics problem.
    %
    % Parameters:
    % r : A 6-element vector representing the state of the spacecraft.
    %     r(1), r(2), r(3) are the spacecraft's Cartesian coordinates (x, y, z) in km.
    %     r(4), r(5), r(6) are the velocities in each direction (vx, vy, vz) in km/s.
    % mu: Gravitational parameter of the central body (Sun, Earth, etc.) in km^3/s^2.
    %
    % Returns:
    % r_dot: A 6-element vector representing the time derivative of the state vector.
    %        The first three elements are the velocities in the x, y, z directions,
    %        and the last three are the accelerations.

    % Calculate the radial distance from the central body
    r_ij = sqrt(r(1)^2 + r(2)^2 + r(3)^2);

    % Initialize the derivative vector
    r_dot = zeros(6, 1);

    % The time derivative of the position vector is the velocity vector
    r_dot(1) = r(4); % dx/dt = vx
    r_dot(2) = r(5); % dy/dt = vy
    r_dot(3) = r(6); % dz/dt = vz

    % Calculate the acceleration components using Newton's law of gravitation
    % a = -mu * r / |r|^3
    r_dot(4) = -(mu / (r_ij^3)) * r(1); % d^2x/dt^2 = ax
    r_dot(5) = -(mu / (r_ij^3)) * r(2); % d^2y/dt^2 = ay
    r_dot(6) = -(mu / (r_ij^3)) * r(3); % d^2z/dt^2 = az
end
