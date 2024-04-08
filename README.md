# Trajectory Design for Mars 2020

## Overview

This project presents a design for a Mars 2020 spacecraft trajectory, utilizing a simplified solar system model. By focusing on the heliocentric two-body orbits of Earth and Mars, described through Keplarian elements, it simplifies the interplanetary transfer. This approach employs the zero-radius Sphere of Influence (SOI) patched conics approximation, allowing for the analysis of spacecraft motion primarily influenced by the Sun's gravity, thereby simplifying the trajectory design between Earth and Mars. The use of the Heliocentric Ecliptic J2000 frame provides a consistent and accurate reference for all vectors and coordinate calculations.

## Methodology

The design and optimization of the spacecraft's trajectory to Mars are conducted through a series of calculated steps, using MATLAB for precise computation and iteration. 

1. **Orbit Determination**: Using MATLAB's `ode45` solver, calculate the initial and final orbits of the spacecraft relative to Earth and Mars. This forms the baseline trajectory.
2. **Trajectory Design**: The trajectory design process begins with the spacecraft's departure from Earth's orbit and its subsequent capture into Mars' orbit, calculated to ensure precise position and velocity for successful mission execution. Following this, a Hohmann transfer orbit is planned, characterized by its efficiency in terms of fuel consumption and time, making it an ideal choice for connecting concentric, coplanar orbits between Earth and Mars.
3. **Optimization**: Optimization of the spacecraft's trajectory focuses on three key areas: velocity matching, to ensure the spacecraft's speed aligns with the that required for Mars orbit capture; plane adjustment, to adjust the transfer orbit so it matches the orbital planes of both Earth and Mars, correcting any inclination discrepancies; and fuel efficiency, where the final trajectory adjustments are made to minimize the delta-v required for direction changes, critical for overall mission viability and success. 

## User Guide

### Simulation Setup and Execution

- **Prerequisites**: Install MATLAB and configure it with the Aerospace Toolbox. Clone the repository and ensure the `state_vector.m` file is included in the MATLAB working directory.
- **Running Simulations**: Execute `main.m` to start the trajectory design process.

### Utilities

- **State Vector and Orbital Element Conversion**: Utilize `rv_oe.m` and `oe_rv.m` for converting between state vectors and Keplarian orbital elements.
- **Hohmann Transfer Calculations**: The `hohmann_transfer.m` function calculates the velocities for a Hohmann transfer, providing initial delta-v estimates for the mission.
