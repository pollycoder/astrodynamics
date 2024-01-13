# Tsinghua University space dynamics course learning warehouse

## Content
1. Two body problem: the conversion of six radical number of orbit and Descartes position velocity, the conversion and evolution of three near-point angles, the initial value of orbit, Newton iteration.
Note: Support ellipse, hyperbola, parabola three cases.
2. Orbital maneuver: Hohmann transfer, double ellipse triple pulse transfer, including coplanar and non-coplanar cases, the advantages and disadvantages of the two transfer methods under different working conditions.
3. Orbital perturbation: J2 perturbation average method, spherical harmonic function method of Earth non-spherical gravitational perturbation, J2 perturbation of two-body problem differential equation solution. Other perturbative processes will be added later.
4. Coordinate and time systems: calendar and Earth coordinate system conversion.
Interplanetary flight: Trajectory optimization of double pulse transfer in arbitrary orbits, influencing spherical equivalent pulse model.
5. Detailed information is provided in the READEME file of each part.


## Requirements
1. The main program is completed by MATLAB
2. The coordinate and time system part is written in c++ and compiled using cmake. Please install The toolkit SOFA released by IAU and compile the static link library in advanced according to the official documentation.
3. Please run the code in computer with x86 architecture because the static link library is not compatible with arm64 architecture.
