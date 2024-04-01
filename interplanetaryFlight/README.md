# Interplanetary flight

## Tsiolkowski equation:

```matlab
[mf, dm] = impulseFuel(m0, dv, Isp, g0)
```

Fuel consumption is calculated according to Tsiolkovsky equation.

$m_f$ represents final mass and $dm$ represents fuel consumption.

$m_0$ represents the initial mass, $dv$ represents the pulse, $I_{SP}$ represents the specific impulse, and $g_0$ represents the acceleration of gravity (typically $9.8m/s^2$).

## Gravity assist - equivalent pulse model

```matlab
[vMid1, vMid2, dvGA] = SOI_opt(v1, v2, vPlanet, muPlanet, rp, phi)
```

$v_{Mid_1}$ represents the remaining velocity of the influence sphere, $v_{Mid_2}$ represents the remaining velocity of the influence sphere, and $dv_{GA}$ represents the total pulsing force.

$v_1,v_2$represent the velocities of the incoming and outgoing spheres in the solar reference system;

$v_{Planet}$ represents the velocity of the Planet when it pulls, $\mu_{Planet}$ represents the gravitational coefficient of the planet;

$r_p$ indicates the borrowing altitude and $\phi$indicates the Angle of the orbital plane. The detailed model can be found in Chapter 3 of the Dynamics and Control of Deep Space Exploration.

## Objective function

'biGA_obj', 'impulse_obj', already covered in the root README.

These two objective functions have been processed into unconstrained objective functions, and the processing method is penalty function method, which is relatively simple, but the effect is not optimal. It is recommended to use the PSO optimizer provided by MATLAB for optimization. (You can also use fmincon for constrained optimization for better results)