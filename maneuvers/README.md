# Maneuver - Special transfer orbit

## Homan transfer

```matlab
[dv1, dv2, dv, dt, coeh] = hohmann(r1, r2, i1, i2, mu)
```

$dv_1$ and $dv_2$ represent the departure and arrival pulses respectively, $dv$ represents the characteristic velocity;

$dt$ represents the transfer time;

$coeh$ represents the three roots of the transfer orbit: $a,e,i$; 

$r_1, r_2$ represent the radius of initial orbit and final orbit (both circular orbits); 

$i_1,i_2$ represent the inclination of initial orbit and final orbit, respectively. $\mu$indicates the gravitational coefficient of the central object.

## Double elliptic triple pulse transfer

```matlab
[DV, dv, dt, at, et, it] = double_ellipse(r1, r2, i1, i2, mu)
```

$DV(len=3)$represents the size of the starting, middle and arriving pulses respectively.

$dv$ represents the characteristic speed and $dt$represents the transfer time.

$a_t,e_t$, and $i_t$ (all with $len=2$) represent the semi-major axis, eccentricity, and orbital inclination of the two transfer orbits, respectively.

$\textcolor{red}{Note:~orbital~inclination~adds~a~step~to~the~ optimization~ process~ when~ non-coplanar. } $