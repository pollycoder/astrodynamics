# Calculation of track roots

1. Near Angle $\iff$Flat near Angle

```matlab
M = E2M(E, e)
E = M2E(M, e)
```

$M$ indicates the flat approach Angle, $E$ indicates the off-approach Angle, and $e$ indicates the eccentricity.

2. Close Angle $\iff$True close Angle

```matlab
f = E2f(E, e)
E = f2E(f, e)
```

$f$ represents the true near-point Angle, $E$ represents the partial near-point Angle, and $e$ represents the eccentricity.

3. Evolution of true peripoint Angle

```matlab
ft = f0dt2ft(f0, dt, a, e, mu)
dt = f0ft2dt(f0, ft, a, e, mu)
```

$f_0$ represents the initial true near-point Angle, $f_t$ represents the final true near-point Angle, $a$represents the semi-major axis, $e$ represents the eccentricity, and $\mu$represents the gravitational coefficient of the central body.

4. Transformation & evolution of orbital roots and Cartesian coordinate system

```matlab
coe = rv2coe(CartesianR, CartesianV, mu)
[CartesianR, CartesianV] = coe2rv(coe, mu, tol)
```

coe represents the group of orbital roots ($len=6$) in the order $a, e, i, \Omega, \omega, f$

CartesianR and CartesianV represent the position and velocity in Cartesian coordinates respectively;

$\mu$ represents the gravitational coefficient of the central object, and $tol$represents the singularity tolerance.

5. Two-body orbit initial value problem solver

```matlab
[rt, vt] = rv02rvf(r0, v0, dt, mu)
```

$r_t,v_t$ represents the final position velocity, $r_0,v_0$ represents the initial position velocity, dt represents the motion time, and $\mu$ represents the gravitational coefficient of the central body.

6. Newton iteration method

```matlab
[x,steps]=newton(obj_fun,x0,tol,max_iter)
```

Details can inquire numerical analysis learning library: [Newton Iteration] (https://github.com/pollycoder/numerical_analysis/tree/main/exp2_nonLinear)