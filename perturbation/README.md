# Orbital perturbation - J2 perturbation

## Original two-body orbital motion equation (ODE)

```matlab
dRVdt = twoBodyOde(t, RV, mu)
```

## $J_2$ equation of two body orbital motion under perturbation

1. Average Method:

```matlab
[rt, vt] = rv02rvf_aveJ2(r0, v0, dt, mu, Re)
```

2. Exact $J_2$ perturbation (ODE)

```matlab
dRVdt = twoBodyJ2Ode(t, RV, mu, Re)
```

The ODE equations please use MATLAB to solve the ODE solver [ode45](https://www.mathworks.com/help/matlab/ref/ode45.html) (recommended).