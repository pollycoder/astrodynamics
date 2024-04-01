# 轨道摄动 - J2摄动

## 原始二体轨道运动方程(ODE)

```matlab
dRVdt = twoBodyOde(t, RV, mu)
```

## $J_2$摄动下的二体轨道运动方程

1. 平均法：

```matlab
[rt, vt] = rv02rvf_aveJ2(r0, v0, dt, mu, Re)
```

2. 精确$J_2$摄动（ODE）

```matlab
dRVdt = twoBodyJ2Ode(t, RV, mu, Re)
```

ODE方程请使用MATLAB的ODE求解器进行求解（推荐[ode45](https://www.mathworks.com/help/matlab/ref/ode45.html)）。