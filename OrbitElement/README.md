# 轨道根数计算

1. 偏近点角$\iff$平近点角

```matlab
M = E2M(E, e)
E = M2E(M, e)
```

$M$表示平近点角，$E$表示偏近点角，$e$表示偏心率。

2. 偏近点角$\iff$真近点角

```matlab
f = E2f(E, e)
E = f2E(f, e)
```

$f$表示真近点角，$E$表示偏近点角，$e$表示偏心率。

3. 真近点角的演化

```matlab
ft = f0dt2ft(f0, dt, a, e, mu)
dt = f0ft2dt(f0, ft, a, e, mu)
```

$f0$表示初始真近点角，$ft$表示终末真近点角，$a$表示半长轴，$e$表示偏心率，$\mu$表示中心天体引力系数。

4. 轨道根数和笛卡尔坐标系的转换&演化

```matlab
coe = rv2coe(CartesianR, CartesianV, mu)
[CartesianR, CartesianV] = coe2rv(coe, mu, tol)
```

coe表示轨道根数数组（$len=6$），顺序为$a, e, i, \Omega, \omega, f$

CartesianR，CartesianV分别表示笛卡尔坐标系下的位置和速度；

$\mu$表示中心天体引力系数，$tol$表示奇异性容忍度。

5. 二体轨道初值问题求解器

```
[rt, vt] = rv02rvf(r0, v0, dt, mu)
```

$r_t,v_t$表示终末位置速度，$r_0,v_0$表示初始位置速度，dt表示运动时间，$\mu$表示中心天体引力系数。

6. 牛顿迭代法

```
[x,steps]=newton(obj_fun,x0,tol,max_iter)
```

详情可查询数值分析学习库：[Newton Iteration](https://github.com/pollycoder/numerical_analysis/tree/main/exp2_nonLinear)