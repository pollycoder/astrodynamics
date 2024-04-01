# Astrodynamics - Code Material Repository

The repository is for Astrodynamics course in Tsinghua University, it is designed to help learners of this course to do experiements with what is taught in course through visualization and examples, so as to have a more intuitive understanding of the course.

There are READMEs under each directory, please read them before running the code.

If there are any problems, please raise your questions in Issue column.

Or contact me through email: pollyjoe2003@gmail.com

## Environment requirements

Please make sure that [Global Optimization Toolbox](https://www.mathworks.com/products/global-optimization.html) is installed on your MATLAB, tutorials could be found on the official website of MATLAB.

## Structure：

```shell
.
├── LICENSE
├── OrbitElement
├── README.md
├── interplanetaryFlight
├── maneuvers
├── perturbation
├── plot
├── relativeMotion
├── summary.pdf
└── test_examples
```

| File Name             | Content                                                      |
| --------------------- | ------------------------------------------------------------ |
| /OrbitElement         | The transformation of six radical number of orbit and Cartesian position velocity, the transformation and evolution of three near-point angles, the initial value of orbit, Newton iteration. <br /> **Note: Supports ellipse, hyperbola, and parabola.** |
| /maneuvers            | Orbit maneuver: Hohmann transfer, double ellipse triple pulse transfer, compatible with coplanar and non-coplanar cases, and the advantages and disadvantages of two transfer methods under different working conditions (orbit radius ratio $r_2/r_1$, inclination difference $i_2-i_1$) are compared. |
| /relativeMotion       | Relative motion: Comparison of the trajectory and the $X,Y,Z$coordinates with time under periodic and general conditions |
| /perturbation         | Orbital perturbation: average method of $J_2$ perturbation, solution of two-body differential equation under $J_2$ perturbation. Other perturbative processes will be added later. |
| /interplanetaryFlight | Interplanetary flight: Influence spherical equivalent pulse model, Tsiolkowski equation solving fuel consumption, double pulse trajectory optimization, pulse trajectory optimization with gravity assistance. |
| /plot                 | Plotting module                                              |
| /test_examples        | Examples                                                     |
| summary.pdf           | Summary of the course                                        |

## Introduction of Examples

Content:

```shell
.
├── biE_test.m
├── biGA_test.m
├── hohmannVSbiE.m
├── hohmann_test.m
├── impulse_test.m
├── perturbation_test.m
└── relative_test.m
```

| Name              | Content                                                      |
| ----------------- | ------------------------------------------------------------ |
| hohmann_test      | Hohmann transfer, compatible with both coplanar and non-coplanar cases. Parameters $r_1,r_2,i_1,i_2$can be modified to observe the shape of the transfer orbit. |
| biE_test          | Double elliptic triple pulse transfer, compatible with coplanar and non-coplanar cases. Parameters $r_1,r_2,i_1,i_2$can be modified to observe the shape of the transfer orbit. |
| hohmannVSbiE      | Comparison of Hohmann transfer and double elliptic triple pulse transfer under different **orbital radius ratio/orbital inclination difference**. When the ratio of orbit radius is too large or the difference of inclination is too large, it is concluded that the double ellipse triple pulse transfer is better. |
| perturbation_test | The average method of $J_2$ perturbation, the comparison of trajectories after solving the differential equations of the two-body problem under $J_2$ perturbation, and the change of angular momentum. It should be observed that the surface swept by the angular momentum under the perturbation is a cone, corresponding to the phenomenon of the ascending node regression. <br/> In addition, we can also observe the end point position of the average method and the solution result of the differential equation, from which we can find that only the true near-point Angle is a fast variable in the orbital root number, and the rest are slow variables. |
| relative_test     | The comparison of the trajectories of C-W equation under the condition of periodic solution with the general case, as well as the variation of $X,Y and Z$over time, the phenomena should be consistent with what we have said in class. |
| impulse_test      | The basic method of trajectory optimization is demonstrated with the example of dual-pulse ground fire transfer. The corresponding objective function location: `. / interplanetaryFlight impulse_obj. M `, can be observed that the optimal trajectory and horman transfer trajectory. |
| biGA_test         | The trajectory optimization method with planetary force under the equivalent pulse model is demonstrated with the example of Mars double force mining. The corresponding objective function location: `. / interplanetaryFlight biGA_obj. M ` |

## Acknowledgement

Special thanks to our teacher Professor Fanghua Jiang from Department of Aeronautics and Astronautics, Tsinghua University, our TA Zipeng Wu , and Zhong Zhang from Tsinghua Aerospace Dynamics and Control Laboratory. Prof. Jiang and Zi Peng Wu gave me a lot of help and guidance in the course learning process, and Zhong Zhang provided a lot of useful suggestions for the design of the library. At the same time, some of the programs in this learning program library are also from their hands, among which the Lambert solver program of the interstellar flight part is provided by Prof. Jiang.

