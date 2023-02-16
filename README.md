# Runge–Kutta method RK4 for second order differential equations

Differential equations of the type:
y'' = f(x, y, y')

Due to the implementation with pointers to arrays, all of the variables are stored in arrays.

Console Output:
~~~text
------------------------------------------
x: 0.000, y: 0.000, yd: 4.000, ydd: -8.000
------------------------------------------
x: 0.020, y: 0.078, yd: 3.846, ydd: -8.000
------------------------------------------
x: 0.040, y: 0.154, yd: 3.703, ydd: -7.416
------------------------------------------
x: 0.060, y: 0.227, yd: 3.571, ydd: -6.865
------------------------------------------
x: 0.080, y: 0.297, yd: 3.449, ydd: -6.342
------------------------------------------
x: 0.100, y: 0.365, yd: 3.337, ydd: -5.848
------------------------------------------
x: 0.120, y: 0.430, yd: 3.234, ydd: -5.380
------------------------------------------
x: 0.140, y: 0.494, yd: 3.139, ydd: -4.937
------------------------------------------
x: 0.160, y: 0.556, yd: 3.053, ydd: -4.517
------------------------------------------
~~~

The code is ready to run. You can modify the function inside the first function

References:\
[1] Papula, Lothar. (2015, p. 485f). Mathematik für Ingenieure und Naturwissenschaftler - Band 2. 14th edition.
