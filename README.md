# PorousMedium1D
A code solving the 1D porous medium equation
$$
\begin{equation*}
    \frac{\partial u}{\partial t}\left(t, x\right) = \alpha\frac{\partial}{\partial x^2}\left(u^{\gamma + 1}\right)\left(t, x\right)\;,
\end{equation*}
$$
where $\alpha > 0$ and $\gamma > 0$, with Dirichlet boundary conditions
$$
\begin{align*}
    u\left(t, 0\right) &= u_0\\
    u\left(t, L\right) &= u_L
\end{align*}
$$
by implicit time integration (backward Euler).

**NOTE:** we must have $u\geq 0$ at all times for $u^\gamma$ to be well-defined
