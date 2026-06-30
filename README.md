# PorousMedium1D
A code solving the 1D porous medium equation with Dirichlet boundary conditions on a uniform grid by implicit time integration (backward Euler)


## Description
The porous medium equation in one spatial dimension with Dirichlet boundary conditions is

$$
\begin{align*}
    \frac{\partial u}{\partial t}\left(t, x\right) &= \alpha\frac{\partial}{\partial x^2}\left(u^{\gamma + 1}\right)\left(t, x\right) \\
    u\left(t, 0\right) &= u_0 \\
    u\left(t, L\right) &= u_L
\end{align*}
$$

where $\alpha > 0$ and $\gamma > 0$. This code solves this equation by implicit time integration (backward Euler) on a uniform grid.

**NOTE:** we must have $u\geq 0$ at all times for $u^\gamma$ to be well-defined.


## Minimal requirements
- A C++ compiler supporting the `c++17` standard
- `cmake`
- The HDF5 library
- `python3` with `numpy`, `matplotlib`, and `h5py` to generate evolution snapshots
- `ffmpeg` to generate movies


## Usage
1. Tune the parameters in `Parameters.hh`
2. Compile:
   ```
   ./build.sh
   ```
3. Run:
   ```
   ./install/bin/PorousMedium1D_exe
   ```
4. To plot evolution snapshots and make movies, see `Utils/README.md`
