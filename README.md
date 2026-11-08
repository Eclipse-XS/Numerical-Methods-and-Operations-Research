# Numerical Methods and Operations Research
## Academic Projects | 2024–2025 Academic Year

### Overview
This repository contains academic software implementations developed during the 2024–2025 academic year, focusing on fundamental computational mathematics, numerical analysis, and mathematical optimization (Operations Research). 

The repository preserves multiple implementations, problem variations, and reference solutions for educational, archival, and comparative analysis purposes.

### Academic Context
These projects were developed during the 2024–2025 academic year.

### Numerical Methods (01_Numerical_Methods)
The Numerical Methods section covers algorithmic solutions to continuous mathematical problems:

1. **Nonlinear Equations and Root Finding (`01_Nonlinear_Equations/`):**
   - Fixed-point simple iteration method.
   - Newton-Raphson method with quadratic convergence.
   - Secant and chord methods eliminating analytical derivative requirements.
   - Analytical derivative extraction using symbolic computation.

2. **Polynomial Root Isolation (`02_Polynomial_Root_Isolation/`):**
   - Real root counting using Sturm chains and polynomial remainder sequences.
   - Sign change variation counting over bounded intervals.

3. **Linear Systems (`03_Linear_Systems/`):**
   - Direct solution of systems of linear algebraic equations using Gauss-Jordan elimination with partial pivoting.

4. **Nonlinear Systems (`04_Nonlinear_Systems/`):**
   - Multidimensional Newton-Raphson method for systems of nonlinear equations.
   - Dynamic analytical Jacobian matrix derivation and inversion.

5. **Numerical Integration (`05_Numerical_Integration/`):**
   - Composite quadrature approximations: Left, Right, and Midpoint Rectangle rules.
   - Composite Trapezoidal rule.
   - Composite Simpson's 1/3 rule with fourth-order convergence.

6. **Interpolation (`06_Interpolation/`):**
   - Natural cubic spline interpolation.
   - Continuous piecewise polynomial construction and boundary curvature constraints.

7. **Differential Equations (`07_Differential_Equations/`):**
   - Fourth-order classical Runge-Kutta (RK4) method for initial value problems (IVP).
   - Finite-difference central approximation method for second-order two-point boundary value problems (BVP).

### Operations Research (02_Operations_Research)
The Operations Research section encompasses discrete and continuous mathematical programming, linear optimization, network routing, and non-linear search algorithms:

1. **Linear Programming (`01_Linear_Programming/`):**
   - **Graphical Method:** Two-variable polytope intersection, feasible region isolation, and objective gradient tracing.
   - **Simplex Method:** Canonical tableau transformations, pivot selection using Dantzig's rule, and minimum ratio tests.
   - **Duality and Big-M Method:** Primal-to-dual transformation, dual simplex iteration, and artificial basis penalty formulation.

2. **Transportation Problems (`02_Transportation_Problem/`):**
   - Balanced and unbalanced transportation models with automated dummy node augmentation.
   - Initial basic feasible solution generation via the North-West Corner method.
   - Optimality verification and cycle re-distribution via the Method of Potentials ($u_i + v_j = c_{ij}$).

3. **Integer Programming (`03_Integer_Programming/`):**
   - Gomory's fractional cutting-plane algorithm.
   - Dual simplex re-optimization following cut addition.
   - High-precision symbolic rational arithmetic using exact fractions to prevent truncation errors.

4. **Fractional-Linear Programming (`04_Fractional_Programming/`):**
   - Algebraic linearization using the Charnes-Cooper transformation.
   - Geometric extreme-point evaluation in two-variable fractional problem spaces.

5. **Dynamic Programming (`05_Dynamic_Programming/`):**
   - Discrete multi-stage resource allocation governed by Bellman's principle of optimality.
   - Backward state-decision recurrence relations.

6. **Nonlinear Optimization (`06_Nonlinear_Optimization/`):**
   - **Nelder-Mead Method:** Zero-order geometric simplex search with reflection, expansion, and contraction.
   - **Fletcher-Reeves Method:** First-order conjugate gradient optimization with analytical gradient calculations.
   - **Davidon-Fletcher-Powell (DFP) Method:** Second-order quasi-Newton optimization with rank-two inverse Hessian updates.
   - **Parabolic Interpolation:** One-dimensional line search via three-point quadratic approximations.

### Repository Structure

```text
Numerical-Methods-and-Operations-Research/
│
├── 01_Numerical_Methods/
│   ├── 01_Nonlinear_Equations/         # Single-variable root finding
│   ├── 02_Polynomial_Root_Isolation/   # Sturm chain polynomial analysis
│   ├── 03_Linear_Systems/              # Gauss-Jordan elimination
│   ├── 04_Nonlinear_Systems/           # Multidimensional Newton-Raphson
│   ├── 05_Numerical_Integration/       # Composite numerical quadrature
│   ├── 06_Interpolation/               # Cubic spline interpolation
│   └── 07_Differential_Equations/      # RK4 integrator & finite-difference BVP
│
├── 02_Operations_Research/
│   ├── 01_Linear_Programming/
│   │   ├── Graphical_Method/           # 2D linear inequality solvers
│   │   ├── Simplex_Method/             # Standard simplex tableau solvers
│   │   └── Duality_and_Big_M/          # Duality and artificial basis methods
│   ├── 02_Transportation_Problem/      # Potentials method and distribution models
│   ├── 03_Integer_Programming/         # Gomory fractional cutting-plane methods
│   ├── 04_Fractional_Programming/      # Charnes-Cooper algebraic & graphical solvers
│   ├── 05_Dynamic_Programming/         # Bellman multi-stage capital allocation
│   └── 06_Nonlinear_Optimization/      # Nelder-Mead, Fletcher-Reeves, DFP, Parabolic
│
├── requirements.txt                    # Python runtime package dependencies
├── README.md                           # Repository documentation
└── .gitignore                          # Visual Studio and Python build exclusions
```

### Technologies
The implementations in this repository utilize:
- **Python 3:** Primary language for mathematical models and scripting.
- **C++ (MSVC / ISO C++17):** High-performance procedural implementations of linear programming solvers (`do`, `DO_BIG-M`).
- **SymPy:** Computer algebra system for analytical derivatives, exact rational arithmetic, and symbolic Jacobian construction.
- **NumPy & SciPy:** Dense matrix linear algebra, vectorization, and spline interpolation routines.
- **Matplotlib:** Graphical visualization of feasible regions, objective function contours, and differential solution curves.
- **Tabulate:** Formatted ASCII/Markdown console reporting for simplex and quadrature steps.
- **python-docx:** Automated report generation and table export.

### Project Organization
Multiple implementations of individual algorithms exist across subdirectories. This architecture intentionally preserves:
- Canonical reference solutions.
- Extended calculators supporting generalized matrix dimensions.
- Historical variations demonstrating differing input parameters, constraint structures, or output formats.
- Language-specific implementations comparing Python symbolic workflows against C++ procedural memory models.

Detailed cluster documentation is provided in the `README.md` file located inside each respective subfolder.

### Running Projects

#### Python Solutions
1. Ensure Python 3.10 or later is installed.
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Navigate to the project directory and execute the primary script:
   ```bash
   python <ScriptName>.py
   ```

#### Visual Studio Solutions
1. Open the `.sln` file in Visual Studio 2022 (with Python Development or Desktop Development with C++ workloads installed).
2. Set the startup project if prompted.
3. Build and run (Ctrl+F5).

### Limitations
These implementations were authored for academic demonstration and computational research. They are not intended as hardened production solvers:
- Floating-point implementations without symbolic rational arithmetic may experience rounding drift on ill-conditioned matrices.
- Command-line interfaces generally rely on pre-configured input matrices within the script or formatted console prompts rather than external file streaming.
- Error handling prioritizes mathematical transparency over comprehensive input fault tolerance.
