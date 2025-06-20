# Cantilever Beam Finite Element Analysis (FEA) in MATLAB

## Project Overview

This project implements a Finite Element Analysis (FEA) framework in MATLAB for analyzing the deflection and rotation of a cantilever beam under a point load at its free end. Developed as part of a university assignment for "Finite Element Method for the Analysis of Structures" (FEA II) in an MEng Mechanical Engineering and Aeronautics program, the project focuses on comparing the behavior of beam elements derived from different structural theories: Euler-Bernoulli and Timoshenko (with both Full and Reduced Integration).

The analysis incorporates realistic material properties and a specific cross-sectional geometry, allowing for comprehensive insights into beam mechanics under various discretization and slenderness conditions.

## Features

*   **Material Properties:** User-defined Young's Modulus (E), Poisson's Ratio (ν), and derived Shear Modulus (G).
*   **Geometric Input:** User-specified beam length (L) and cross-sectional parameters.
*   **Cross-Sectional Properties:** Automated calculation of Area (A) and Area Moment of Inertia (I) for a **solid semicircle** cross-section.
*   **Beam Element Formulations:**
    *   **Euler-Bernoulli Beam Element:** Classical theory neglecting shear deformation.
    *   **Timoshenko Beam Element (Full Integration):** Accounts for shear deformation with standard numerical integration.
    *   **Timoshenko Beam Element (Reduced Integration):** Accounts for shear deformation using reduced integration, which can alleviate shear locking issues in slender beams.
*   **Global Stiffness Matrix Assembly:** Robust assembly of element stiffness matrices into the global system matrix.
*   **Boundary Conditions & Loading:** Imposition of fixed boundary conditions at the clamped end and application of a concentrated load at the free end.
*   **Solution:** Solves the system of linear equations to obtain nodal displacements and rotations.
*   **Parametric Studies:**
    *   **Element Discretization:** Analysis of deflection and rotation profiles for varying numbers of finite elements (e.g., 8, 24, 144 elements).
    *   **Beam Slenderness (L/h Ratio):** Investigation of tip deflection behavior as a function of the beam's length-to-height ratio.
*   **Visualization:** Generates informative plots comparing results from different beam theories, including:
    *   Transverse displacement along the beam length.
    *   Rotation along the beam length.
    *   Tip displacement vs. L/h ratio.
    *   Normalized tip displacement (w/w_Euler) vs. L/h ratio.

## Cross-Sectional Geometry (Solid Semicircle)

The project specifically handles a **solid semicircle** cross-section for its area (A) and moment of inertia (I). The `Area_and_Moment_of_Inertia.m` function calculates these properties.

**Note on Parameterization:** The current implementation of `Area_and_Moment_of_Inertia.m` uses parameters `h` and `a` to derive the semicircle's radius `r` (via `r = h*sin(a/2)`). While mathematically functional, this parametrization is non-standard for defining a simple semicircle (which is typically defined by its radius or diameter directly). If the original assignment intended a composite shape (e.g., the semicircle on top of a triangle as shown in Figure 2, Shape 8 of the assignment PDF), then the `Area_and_Moment_of_Inertia.m` function would need to be expanded to compute properties for such a composite section, including its centroid. For this project, it is assumed the analysis focuses solely on the solid semicircle geometry as implied by the specific image provided.

## Methodology

The core FEA workflow implemented follows standard procedures:

1.  **Element Stiffness Matrices:** For each chosen beam theory (Euler-Bernoulli, Timoshenko Full, Timoshenko Reduced), the `Local_Stifness.m` function calculates the 4x4 stiffness matrix for a single 2-node beam element (with 2 degrees of freedom: transverse displacement and rotation per node).
2.  **Global Assembly:** The `Global_Stifness.m` function systematically assembles these local stiffness matrices into a larger global stiffness matrix for the entire beam, respecting the nodal connectivity and degrees of freedom mapping.
3.  **Boundary Conditions & Load Vector:** The global stiffness matrix is condensed by applying the clamped boundary conditions (zero transverse displacement and rotation at the fixed end). The external point load is applied to the appropriate degree of freedom in the global force vector.
4.  **Solution:** The resulting system of linear algebraic equations (`K * U = F`) is solved using MATLAB's backslash operator (`\`) to obtain the unknown nodal displacements and rotations for the active (unconstrained) degrees of freedom.

## Usage

To run this FEA project:

1.  **Clone the repository** to your local machine.
2.  **Open MATLAB**.
3.  **Navigate to the `src/` directory** within the cloned repository.
4.  **Run the `main.m` script.** You can do this by typing `main` in the MATLAB Command Window and pressing Enter, or by clicking the "Run" button in the MATLAB editor.
5.  **Follow the prompts** in the Command Window to input the required beam parameters (length, cross-section parameters, number of elements, applied force).
6.  The script will then process the analysis and generate various plots illustrating the beam's behavior under different conditions.

## Code Quality & Possible Improvements

This project, typical of university assignments, prioritizes functional correctness. However, there's always room for enhancement!

*   **Input Validation:** Implement more robust input validation (e.g., using `validateattributes` or `inputParser`) to ensure user inputs are within reasonable and expected ranges.
*   **Function Encapsulation:** Consider refactoring the main script (`main.m`) into a function. This would improve reusability, make testing easier, and allow for cleaner execution with different sets of parameters without user interaction.
*   **Plotting Refinements:** While the `yyaxis` function effectively allows comparison of data on different scales, crafting a perfectly aesthetic and informative plot with multiple axes and legends across various scenarios can feel like a delightful challenge! Further customization of plot aesthetics (e.g., dynamic color schemes, more fine-grained axis control) could enhance clarity.
*   **Error Handling:** The existing code includes basic checks for `NaN`/`Inf` results from the solver, which is a good step towards robustness. More advanced error handling or warnings could be implemented for other potential issues (e.g., singular stiffness matrices, non-physical inputs).
*   **Comments in `Global_Stifness.m`:** While functional, adding more explicit comments within the `Global_Stifness.m` loop detailing the node-to-DOF mapping logic could further enhance readability for those less familiar with FEA assembly.

## Dependencies

*   **MATLAB:** This project requires a working installation of MATLAB.

## Author

Kampouropoulos Georgios
Student, MEng Mechanical Engineering and Aeronautics
University of Patras

## License

This project is open-sourced under the MIT License. See the `LICENSE` file for details.
