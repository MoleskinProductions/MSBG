# MSBG Houdini HDK Integration

This folder contains the source code for the `msbg_solver` SOP node. The node
wraps the reusable solver interface introduced in `src/msbg_solver_module.*` so
that Houdini can invoke the same sparse reconstruction logic that is exercised
by the headless test harness.

## Building

1. Install the Houdini Development Kit (HDK) that matches your Houdini build.
2. Open an HDK-enabled shell (e.g. `houdini_setup` on Linux) so that the `hcustom`
   helper is available.
3. Compile the plugin:

   ```bash
   hcustom -I.. -L.. -lmsbg SOP_MSBGSolver.C
   ```

   The command assumes the MSBG static library has been built via `../mk`.

4. Copy the resulting DSOs into your Houdini `dso` search path and restart
   Houdini. The node will appear as **MSBG Solver** in the SOP tab menu.

## Parameters

The node exposes the same tunable values as the headless executable:

| Parameter | Description |
| --- | --- |
| Block Size | MSBG base block size (16 or 32). |
| Resolution | Effective grid resolution in voxels. |
| PDE Iterations | Number of mean-curvature smoothing steps. |
| PDE Time Step | Optional override for the PDE solver time step. |
| Test Case | Matches the presets from the original `msbg_demo`. |
| Visualization Output | Destination directory for bitmap dumps. |
| Dump Visualizations | Enable/disable bitmap output. |

On cook the node copies the incoming point geometry, launches the solver via
`MSBG::Solver::RunPointCloudSolver`, and stores the execution status in a
detail-string attribute called `msbg_status`.

