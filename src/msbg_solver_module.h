#pragma once

#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "VEC3.h"
#include "msbg.h"

namespace MSBG
{
namespace Solver
{

struct SolverSettings
{
  std::string pointCloudFile;          ///< Optional path to PLY file containing the source point cloud.
  std::vector<Vec3Float> pointCloud;   ///< Optional in-memory copy of the source points (copied before use).

  int blockSize = 16;                  ///< Base block size (16 or 32).
  int resolutionX = 1024;              ///< Effective domain resolution in X.
  int resolutionY = 1024;              ///< Effective domain resolution in Y.
  int resolutionZ = 1024;              ///< Effective domain resolution in Z.

  int testCase = 1;                    ///< Keeps compatibility with existing presets from msbg_demo.cpp.
  float particleRadius = 2.0f;         ///< Particle splat radius used for point rasterisation.
  float neighborDistance = 2.0f;       ///< Additional radius for neighbor search when activating blocks.
  int pdeIterations = 40;              ///< Number of mean curvature smoothing iterations.
  float timeStep = 0.05f;              ///< Time step for the PDE solver.

  std::string outputDirectory = "out_msbg_solver"; ///< Directory used for optional visualization dumps.
  bool enableVisualizationOutput = false;           ///< Enable bitmap dumps via MultiresSparseGrid.
  uint64_t maxInstances = 0;            ///< Optional cap on the number of replicated instances.
  size_t maxBasePoints = 0;             ///< Optional limit for the number of source samples used.
};

struct SolverCallbacks
{
  std::function<void(MultiresSparseGrid*, int)> afterParticleSplat; ///< Called after particles have been rasterised.
  std::function<void(MultiresSparseGrid*, int)> afterPdeSolve;       ///< Called once the PDE solve finished.
};

/**
 * Execute the MSBG sparse point-cloud reconstruction pipeline.
 *
 * The solver mirrors the behaviour of the original msbg_demo test case but is
 * exposed via a reusable interface that works in both headless tools and
 * DCC integrations such as Houdini's HDK.
 *
 * @returns 0 on success, negative values on failure.
 */
int RunPointCloudSolver(const SolverSettings &settings,
                        const SolverCallbacks &callbacks = {});

} // namespace Solver
} // namespace MSBG

