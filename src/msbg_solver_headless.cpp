/******************************************************************************
 *
 * Headless driver for the reusable MSBG solver module.
 *
 * This tool performs the sparse point-cloud reconstruction without any UI so
 * that the pipeline can be validated in isolation (Step 1 of the integration
 * plan before wiring the solver into Houdini's HDK).
 *
 ******************************************************************************/

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <iostream>

#include "msbg_solver_module.h"
#include "pnoise.h"
#include "thread.h"
#include "util.h"

int nMaxThreads = -1;

namespace
{

void showUsage()
{
  std::cout << "Usage: msbg_solver_headless -f <points.ply> [options]\n"
            << "Options:\n"
            << "  -f <file>          Input point cloud in PLY format\n"
            << "  -b <size>          Base block size (16 or 32, default: 16)\n"
            << "  -r <Nx,Ny,Nz>      Effective resolution (default: 1024,1024,1024)\n"
            << "  -c <case>          Test case preset (default: 1)\n"
            << "  -i <iterations>    Mean curvature iterations (default: 5)\n"
            << "  -t <dt>            PDE time step override\n"
            << "  -o <dir>          Visualization output directory\n"
            << "  -N <instances>    Limit replicated base point instances (default: all)\n"
            << "  -P <points>       Limit number of source points sampled (default: all)\n"
            << "  -V                Enable visualization bitmaps\n"
            << "  -v <level>        Log level (0-3, default: 2)\n"
            << std::endl;
}

} // namespace

int main(int argc, char **argv)
{
  using namespace MSBG::Solver;

  SolverSettings settings;
  settings.blockSize = 16;
  settings.resolutionX = settings.resolutionY = settings.resolutionZ = 1024;
  settings.testCase = 1;
  settings.pdeIterations = 5;

  const char *inputFile = nullptr;

  int opt;
  while ((opt = getopt(argc, argv, "hf:b:r:c:i:t:o:N:P:Vv:")) != -1)
  {
    switch (opt)
    {
    case 'h':
      showUsage();
      return 0;
    case 'f':
      inputFile = optarg;
      break;
    case 'b':
      settings.blockSize = atoi(optarg);
      break;
    case 'r':
    {
      UT::StrTok t(optarg, ",");
      if (t.n() != 3)
      {
        TRCERR((("-r expects three comma separated values")));
        return 1;
      }
      settings.resolutionX = static_cast<int>(t.getFloat(0));
      settings.resolutionY = static_cast<int>(t.getFloat(1));
      settings.resolutionZ = static_cast<int>(t.getFloat(2));
      break;
    }
    case 'c':
      settings.testCase = atoi(optarg);
      break;
    case 'i':
      settings.pdeIterations = atoi(optarg);
      break;
    case 't':
      settings.timeStep = static_cast<float>(atof(optarg));
      break;
    case 'o':
      settings.outputDirectory = optarg;
      break;
    case 'N':
      settings.maxInstances = static_cast<uint64_t>(atoll(optarg));
      break;
    case 'P':
    {
      long long cap = atoll(optarg);
      settings.maxBasePoints = cap > 0 ? static_cast<size_t>(cap) : 0;
      break;
    }
    case 'V':
      settings.enableVisualizationOutput = true;
      break;
    case 'v':
      LogLevel = atoi(optarg);
      break;
    default:
      showUsage();
      return 1;
    }
  }

  if (!inputFile)
  {
    TRCERR((("No input point cloud specified")));
    showUsage();
    return 1;
  }

  settings.pointCloudFile = inputFile;

  ThrInit();

  int nMaxOmpThreadsAct = 0;
  int nMaxTbbThreadsAct = 0;
  ThrGetMaxNumberOfThreads(&nMaxOmpThreadsAct, &nMaxTbbThreadsAct);
  int nMaxThreadsAct = std::max(nMaxOmpThreadsAct, nMaxTbbThreadsAct);

  nMaxThreads = nMaxThreadsAct > 0 ? nMaxThreadsAct : 1;
  ThrSetMaxNumberOfTBBThreads(nMaxThreads);

  TRCP(("Using %d threads (OMP=%d, TBB=%d)\n", nMaxThreads, nMaxOmpThreadsAct,
        nMaxTbbThreadsAct));

  int noiseInit = PNS_Init(1234);
  if (noiseInit)
  {
    std::cerr << "PNS_Init failed with error code " << noiseInit << std::endl;
    return noiseInit;
  }

  char *logFile = (char *)"msbg_solver_headless.log";
  if (UtOpenLogFile(logFile, FALSE, FALSE))
  {
    std::cerr << "Failed to open log file: " << logFile << std::endl;
    return 1;
  }

  TRCP((("Running headless MSBG solver...\n")));

  int rc = RunPointCloudSolver(settings, {});
  if (rc != 0)
  {
    std::cerr << "RunPointCloudSolver failed with error code " << rc << std::endl;
    return rc;
  }

  TRCP((("Done.\n")));
  return 0;
}

