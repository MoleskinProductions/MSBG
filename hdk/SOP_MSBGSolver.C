/******************************************************************************
 *
 * Houdini HDK SOP node that wraps the reusable MSBG solver module.
 *
 * The node expects a point-cloud input and runs the headless solver in-place,
 * allowing Houdini artists to trigger the MSBG reconstruction pipeline inside
 * the DCC while keeping the heavy lifting inside the standalone module.
 *
 ******************************************************************************/

#include <GA/GA_Handle.h>
#include <GA/GA_Iterator.h>
#include <GA/GA_Names.h>
#include <GU/GU_Detail.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Default.h>
#include <PRM/PRM_Template.h>
#include <SOP/SOP_Node.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_StringHolder.h>
#include <vector>

#include "msbg_solver_module.h"

namespace
{

constexpr const char *kParmBlockSize = "blocksize";
constexpr const char *kParmResolution = "resolution";
constexpr const char *kParmIterations = "iterations";
constexpr const char *kParmTimeStep = "timestep";
constexpr const char *kParmTestCase = "testcase";
constexpr const char *kParmOutputDir = "outputdir";
constexpr const char *kParmEnableVis = "enablevis";

PRM_Template *buildTemplates()
{
  static PRM_Name blocksizeName(kParmBlockSize, "Block Size");
  static PRM_Name resolutionName(kParmResolution, "Resolution");
  static PRM_Name iterationsName(kParmIterations, "PDE Iterations");
  static PRM_Name timestepName(kParmTimeStep, "PDE Time Step");
  static PRM_Name testcaseName(kParmTestCase, "Test Case");
  static PRM_Name outputdirName(kParmOutputDir, "Visualization Output");
  static PRM_Name enablevisName(kParmEnableVis, "Dump Visualizations");

  static PRM_Default blocksizeDefault(16);
  static PRM_Default iterationsDefault(5);
  static PRM_Default timestepDefault(0.0f);
  static PRM_Default testcaseDefault(1);
  static PRM_Default resolutionDefault[3] = {
      PRM_Default(1024), PRM_Default(1024), PRM_Default(1024)};

  static PRM_Template theTemplates[] = {
      PRM_Template(PRM_INT, 1, &blocksizeName, &blocksizeDefault),
      PRM_Template(PRM_INT, 3, &resolutionName, resolutionDefault),
      PRM_Template(PRM_INT, 1, &iterationsName, &iterationsDefault),
      PRM_Template(PRM_FLT, 1, &timestepName, &timestepDefault),
      PRM_Template(PRM_INT, 1, &testcaseName, &testcaseDefault),
      PRM_Template(PRM_STRING, 1, &outputdirName, nullptr),
      PRM_Template(PRM_TOGGLE, 1, &enablevisName, PRMoneDefaults),
      PRM_Template()
  };

  return theTemplates;
}

class SOP_MSBGSolver : public SOP_Node
{
public:
  static OP_Node *constructor(OP_Network *net, const char *name, OP_Operator *op)
  {
    return new SOP_MSBGSolver(net, name, op);
  }

  static PRM_Template *buildTemplateList() { return buildTemplates(); }

protected:
  SOP_MSBGSolver(OP_Network *net, const char *name, OP_Operator *op)
      : SOP_Node(net, name, op)
  {
  }

  ~SOP_MSBGSolver() override = default;

  OP_ERROR cookMySop(OP_Context &context) override
  {
    duplicateSource(0, context);

    fpreal t = context.getTime();
    int blockSize = evalInt(kParmBlockSize, 0, t);
    int resolutionX = evalInt(kParmResolution, 0, t);
    int resolutionY = evalInt(kParmResolution, 1, t);
    int resolutionZ = evalInt(kParmResolution, 2, t);
    int iterations = evalInt(kParmIterations, 0, t);
    int testCase = evalInt(kParmTestCase, 0, t);
    fpreal timeStep = evalFloat(kParmTimeStep, 0, t);
    bool enableVis = evalInt(kParmEnableVis, 0, t) != 0;

    UT_StringHolder outputDir;
    evalString(outputDir, kParmOutputDir, 0, t);

    GU_DetailHandleAutoReadLock inputLock(inputGeoHandle(0));
    const GU_Detail *inputGeo = inputLock.getGdp();
    if (!inputGeo)
    {
      addWarning(SOP_MESSAGE, "No input geometry to process.");
      return UT_ERROR_NONE;
    }

    GA_ROHandleV3 hP(inputGeo, GA_ATTRIB_POINT, GA_Names::P);
    if (!hP.isValid())
    {
      addError(SOP_MESSAGE, "Input geometry has no position (P) attribute.");
      return UT_ERROR_ABORT;
    }

    std::vector<Vec3Float> pointCloud;
    pointCloud.reserve(static_cast<size_t>(inputGeo->getNumPoints()));
    for (GA_Iterator it(inputGeo->getPointRange()); !it.atEnd(); ++it)
    {
      const GA_Offset offset = it.getOffset();
      const UT_Vector3 pos = hP.get(offset);
      pointCloud.emplace_back(static_cast<float>(pos.x()),
                              static_cast<float>(pos.y()),
                              static_cast<float>(pos.z()));
    }

    MSBG::Solver::SolverSettings settings;
    settings.pointCloud = pointCloud;
    settings.blockSize = blockSize;
    settings.resolutionX = resolutionX;
    settings.resolutionY = resolutionY;
    settings.resolutionZ = resolutionZ;
    settings.testCase = testCase;
    settings.pdeIterations = iterations > 0 ? iterations : 5;
    settings.timeStep = static_cast<float>(timeStep);
    settings.enableVisualizationOutput = enableVis;
    if (enableVis && outputDir.isstring())
      settings.outputDirectory = outputDir.c_str();

    MSBG::Solver::SolverCallbacks callbacks;
    callbacks.afterPdeSolve = [this](MSBG::MultiresSparseGrid *, int) {
      addMessage(SOP_MESSAGE, "MSBG solver finished PDE stage.");
    };

    int rc = MSBG::Solver::RunPointCloudSolver(settings, callbacks);
    GA_RWHandleS statusHandle(
        gdp->addStringTuple(GA_ATTRIB_DETAIL, "msbg_status", 1));
    if (statusHandle.isValid())
    {
      UT_StringHolder status(rc == 0 ? "success" : "error");
      statusHandle.set(0, status);
    }

    if (rc != 0)
    {
      addError(SOP_MESSAGE, "MSBG solver failed - see console for details.");
      return UT_ERROR_ABORT;
    }

    return UT_ERROR_NONE;
  }
};

} // namespace

void newSopOperator(OP_OperatorTable *table)
{
  table->addOperator(new OP_Operator("msbg_solver",
                                      "MSBG Solver",
                                      SOP_MSBGSolver::constructor,
                                      SOP_MSBGSolver::buildTemplateList(),
                                      1,
                                      1));
}

// Required Houdini DSOVersion export.
extern "C" {
  DLLEXPORT void HoudiniDSOVersion(int *major, int *minor, int *build)
  {
    *major = UT_DSOVersion::getMajorVersion();
    *minor = UT_DSOVersion::getMinorVersion();
    *build = UT_DSOVersion::getBuildVersion();
  }
}

