/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#include "msbg_solver_module.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <tuple>

#include <tbb/tbb.h>

#include "thread.h"
#include "util.h"

namespace MSBG
{
namespace Solver
{

namespace
{

void limitBasePointCount(std::vector<Vec3Float> &points, size_t maxPoints)
{
  if (maxPoints == 0 || points.size() <= maxPoints)
  {
    return;
  }

  const size_t originalSize = points.size();

  std::vector<Vec3Float> limited;
  limited.reserve(maxPoints);

  double step = static_cast<double>(originalSize) / static_cast<double>(maxPoints);
  double cursor = 0.0;
  for (size_t i = 0; i < maxPoints; ++i)
  {
    size_t idx = static_cast<size_t>(cursor);
    if (idx >= originalSize)
    {
      idx = originalSize - 1;
    }
    limited.push_back(points[idx]);
    cursor += step;
  }

  limited.back() = points.back();

  TRCP(("Limiting base point count from %zu to %zu\n", originalSize, limited.size()));

  points.swap(limited);
}

using ParticleIdx = uint32_t;

struct Particle
{
  Vec3Float pos;
  ParticleIdx idxNext;
};

inline ParticleIdx getParticleIdx(int nInst, int iInst, int i)
{
  return nInst * static_cast<ParticleIdx>(iInst) + static_cast<ParticleIdx>(i) + 1;
}

inline Particle *getParticle(ParticleIdx ixp, Particle *base, uint64_t nMaxParticles)
{
  UT_ASSERT2(ixp > 0 && ixp < nMaxParticles + 1);
  return &base[ixp];
}

inline void pushToListMonotonic(Particle *particle, ParticleIdx idxParticle,
                                ParticleIdx *listRoot)
{
  ParticleIdx idx =
      InterlockedExchange(reinterpret_cast<LONG volatile *>(listRoot), idxParticle);
  particle->idxNext = idx;
}

std::tuple<std::unique_ptr<std::vector<Vec3Float>>, Vec3Float, Vec3Float>
readVerticesFromPLY(const char *filename)
{
  std::ifstream file(filename);
  if (!file)
  {
    TRCERR(("Unable to open PLY file '%s'\n", filename));
    return {std::make_unique<std::vector<Vec3Float>>(), Vec3Float(0.f), Vec3Float(0.f)};
  }

  std::string line;
  auto vertices = std::make_unique<std::vector<Vec3Float>>();
  Vec3Float bbMin(1e20f), bbMax(-1e20f);

  bool parsingHeader = true;
  int vertexCount = 0;
  int verticesParsed = 0;

  while (std::getline(file, line))
  {
    std::istringstream iss(line);
    if (parsingHeader)
    {
      if (line.substr(0, 10) == "element ve")
      {
        iss >> line >> line >> vertexCount;
        TRCP(("vertex count = %d\n", vertexCount));
      }
      else if (line == "end_header")
      {
        parsingHeader = false;
      }
    }
    else if (verticesParsed < vertexCount)
    {
      Vec3Float vertex;
      iss >> vertex.x >> vertex.y >> vertex.z;
      vertices->push_back(vertex);
      for (int k = 0; k < 3; ++k)
      {
        bbMin[k] = std::min(bbMin[k], vertex[k]);
        bbMax[k] = std::max(bbMax[k], vertex[k]);
      }
      verticesParsed++;
    }
  }

  return {std::move(vertices), bbMin, bbMax};
}

void computeBoundingBox(const std::vector<Vec3Float> &points, Vec3Float &bbMin,
                        Vec3Float &bbMax)
{
  bbMin = Vec3Float(1e20f);
  bbMax = Vec3Float(-1e20f);
  for (const Vec3Float &p : points)
  {
    for (int k = 0; k < 3; ++k)
    {
      bbMin[k] = std::min(bbMin[k], p[k]);
      bbMax[k] = std::max(bbMax[k], p[k]);
    }
  }
}

} // namespace

int RunPointCloudSolver(const SolverSettings &settings, const SolverCallbacks &callbacks)
{
  float rParticle = settings.particleRadius;
  float nbDist = settings.neighborDistance;

  UtTimer tm;

  std::unique_ptr<std::vector<Vec3Float>> basePoints;
  Vec3Float basePointsBMin, basePointsBMax;

  if (!settings.pointCloud.empty())
  {
    basePoints = std::make_unique<std::vector<Vec3Float>>(settings.pointCloud);
    limitBasePointCount(*basePoints, settings.maxBasePoints);
    computeBoundingBox(*basePoints, basePointsBMin, basePointsBMax);
  }
  else if (!settings.pointCloudFile.empty())
  {
    std::tie(basePoints, basePointsBMin, basePointsBMax) =
        readVerticesFromPLY(settings.pointCloudFile.c_str());
    limitBasePointCount(*basePoints, settings.maxBasePoints);
    if (basePoints && !basePoints->empty())
    {
      computeBoundingBox(*basePoints, basePointsBMin, basePointsBMax);
    }
  }
  else
  {
    TRCERR(("RunPointCloudSolver: no input point cloud specified."));
    return -1;
  }

  if (!basePoints || basePoints->empty())
  {
    TRCERR(("RunPointCloudSolver: the supplied point cloud is empty."));
    return -2;
  }

  Vec3Float basePointsSpan;
  float basePointsSpanMax = -1e20f;
  for (int k = 0; k < 3; ++k)
  {
    basePointsSpan[k] = basePointsBMax[k] - basePointsBMin[k];
    UT_ASSERT0(std::isfinite(basePointsSpan[k]) && basePointsSpan[k] >= 0.f);
    basePointsSpanMax = std::max(basePointsSpanMax, basePointsSpan[k]);
  }

  TRCP(("read %zu points. bbox = %g,%g,%g, spanMax=%g\n", basePoints->size(),
         basePointsSpan.x, basePointsSpan.y, basePointsSpan.z, basePointsSpanMax));

  tbb::parallel_sort(
      basePoints->begin(), basePoints->end(),
      [&](const Vector3Dim<float> &p1_, const Vector3Dim<float> &p2_) -> bool
      {
        Vec3Float p1 = (p1_ - basePointsBMin) / basePointsSpan,
                  p2 = (p2_ - basePointsBMin) / basePointsSpan;
        Vec4i ipos1 = Vec4i(p1.x, p1.y, p1.z, 0),
              ipos2 = Vec4i(p2.x, p2.y, p2.z, 0);
        uint64_t m1 = FMA::morton64_encode(ipos1),
                 m2 = FMA::morton64_encode(ipos2);
        return m1 < m2;
      });

  MSBG::MultiresSparseGrid *msbg = MSBG::MultiresSparseGrid::create(
      "MSBG_SOLVER", settings.resolutionX, settings.resolutionY, settings.resolutionZ,
      settings.blockSize, -1, 0, -1,
      MSBG::OPT_SINGLE_LEVEL | MSBG::OPT_SINGLE_CHANNEL_FLOAT);

  if (!msbg)
  {
    TRCERR(("RunPointCloudSolver: failed to create MSBG grid"));
    return -3;
  }

  if (settings.enableVisualizationOutput)
  {
    msbg->setVisOutDir(2, const_cast<char *>(settings.outputDirectory.c_str()));
  }

  LONG *blockActive = NULL;
  ParticleIdx *particlesPerBlock = NULL;
  Particle *particles = NULL;

  ALLOCARR0_(blockActive, LONG, msbg->nBlocks());
  ALLOCARR0_(particlesPerBlock, ParticleIdx, msbg->nBlocks());

  int chan = CH_UINT16_1;
  SparseGrid<uint16_t> *sg0 = msbg->getUint16Channel(chan, 0);

  uint64_t nBasePoints = basePoints->size();
  uint64_t nInst = settings.maxInstances > 0
                        ? std::min<uint64_t>(settings.maxInstances, nBasePoints)
                        : nBasePoints;
  uint64_t nMaxParticles = nInst * static_cast<uint64_t>(nBasePoints);

  UT_ASSERT0(nMaxParticles + 1 < 4 * (static_cast<uint64_t>(ONE_GB)));

  TRCP(("nBasePoints=%.0f  nInst=%.0f nMaxParticles=%.0f\n", (double)nBasePoints,
         (double)nInst, (double)nMaxParticles));

  float rScan = rParticle + nbDist;
  TRCP(("rParticle=%g, rScan=%g\n", rParticle, rScan));

  ALLOCARR_(particles, Particle, nMaxParticles + 1);

  uint64_t nActParticles = 0;

  TIMER_START(&tm);
  {
    int bxMax = sg0->nbx() - 1, byMax = sg0->nby() - 1, bzMax = sg0->nbz() - 1;
    const float scale2DestBlockGrid = 1.0f / static_cast<float>(sg0->bsx());
    const float rScanBsx = rScan / static_cast<float>(sg0->bsx());
    Vec4i bposMax(bxMax, byMax, bzMax, INT32_MAX);

    using ThreadLocals = struct
    {
      size_t nActParticles;
    };

    ThrRunParallel<ThreadLocals>(
        nInst, nullptr,
        [&](ThreadLocals &tls, int /*tid*/, size_t iInst)
        {
          Vec4f baseMin(basePointsBMin.x, basePointsBMin.y, basePointsBMin.z, 0);
          float baseScale = 1.0f / basePointsSpanMax, scale;
          Vec4f pos0;
          scale = sg0->sxyzMin() * (settings.testCase == 2 ? 0.005f : 0.01f);
          Vec3Float &pos_ = (*basePoints)[iInst];
          Vec4f pos = sg0->sxyzMax() * baseScale *
                      (Vec4f(pos_.x, pos_.y, pos_.z, 0.f) - baseMin);
          pos0 = sg0->sxyzMax() * .2f + 0.6f * pos;

          for (int i = 0; i < static_cast<int>(nBasePoints); i++)
          {
            Vec3Float &p = (*basePoints)[i];
            Vec4f posSample = baseScale * (Vec4f(p.x, p.y, p.z, 0.f) - baseMin);
            posSample = pos0 + scale * posSample;
            if (!msbg->isInDomainRange(posSample))
              continue;
            Vec4i ipos = truncate_to_int(posSample);
            if (!sg0->inRange(ipos))
              continue;

            ParticleIdx ixp = getParticleIdx(nInst, static_cast<int>(iInst), i);
            Particle *part = getParticle(ixp, particles, nMaxParticles);
            part->pos = Vec3Float(posSample);

            Vec4f bpos = posSample * scale2DestBlockGrid;
            Vec4i bpos1 = max(truncate_to_int(bpos - rScanBsx), 0),
                  bpos2 = min(truncate_to_int(bpos + rScanBsx), bposMax);
            int bx1 = vget_x(bpos1), by1 = vget_y(bpos1), bz1 = vget_z(bpos1);
            int bx2 = vget_x(bpos2), by2 = vget_y(bpos2), bz2 = vget_z(bpos2);
            for (int bz = bz1; bz <= bz2; bz++)
              for (int by = by1; by <= by2; by++)
                for (int bx = bx1; bx <= bx2; bx++)
                {
                  int bid = sg0->getBlockIndex(bx, by, bz);
                  UT_ASSERT2(bid >= 0 && bid < sg0->nBlocks());
                  InterlockedIncrement(
                      reinterpret_cast<volatile LONG *>(&blockActive[bid]));
                }

            int bid = sg0->getBlockIndex(sg0->getBlockCoords(ipos));
            pushToListMonotonic(part, ixp, &particlesPerBlock[bid]);
            tls.nActParticles++;
          }
        },
        [&](ThreadLocals &tls, int /*tid*/)
        { nActParticles += tls.nActParticles; });
  }

  TIMER_STOP(&tm);
  TRCP(("CPU (determine active blocks) %.2f sec, %.0f points/sec)\n",
         (double)TIMER_DIFF_MS(&tm) / 1000.,
         ((double)nActParticles) / (double)(TIMER_DIFF_MS(&tm) / 1000.0)));

  TRCP(("nActParticles = %.0f\n", (double)nActParticles));

  LongInt nActiveBlocks = 0;
  {
    TIMER_START(&tm);

    std::unique_ptr<int[]> blockLevels(new int[msbg->nBlocks()]);
    for (int i = 0; i < msbg->nBlocks(); i++)
    {
      int level = msbg->getNumLevels() - 1;
      if (blockActive[i])
      {
        level = 0;
        nActiveBlocks++;
      }
      blockLevels[i] = level;
    }
    msbg->regularizeRefinementMap(blockLevels.get());
    msbg->setRefinementMap(blockLevels.get(), NULL, -1, NULL, false);

    TIMER_STOP(&tm);
    TRCP(("CPU (Set refinement map) %.2f sec, %.0f points/sec)\n",
           (double)TIMER_DIFF_MS(&tm) / 1000.,
           ((double)nActParticles) / (double)(TIMER_DIFF_MS(&tm) / 1000.0)));
  }

  TRCP(("%s: nActiveBlocks=%ld/%d %.2f%%\n", UT_FUNCNAME, nActiveBlocks,
         msbg->nBlocks(), 100. * nActiveBlocks / (double)msbg->nBlocks()));

  sg0->reset();
  sg0->prepareDataAccess(chan);
  sg0->setEmptyValue(renderDensFromFloat_<RenderDensity>(0.0f));
  sg0->setFullValue(renderDensFromFloat_<RenderDensity>(1.0f));

  std::vector<int> activeBlocks;
  {
    using ThreadLocals = struct
    {
      std::vector<int> activeBlocks;
    };

    ThrRunParallel<ThreadLocals>(
        sg0->nBlocks(),
        [&](ThreadLocals &tls, int /*tid*/)
        { tls.activeBlocks.reserve(256); },
        [&](ThreadLocals &tls, int /*tid*/, int bid)
        {
          BlockInfo *bi = msbg->getBlockInfo(bid);
          FLAG_RESET(bi->flags, BLK_EXISTS);
          if (bi->level > 0)
          {
            sg0->setEmptyBlock(bid);
            return;
          }

          tls.activeBlocks.push_back(bid);

          FLAG_SET(bi->flags, BLK_EXISTS);
          RenderDensity *data = sg0->getBlockDataPtr(bid, 1, 0);
          for (int i = 0; i < sg0->nVoxelsInBlock(); i++)
            data[i] = renderDensFromFloat_<RenderDensity>(1.0f);
        },
        [&](ThreadLocals &tls, int /*tid*/)
        { UT_VECTOR_APPEND(activeBlocks, tls.activeBlocks); },
        0, true);
  }

  TIMER_STOP(&tm);
  TRCP(("CPU (allocate & initialize %ld active blocks) %.2f sec, %.0f points/sec)\n",
         (long)activeBlocks.size(), (double)TIMER_DIFF_MS(&tm) / 1000.,
         ((double)nActParticles) / (double)(TIMER_DIFF_MS(&tm) / 1000.0)));

  UT_ASSERT0(static_cast<size_t>(nActiveBlocks) == activeBlocks.size());

  TRCP(("Splatting particles.\n"));

  TIMER_START(&tm);

  const float distSqMax = MT::sqf(rParticle + nbDist);
  const float distSqMaxInv = 1.0f / distSqMax;

  std::vector<int> activeBlocksPerCol[8];
  for (int bid : activeBlocks)
  {
    Vec4i bpos = sg0->getBlockCoordsById(bid);
    int bx = vget_x(bpos), by = vget_y(bpos), bz = vget_z(bpos);
    int icol = getBlockColor8(bx, by, bz);
    UT_ASSERT0(icol >= 0 && icol < 8);
    activeBlocksPerCol[icol].push_back(bid);
  }

  for (int icol = 0; icol < 8; icol++)
  {
    std::vector<int> *blockList = &activeBlocksPerCol[icol];
    if (blockList->empty())
      continue;

    ThrRunParallel<int>(
        blockList->size(), nullptr,
        [&](int &, int /*tid*/, int ibid)
        {
          int bid = (*blockList)[ibid];
          Vec4i bpos = sg0->getBlockCoordsById(bid);
          int bx = vget_x(bpos), by = vget_y(bpos), bz = vget_z(bpos);
          UT_ASSERT0(getBlockColor8(bx, by, bz) == icol);

          Particle *p = NULL;
          for (ParticleIdx ixp = particlesPerBlock[bid]; ixp; ixp = p->idxNext)
          {
            p = getParticle(ixp, particles, nMaxParticles);

            Vec4f pos0 = Vec4f(p->pos.x, p->pos.y, p->pos.z, 0.0f);
            Vec4i ipos1 = max(truncate_to_int(ceil(pos0 - rScan - 0.5f)), 0),
                  ipos2 = min(truncate_to_int(floor(pos0 + rScan - 0.5f)), sg0->v4iDomMax());

            int ix1 = vget_x(ipos1), ix2 = vget_x(ipos2),
                iy1 = vget_y(ipos1), iy2 = vget_y(ipos2),
                iz1 = vget_z(ipos1), iz2 = vget_z(ipos2);
            Vec4f pshift = pos0 - 0.5f;
            float x0 = vfget_x(pshift), y0 = vfget_y(pshift), z0 = vfget_z(pshift);

            const int bsxLog2 = sg0->bsxLog2(),
                      bsx2Log2 = sg0->bsx2Log2(),
                      bsxMask = sg0->bsx() - 1,
                      nx = sg0->nbx(),
                      nxy = sg0->nbxy();

            for (int iz = iz1; iz <= iz2; iz++)
            {
              float dz = iz - z0;

              int ibz = (iz >> bsxLog2) * nxy,
                  ivz = (iz & bsxMask) << bsx2Log2;

              for (int iy = iy1; iy <= iy2; iy++)
              {
                float dy = iy - y0,
                      distSqZY = dz * dz + dy * dy;

                int ibzy = ibz + (iy >> bsxLog2) * nx,
                    ivzy = ivz + ((iy & bsxMask) << bsxLog2);

                for (int ix = ix1; ix <= ix2; ix++)
                {
                  float dx = ix - x0,
                        distSq = distSqZY + dx * dx;

                  if (distSq > distSqMax)
                    continue;

                  int ib = ibzy + (ix >> bsxLog2),
                      iv = ivzy + (ix & bsxMask);

                  UT_ASSERT2(ib >= 0 && ib < sg0->nBlocks());
                  UT_ASSERT2(iv >= 0 && iv < sg0->nVoxelsInBlock());

                  SBG::Block<RenderDensity> *block = sg0->getBlock(ib);

                  if (!sg0->isValueBlock(block))
                  {
                    goto loop_exit;
                  }

                  uint16_t val =
                      renderDensFromFloat_<RenderDensity>(distSq * distSqMaxInv);
                  block->_data[iv] = std::min(block->_data[iv], val);
                }
              }
            }
          loop_exit:;
          }
        });
  }

  TIMER_STOP(&tm);
  TRCP(("CPU (splatting particles) %.2f sec, %.0f points/sec)\n",
         (double)TIMER_DIFF_MS(&tm) / 1000.,
         ((double)nActParticles) / (double)(TIMER_DIFF_MS(&tm) / 1000.0)));

  FREEMEM(particlesPerBlock);
  FREEMEM(particles);

  TRCP(("Finalizing splatting pass.\n"));

  TIMER_START(&tm);

  ThrRunParallel<int>(
      activeBlocks.size(), nullptr,
      [&](int &, int /*tid*/, int ibid)
      {
        int bid = activeBlocks[ibid];
        UT_ASSERT0(sg0->isValueBlock(bid));
        RenderDensity *data = sg0->getBlockDataPtr(bid);

        for (int vid = 0; vid < sg0->nVoxelsInBlock(); vid++)
        {
          float distSq = renderDensToFloat_(data[vid]);
          uint16_t uiVal = 0;
          if (distSq < 0.999f)
          {
            distSq *= distSqMax;
            float dist = sqrtf(distSq) - rParticle;
            float f = 1.0f - MT_LINSTEPF(-rParticle, (rParticle + nbDist), dist);
            uiVal = renderDensFromFloat_<RenderDensity>(f);
          }
          data[vid] = uiVal;
        }
      });

  TIMER_STOP(&tm);

  TRCP(("CPU (finalizing splatting) %.2f sec\n", (double)TIMER_DIFF_MS(&tm) / 1000.));

  if (callbacks.afterParticleSplat)
  {
    callbacks.afterParticleSplat(msbg, chan);
  }

  int nIterPde = settings.pdeIterations > 0 ? settings.pdeIterations : 5;

  TRCP(("Applying mean curvature smoothing PDE (%d iterations)\n", nIterPde));

  TIMER_START(&tm);

  float timeStep = settings.timeStep > 0.f ? settings.timeStep
                                           : (settings.testCase == 2 ? 0.1f : 0.05f);

  msbg->applyChannelPdeFast<RenderDensity>(
      chan, CH_NULL, CH_NULL, &activeBlocks,
      -(PDE_MEAN_CURVATURE + OPT_8_COLOR_SCHEME), TRUE, 0, 0, 0, nIterPde, 1.0f,
      timeStep, 0, NULL, 0, 0, 0, 0, 0, 0, 0, 0);

  TIMER_STOP(&tm);

  TRCP(("CPU (Mean Curvature PDE) %.2f seconds.\n", (double)TIMER_DIFF_MS(&tm) / 1000.));

  if (callbacks.afterPdeSolve)
  {
    callbacks.afterPdeSolve(msbg, chan);
  }

  FREEMEM(blockActive);
  MSBG::MultiresSparseGrid::destroy(msbg);

  return 0;
}

} // namespace Solver
} // namespace MSBG

