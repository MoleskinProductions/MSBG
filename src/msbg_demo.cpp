/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <tuple> 
#include <iostream>
#include <fstream>
#include <tbb/tbb.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <errno.h>
#include <unistd.h>

#include "msbg.h"
#include "msbg_solver_module.h"

#include "render.h"

#define VIS_OUTPUT_DIR "out_msbg_demo"

float camPos[3] = {.5f,.8f,.7f},
      camLookAt[3] = {.5,.5,.5},
      camLight[3] = {.5,5,5},
      camZoom =1.0f;
int camRes[2] = {3840,2160};

static void render_scene( int testCase, MSBG::MultiresSparseGrid *msbg, int chan, 
     			  PnlPanel *pnl )
{
  SBG::SparseGrid<uint16_t> *sg = msbg->getUint16Channel(chan,0);

  UtTimer tm;
  int sx=camRes[0], sy=camRes[1];

  TRCP(("Rendering image %dx%d.\n",sx,sy));

  RDR::RaymarchRenderer renderer(sx,sy, sg);
#if 1

  #if 1
  renderer.setCamera( {camPos[0],camPos[1],camPos[2]}, 
      		      {camLookAt[0],camLookAt[1],camLookAt[2]}, camZoom); 


//  renderer.setCamera({.5,.8,.7}, {.5,.5,.5}, 1); 
  renderer.setSunLight({camLight[0],camLight[1],camLight[2]});
  renderer.setSurfaceColor({0.8f,0.6f,0.4f});
  #else

  renderer.setCamera({.5,.8,.8}, {.5,.5,.5}, 1.f);
//  renderer.setCamera({.5,1,1}, {.5,.5,.5}, 1.5f);
  

  renderer.setSunLight({.5,5,5});
  renderer.setSurfaceColor({0.8f,0.6f,0.4f});
  #endif

#else
  renderer.setCamera({.5,.5,.65}, {.5,.5,.5}, 0.8f);
//  renderer.setCamera({0,.5,1}, {0.4,0.5,0.5}, 2.5f);
  renderer.setSunLight({.5,5,5});
  renderer.setSurfaceColor({0.8f,0.6f,0.4f});
#endif


  TIMER_START(&tm);

  BmpBitmap *B = renderer.render();

  TIMER_STOP(&tm);
  TRCP(("CPU (rendering) %.2f sec, %.0f pixels/sec)\n",
    (double)TIMER_DIFF_MS(&tm)/1000.,
    ((double)sy*sy)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));	

  PnlShowBitmap2( pnl, B );  

  pnl->totalTime = 0;
  msbg->saveVisualizationBitmap( B, pnl->title );

  TRCP(("Output images saved to '%s/'.\n",VIS_OUTPUT_DIR));
  
  BmpDeleteBitmap(&B);
}

/////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////
int msbg_test_multires(int bsx0, int sx, int sy, int sz)
{
  using namespace MSBG;
  using namespace SBG;
  using MSBG::MultiresSparseGrid;

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 	Create MSBG grid
  //
  /////////////////////////////////////////////////////////////////////////////////

  MSBG::MultiresSparseGrid *msbg = 
    MSBG::MultiresSparseGrid::create( "TEST", sx,sy,sz, bsx0,
				      -1, 0,-1, 
				      MSBG::OPT_SINGLE_LEVEL |
					MSBG::OPT_SINGLE_CHANNEL_FLOAT );

  msbg->setVisOutDir( 2, (char*)VIS_OUTPUT_DIR );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 		Free resources 
  //
  /////////////////////////////////////////////////////////////////////////////////

  MSBG::MultiresSparseGrid::destroy(msbg);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////
int msbg_test_sparse(int testCase, const char *basePointsFile,
                     int bsx0, int sx, int sy, int sz)
{
  MSBG::Solver::SolverSettings settings;
  if (basePointsFile)
    settings.pointCloudFile = basePointsFile;
  settings.blockSize = bsx0;
  settings.resolutionX = sx;
  settings.resolutionY = sy;
  settings.resolutionZ = sz;
  settings.testCase = testCase;
  settings.pdeIterations = 5;
  settings.outputDirectory = VIS_OUTPUT_DIR;
  settings.enableVisualizationOutput = true;

  MSBG::Solver::SolverCallbacks callbacks;
  callbacks.afterPdeSolve = [testCase](MSBG::MultiresSparseGrid *grid, int chan) {
    static int pnl = -1;
    render_scene(testCase, grid, chan,
                 MiGetAuxPanel2(&pnl, "scene_after_PDE_smoothing"));
  };

  return MSBG::Solver::RunPointCloudSolver(settings, callbacks);
}

int test_MSBG_0(int opcode, char *fpath, int bsx0, int sx, int sy, int sz)
{
#if 1
  return 1;
#else
  using namespace SBG;
  //int rcThis=0;
  //char *tok;

  //bool doWrite = opcode & 1;

//  bool doLevelSet = true;
  bool doLevelSet = false;

  //#define DENSITY_8_BIT

  std::vector<SparseGrid<float>*> sbgGrids;
  std::vector<char *> sbgNames;
  std::vector<unsigned> sbgOptions;

  int bsx = bsx0;
  MSBG::MultiresSparseGrid *msbg = 
    MSBG::MultiresSparseGrid::create( "TEST", sx,sy,sz, bsx,
				      -1, 0,-1, 
				      MSBG::OPT_SINGLE_LEVEL |
					MSBG::OPT_SINGLE_CHANNEL_FLOAT );
  UT_ASSERT0(msbg);

  int *blockRefinementMap = NULL;
  ALLOCARR_(blockRefinementMap,int,msbg->nBlocks());

 {

  #ifdef DENSITY_8_BIT
  typedef uint8_t DataType;
  int chan = MSBG::CH_GEN_UINT8;
  #else
  typedef float DataType;
  int chan = MSBG::CH_FLOAT_1;
  #endif
  msbg->resetChannel( chan );
  msbg->prepareDataAccess( chan );
  
  #ifdef DENSITY_8_BIT
  SparseGrid<uint8_t> *sbg = msbg->getUint8Channel(chan,0);
  #else
  SparseGrid<float> *sbg = msbg->getFloatChannel(chan,0);
  #endif

  sbg->setEmptyValue(0.0f);
  sbg->setFullValue(1.0f);

  #if 1
  {
    //
    // Fill with sphere surface (Enright)
    //

#if 1
    int fbmBaseShape = 1;
    Vec4f center_(0.5,0.5,.5,0.0f),
	  radius_(0.3,.3,.3,0);
    float resolution = msbg->sg0()->sxyzMax();
    float fbmSeed = 4711,
          fbmDepth_ = .1,
	  fbmFscale_ = 0.05,  
	  fbmMinscale = 0.5, // in fractions of voxel
	  fbmAmpl = .2, //0, //1, //1,
	  fbmAlpha = 1.5,
	  fbmOctLeadin = 1,
	  fbmThr = 0;
#endif
    USE(fbmSeed); USE(fbmFscale); USE(fbmMinscale); USE(fbmAlpha);
    USE(fbmOctLeadin); 
    

#if 0
    Vec4f sphereCenter(0.5,0.5,.5,0.0f);
    float sphereRadius = 0.4,
	  resolution = msbg->sg0()->sxyzMax();
    float fbmSeed = 4711,
	  fbmDepth = 4,
	  fbmFscale = 0.4*resolution,
	  fbmMinscale = 0.5,
	  fbmAmpl = 0, //1, //1,
	  fbmAlpha = 2,
	  fbmOctLeadin = 1,
	  fbmThr = 0;
#endif

#if 0
    Vec4f sphereCenter(0.5,0.5,.5,0.0f);
    float sphereRadius = 0.4,
	  resolution = msbg->sg0()->sxyzMax();
    float fbmSeed = 4711,
	  fbmDepth = 0.1*resolution,
	  fbmFscale = 0.1*resolution,
	  fbmMinscale = 0.5,
	  fbmAmpl = 1, //1, //1,
	  fbmAlpha = 1.5,
	  fbmOctLeadin = 1,
	  fbmThr = 1.2;
#endif

#if 0
    Vec4f sphereCenter(0.5,0.5,.5,0.0f);
    float sphereRadius = 0.2,
	  resolution = msbg->sg0()->sxyzMax();
    float fbmSeed = 4711,
	  fbmDepth = 0.4*resolution,
	  fbmFscale = 0.1*resolution,
	  fbmMinscale = 0.5,
	  fbmAmpl = 0, //1, //1,
	  fbmAlpha = 1.5,
	  fbmOctLeadin = 1,
	  fbmThr = .8;
#endif

#if 0
    Vec4f sphereCenter(0.5,0.5,.5,0.0f);
    float sphereRadius = 0.4,
	  resolution = msbg->sg0()->sxyzMax();
    float fbmSeed = 4711,
	  fbmDepth = 10, //0.2*resolution,
	  fbmFscale = 0.1*resolution,
	  fbmMinscale = 0.5,
	  fbmAmpl = 0, //1, //1,
	  fbmAlpha = 1.5,
	  fbmOctLeadin = 1,
	  fbmThr = .8;
#endif
    sphereCenter *= resolution;
    sphereRadius *= resolution;

    float sphereRadius0 = sphereRadius - fbmDepth,
	  sphereRadius1 = sphereRadius + fbmDepth;

    int sx = sbg->sx(),
	sy = sbg->sx(),
	sz = sbg->sx();
    USE(sx);USE(sy);USE(sz);

    UtTimer tm;

    TIMER_START(&tm);

    LongInt nActVoxels = 0;

    using ThreadLocals = struct
    {
      DataType *dataTmp;
      size_t nActVoxels;
    };

    ThrRunParallel<ThreadLocals>( sbg->nBlocks(), 
      [&]( ThreadLocals &tls, int tid )  // Initialize locals
      {
	ALLOCARR_(tls.dataTmp,DataType,sbg->nVoxelsInBlock());
	tls.nActVoxels = 0;
      },

      [&](ThreadLocals &tls, int tid, int bid) 
      {
        float phiMin=1e20,
	      phiMax=-1e20;
	// 
	// Perlin 'hypertexture'
	//	
	Vec4f center = resolution * center_,
	      radius = resolution * radius_;
	float depth = resolution * depth_;      

	auto baseShape = [&]( const Vec4f& pos ) -> float
	{
	  switch( fbmBaseShape )
	  {
	    case 1:  // Sphere
	    {
	      const float rSq = v3normSq( pos ),
			  r0Sq = sqf(radius[0]),
			  r1Sq = sqf(radius[0]-depth);

	      float f = rSq <= r1Sq ? 1.0f :
			rSq >= r0Sq ? 0.0f :
			    ( r0Sq - rSq ) * 1.0f/(r0Sq - r1Sq);
	      return f;
	    }

	    case 4:  // Ellipsoid (https://iquilezles.org/articles/distfunctions/)
	    {
	      const float k0 = sqrtf(v3normSq( pos * 1.0f/radius )),
		          k1 = sqrtf(v3normSq( pos * 1.0f/(radius*radius) ));
	      const float dist = k0*(k0-1.0f)/k1;
	      float f = 1.0f - MT_LINSTEPF( -depth, 0, dist ); 
	      return f;
	    }

	    case 7:  // Plane
	    {
	      Vec4f orig = center,
		    normal = radius * 1.0f/sqrtf(v3normSq(radius)),
		    tmp = normal * (pos-orig);

	      float dist = vfget_x(tmp)+vfget_y(tmp)+vfget_z(tmp);

	      float f = 1.0f-MT_LINSTEPF( 0, depth, dist ); 
	      return f;
	    }

	    default:
	    UT_ASSERT0(FALSE);
	    return 0;
	    break;
	  }
	};

	unsigned bflags = 0;

	Vec4i bpos = sbg->getBlockCoordsById(bid);
	Vec4f boxCenter = ((vint2float(bpos) + 0.5f) * (float)bsx);
	vzero_w(boxCenter);
	float dist = baseShape( boxCenter );

	blockRefinementMap[bid] = msbg->getNumLevels()-1;

	if( dist > (float)MT_SQRT3*bsx*0.5f )
	{
	  sbg->setEmptyBlock(bid);
	  return;
	}
	if( dist < (float)MT_SQRT3*bsx*0.5f )
	{
	  sbg->setFullBlock(bid);
	  return;
	}




	#if 1
	Vec4i bpos = sbg->getBlockCoordsById(bid);
	Vec4f boxCenter = ((vint2float(bpos) + 0.5f) * (float)bsx);
	vzero_w(boxCenter);
	float dist = sqrtf(v3normSq( boxCenter - sphereCenter ));

	blockRefinementMap[bid] = msbg->getNumLevels()-1;

	if( doLevelSet )
	{
	  if( dist > sphereRadius1 + (float)MT_SQRT3*bsx*0.5f )
	  {
	    sbg->setEmptyBlock(bid);
	    return;
	  }
	  if( dist < sphereRadius0 - (float)MT_SQRT3*bsx*0.5f )
	  {
	    sbg->setFullBlock(bid);
	    return;
	  }
	}
	else
	{
	  if( dist > sphereRadius1 + (float)MT_SQRT3*bsx*0.5f )
	  {
	    sbg->setEmptyBlock(bid);
	    return;
	  }
	  if( dist < sphereRadius0 - (float)MT_SQRT3*bsx*0.5f )
	  {
	    sbg->setFullBlock(bid);
	    return;
	  }
	}
	#endif

	{
	  float phiMin=1e20,
	        phiMax=-1e20;

	  SBG_FOR_EACH_VOXEL_IN_BLOCK(sbg,bid,x,y,z,vid)
	  {
	    Vec4f pos = Vec4f( x+0.5f, y+0.5f, z+0.5f, 0.0f);
	    float dist = sqrtf( v3normSq( pos - sphereCenter )) - sphereRadius;

	    float f = doLevelSet ? 	    
	      1.0f - MT_LINSTEPF( -fbmDepth, fbmDepth, dist ) : 
	      1.0f - MT_LINSTEPF( -fbmDepth, fbmDepth, dist );
	    
	    UT_ASSERT2(!( f<0.0f || f>1.0f || !std::isfinite(f)));
#if 1
	    float u = PNS::genFractalScalar3D( sx, sy, sz,
				pos[0]+sx,pos[1]+sy,pos[2]+sz,
				fbmSeed, fbmFscale, fbmMinscale, fbmOctLeadin, fbmAlpha, 0
				//713647, fscale, 2, 0, 2, 0 
				);
	    //TRCP(("f=%g u=%g\n",f,u));
#else
	    float u = 0;
#endif
	    f += fbmAmpl * u;
	    float phi = std::min(std::max(0.0f,f - fbmThr),1.0f);

	    //phi = 0.5; 
	    //phi = powf(x/(double)sbg->sx(),2.);

	    #ifdef DENSITY_8_BIT
	    tls.dataTmp[vid] = renderDensFromFloat_<DataType,false,false>(phi);
	    /*if(phi>0) 
	    {
	      TRCP(("%g -> %d\n",phi,tls.dataTmp[vid]));
	    }*/
	    #else
	    tls.dataTmp[vid] = phi;
	    #endif

	    phiMin = std::min(phiMin,phi);
	    phiMax = std::max(phiMax,phi);
	  }

	  if( phiMax < MT_NUM_EPS )
	  {
	    sbg->setEmptyBlock(bid);
	  }
	  else
	  {
	    blockRefinementMap[bid] = 0;
	    DataType *data = sbg->getBlockDataPtr(bid,1,0);
	    tls.nActVoxels += sbg->nVoxelsInBlock();
	    for(int i=0;i<sbg->nVoxelsInBlock();i++) data[i] = tls.dataTmp[i];
	  }
	}
      },

      [&]( ThreadLocals &tls, int tid )  // Reduce locals
      {
	FREEMEM(tls.dataTmp);
	nActVoxels += tls.nActVoxels;	
      }
    );

    msbg->regularizeRefinementMap( blockRefinementMap );  
    msbg->setRefinementMap( blockRefinementMap, NULL, -1, NULL,			
			       false, false );
    msbg->_blocksValue[0].clear();

    TIMER_STOP(&tm);
    TRCP(("Create: CPU=%.2f sec, %.0f voxels/sec. (%.0f/%.1f%%) voxels\n",
	  (double)TIMER_DIFF_MS(&tm)/1000.,
	(double)nActVoxels/(double)(TIMER_DIFF_MS(&tm)/1000.0),
	(double)nActVoxels,
	100.*nActVoxels/(double)sbg->nTotVirtVoxels()      
	));
  
    TRCP(("CPU %.2f sec  \n",UT_FUNCNAME,(double)TIMER_DIFF_MS(&tm)/1000.));

    sbg->showMemUsage(1);

    #if 1
    {
      static int pnl=-1;
      msbg->visualizeSlices( MiGetAuxPanel2(&pnl,"TEST"), chan );	
    }
    #endif
    #if 0
    {
      static int pnl=-1;
      msbg->visualizeSlices( MiGetAuxPanel2(&pnl,"TEST"),
				 0,IP_NEAREST,NULL, 0,0,0,0,
				 sbg );	
    }
    #endif
  }
  #endif	


  sbgGrids.push_back( sbg );
  sbgNames.push_back( doLevelSet ? (char*)"LevelSet" : (char*)"density" );
  sbgOptions.push_back( doLevelSet ? 1 : 0 );
  
  int rc = SBG::writeAsVDB(fpath,256,sbgGrids,sbgNames,sbgOptions);
  TRCP(("SBG::writeAsVDB('%s') -> %d\n",fpath,rc));

  }

rcCatch:
  MSBG::MultiresSparseGrid::destroy(msbg);
  FREEMEM(blockRefinementMap);
  return 0;
#endif
}

