// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include <gtest/gtest.h>

#include <pbrt/pbrt.h>

#include <pbrt/bsdf.h>
#include <pbrt/interaction.h>
#include <pbrt/options.h>
#include <pbrt/paramdict.h>
#include <pbrt/shapes.h>
#include <pbrt/util/image.h>
#include <pbrt/util/log.h>
#include <pbrt/util/memory.h>
#include <pbrt/util/parallel.h>
#include <pbrt/util/print.h>
#include <pbrt/util/rng.h>
#include <pbrt/util/sampling.h>
#include <pbrt/util/spectrum.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>

using namespace pbrt;

TEST(SGGX, Logging) {
    Vector3f initial_vec = Vector3f(0, 0, 1);

    printf("Logging SGGX");
    Transform m_transform = RotateY(120);
    Vector3f t = m_transform(initial_vec);
    //Vector3f t = Vector3f(0, 0, 1);
    t = Normalize(t);
    //sggx::Ellipsoid S = sggx::Ellipsoid::fromFiber(t, rough);
    Float rough = 0.01;
    
    sggx::Ellipsoid S = sggx::Ellipsoid::fromSurface(t, rough);
    Vector3f p(0,0,1);
    FILE *f = fopen("../../scripts/SGGX/log_sggx_pbrt.txt", "w");
    SGGXPhaseFunctionNormalDistribution phase(rough,-90,0.0f);
    //SGGXPhaseFunction m_phase(rough,Float(90.f));

    for (int j = 0; j <= 360; j++)
    {
        for (int i = 0; i < 180; i++)
        {
            Float radi = 2 * M_PI*(1.0*i / 360);
            Float radj = 2 * M_PI*(1.0*j / 360);

            Vector3f wi = Vector3f(sinf(radi)*sinf(radj), sinf(radi)*cosf(radj), cosf(radi));
            //Float D_sggx = D(wi, S);

            // Float D_sggx = m_phase.PDF(p,wi,Point2f(0.0f,0.0f) );
            //Float D_sggx = m_phase.PDF(p,wi);
            //Check Projected Area
            //Float D_sggx = projArea(wi, S);
            Float D_sggx= phase.projectedArea(wi,Point2f(0.5f,0.5f));

            fprintf(f, "%f ", D_sggx);
        }

        int i = 180;
        Float radi = 2 * M_PI*(1.0*i / 360);
        Float radj = 2 * M_PI*(1.0*j / 360);

        Vector3f wi = Vector3f(sinf(radi)*sinf(radj), sinf(radi)*cosf(radj), cosf(radi));
        //Float D_sggx = m_phase.PDF(p,wi );
        // Float D_sggx = m_phase.PDF(p,wi,Point2f(0.0f,0.0f) );
        
        //Float D_sggx = D(wi, S);
        
        // Check Projected Area
        //Float D_sggx = projArea(wi, S);
        Float D_sggx= phase.projectedArea(wi,Point2f(0.5f,0.5f));
            
        fprintf(f, "%f\n", D_sggx);

        //fprintf_s(f, ";\n");
        //fprintf(f, "\n");
    }
    
    
    fclose(f);

    EXPECT_TRUE(true);
}