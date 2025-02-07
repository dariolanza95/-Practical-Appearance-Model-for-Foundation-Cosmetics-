// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_MEDIA_H
#define PBRT_MEDIA_H

#include <pbrt/pbrt.h>

#include <pbrt/base/medium.h>
#include <pbrt/interaction.h>
#include <pbrt/paramdict.h>
#include <pbrt/textures.h>
#include <pbrt/util/colorspace.h>
#include <pbrt/util/error.h>
#include <pbrt/util/memory.h>
#include <pbrt/util/parallel.h>
#include <pbrt/util/print.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/scattering.h>
#include <pbrt/util/spectrum.h>
#include <pbrt/util/transform.h>

#include <nanovdb/NanoVDB.h>
#include <nanovdb/util/GridHandle.h>
#include <nanovdb/util/SampleFromVoxels.h>
#ifdef PBRT_BUILD_GPU_RENDERER
#include <nanovdb/util/CudaDeviceBuffer.h>
#endif  // PBRT_BUILD_GPU_RENDERER

#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

#include "sggx.h"
#include <iostream>

namespace pbrt {

// Media Function Declarations
bool GetMediumScatteringProperties(const std::string &name, Spectrum *sigma_a,
                                   Spectrum *sigma_s, Allocator alloc);

// HGPhaseFunction Definition
class HGPhaseFunction {
  public:
    // HGPhaseFunction Public Methods
    HGPhaseFunction() = default;
    PBRT_CPU_GPU
    HGPhaseFunction(Float g) : g(g) {}

    PBRT_CPU_GPU
    Float p(Vector3f wo, Vector3f wi) const { return HenyeyGreenstein(Dot(wo, wi), g); }

    PBRT_CPU_GPU
    pstd::optional<PhaseFunctionSample> Sample_p(Vector3f wo, Point2f u) const {
        Float pdf;
        Vector3f wi = SampleHenyeyGreenstein(wo, g, u, &pdf);
        return PhaseFunctionSample{pdf, wi, pdf};
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi) const { return p(wo, wi); }

    static const char *Name() { return "Henyey-Greenstein"; }

    std::string ToString() const;

  private:
    // HGPhaseFunction Private Members
    Float g;
};

// SGGXPhaseFunction Definition

class SGGXPhaseFunction {
  public:
    // SGGXFunction Public Methods
    SGGXPhaseFunction() = default;
        PBRT_CPU_GPU
    SGGXPhaseFunction(Float sigma,Vector3f dir){
        std::cout<<"The following initialization is not supported anymore"<<std::endl;
        this->sigma = sigma;
    }
    PBRT_CPU_GPU
    SGGXPhaseFunction(Float sigma,Float rot_angle){
        
        this->sigma = sigma;
        printf("Sigma = %f\n", sigma);
        //sigma = 0.1;
        this->SetMicroflakes(sigma, rot_angle);
        //this->m_mean = 80.0f;
        //this->m_std = 1.0f;
        //this->SetMicroflakes(t, sigma);//Fiber-like distributions
        printf("Initialization works");
    }

    PBRT_CPU_GPU
    SGGXPhaseFunction(sggx::Ellipsoid S) : S(S) {}


    Transform RotateZ(Float theta) {
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    SquareMatrix<4> m(cosTheta, -sinTheta, 0, 0,
                      sinTheta,  cosTheta, 0, 0,
                             0,         0, 1, 0,
                             0,         0, 0, 1);
    return Transform(m, Transpose(m));
}
Transform RotateY(Float theta) const{
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    SquareMatrix<4> m( cosTheta, 0, sinTheta, 0,
                              0, 1,        0, 0,
                      -sinTheta, 0, cosTheta, 0,
                              0, 0,        0, 1);
    return Transform(m, Transpose(m));
}
inline Float GaussianSample(Float x, Float mu = 0, Float sigma = 1) const {
        return 1 / std::sqrt(2 * Pi * sigma * sigma) *
               std::exp(-Sqr(x - mu) / (2 * sigma * sigma));
};

    //SGGX initialization from Fiber-like distributions
    
    PBRT_CPU_GPU
    void SetMicroflakes(Float sigma,Float angle ){
        //rotate the vector 
        Vector3f initial_vec(0,0,1);        
        Transform m_transform = RotateY(angle);
        std::cout<<"simga "<<sigma << std::endl;
        Vector3f rotate_vec = m_transform(initial_vec);
        S = sggx::Ellipsoid::fromSurface(rotate_vec, sigma);
        //S = sggx::Ellipsoid::fromFiber(t, sigma);

    }

    PBRT_CPU_GPU
    Float p(Vector3f wo, Vector3f wi) const { 
        std::cout<<"here"<<std::endl;
        return sggx::evalPhaseSpecular(wo, wi, S);
    }

    PBRT_CPU_GPU
    Float p(Vector3f wo, Vector3f wi, Point2f u) const { 
        return sggx::evalPhaseDiffuse(wo, wi, S, u.x, u.y);
    }

    PBRT_CPU_GPU
    pstd::optional<PhaseFunctionSample> Sample_p(Vector3f wo, Point2f u) const {

        Vector3f wi = sggx::samplePhaseSpecular(wo, S, u.x, u.y); 
        Float pdf = sggx::evalPhaseSpecular(wo, wi, S);
        return PhaseFunctionSample{pdf, wi, pdf};
    }

    PBRT_CPU_GPU
    pstd::optional<PhaseFunctionSample> Sample_p(Vector3f wo, Point2f u1, Point2f u2) const {
        Vector3f wi = sggx::samplePhaseDiffuse(wo, S, u1.x, u1.y, u2.x, u2.y); 
        Float pdf = sggx::evalPhaseDiffuse(wo, wi, S, u1.x, u1.y);

        //TO DO: Add diffuse SGGX 

        return PhaseFunctionSample{pdf, wi, pdf};
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi) const {std::cout<<"what"<<std::endl; return p(wo, wi); }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, Point2f u) const { return p(wo, wi, u); }

    static const char *Name() { return "SGGX Phase Function"; }

    std::string ToString() const;

  private:
    // SGGXPhaseFunction Private Members
    Float g;
    Float sigma;
    sggx::Ellipsoid S;
    Float m_std;
    Float m_mean;
};

class SGGXPhaseFunctionNormalDistribution {
  public:
    // SGGXFunction Public Methods
    SGGXPhaseFunctionNormalDistribution() = default;
    
    PBRT_CPU_GPU
    SGGXPhaseFunctionNormalDistribution(Float m_rough,Float m_mean,Float m_std):m_rough(m_rough),m_mean(m_mean),m_std(m_std){
        
//        this->m_rough = m_rough;
//        printf("m_rough = %f\n", m_rough);
//        printf("Initialization works");        
    }

    PBRT_CPU_GPU
    SGGXPhaseFunctionNormalDistribution(Float m_mean,Float m_std): m_mean(m_mean),m_std(m_std) {}

Transform RotateX(Float theta) const{
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    SquareMatrix<4> m(1,        0,         0, 0,
                      0, cosTheta, -sinTheta, 0,
                      0, sinTheta,  cosTheta, 0,
                      0,        0,         0, 1);
    return Transform(m, Transpose(m));
}

Transform RotateY(Float theta) const{
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    SquareMatrix<4> m( cosTheta, 0, sinTheta, 0,
                              0, 1,        0, 0,
                      -sinTheta, 0, cosTheta, 0,
                              0, 0,        0, 1);
    return Transform(m, Transpose(m));
}
PBRT_CPU_GPU
sggx::Ellipsoid SampleMicroflakeOrientation(Point2f u,Float m_rough) const{
        Vector3f initial_vec(0,0,1);
        Float sampled_angle = GaussianSample(u,m_mean,m_std);
        // std::cout<<"Sampled angle: "<< sampled_angle << " m_mean " << m_mean << " u " << u << " m_std "<< m_std << std::endl;

        Transform m_transform = RotateY(sampled_angle);
        Transform m_transform2 = RotateX(sampled_angle);
        auto t_transform = m_transform2;
        //better not concatenate the 2 transformations
        // auto t_transform = m_transform2*m_transform;
        Vector3f rotate_vec = t_transform(initial_vec);
        return sggx::Ellipsoid::fromSurface(rotate_vec, m_rough);
    }
inline Float GaussianSample(Point2f u, Float mu = 0, Float m_std = 1) const {
    if (m_std == 0.0){
        return mu;
    }
    else{
        //catmull-box method to transofrm a constant dirstribution to a normal distribution
        Float z1 = std::sqrt(-2 * std::log(u.x)) * std::cos(2*Pi*u.y);
        return z1 *m_std + mu;
        //z2 = std::sqrtf(-2 * std::log(u.x)) * std::sin(2*PI*u.y);
        //account for a non-normalized std
    }
//        return 1 / std::sqrt(2 * Pi * m_rough * m_rough) *
//               std::exp(-Sqr(x - mu) / (2 * m_rough * m_rough));
};


    PBRT_CPU_GPU
    Float p(Vector3f wo, Vector3f wi, Point2f u) const { 
            // std::cout<<"p Specular"<<std::endl;

        sggx::Ellipsoid S = SampleMicroflakeOrientation(u,m_rough);
        return sggx::evalPhaseSpecular(wo, wi, S);  
    }

    PBRT_CPU_GPU
    Float p(Vector3f wo, Vector3f wi, Point2f v, Point2f u) const { 
            // std::cout<<"p Diffuse"<<std::endl;
            sggx::Ellipsoid S = SampleMicroflakeOrientation(u,m_rough);
            return sggx::evalPhaseDiffuse(wo, wi, S, v.x, v.y);
    }

    PBRT_CPU_GPU
    pstd::optional<PhaseFunctionSample> Sample_p(Vector3f wo,Point2f u,Point2f v) const {
        // std::cout<<"Sample Specular"<<std::endl;

        sggx::Ellipsoid S = SampleMicroflakeOrientation(v,m_rough);
        Vector3f wi = sggx::samplePhaseSpecular(wo, S, u.x, u.y); 
        Float pdf = sggx::evalPhaseSpecular(wo, wi, S);
        return PhaseFunctionSample{pdf, wi, pdf};
    }

    PBRT_CPU_GPU
    Float projectedArea(Vector3f wo,Point2f v) const {

        sggx::Ellipsoid S = SampleMicroflakeOrientation(v,m_rough);
        auto res = sggx::projArea(wo,S);
    
        //std::cout<<"res "<< res << " w " << wo<< std::endl;
        return  res;
        // return PhaseFunctionSample{pdf, wi, pdf};
    }

    PBRT_CPU_GPU
    pstd::optional<PhaseFunctionSample> Sample_p(Vector3f wo, Point2f u1, Point2f u2,Point2f u3) const {
        // std::cout<<"Sample Diffuse"<<std::endl;

        sggx::Ellipsoid S = SampleMicroflakeOrientation(u3,m_rough);
        Vector3f wi = sggx::samplePhaseDiffuse(wo, S, u1.x, u1.y, u2.x, u2.y); 
        Float pdf = sggx::evalPhaseDiffuse(wo, wi, S, u1.x, u1.y);

        return PhaseFunctionSample{pdf, wi, pdf};
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi,Point2f v) const { return p(wo, wi,v); }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, Point2f u,Point2f v) const { return p(wo, wi, u,v); }

    static const char *Name() { return "SGGX Phase Function"; }

    std::string ToString() const;

  private:
    // SGGXPhaseFunctionNormalDistribution Private Members
    Float g;
    Float m_rough;
    //sggx::Ellipsoid S;
    Float m_std;
    Float m_mean;
    RNG rng;
    Float max_iter;
};

class SGGXPhaseFunctionMixedRefl {
  public:
    // SGGXFunction Public Methods
    SGGXPhaseFunctionMixedRefl() = default;
    
    PBRT_CPU_GPU
    SGGXPhaseFunctionMixedRefl(Float m_rough,Float m_mean,Float m_std,Float diff_particles_perc,Float g1,Float g2,Float w_g,SampledSpectrum albedo_flake,SampledSpectrum albedo_diff):m_rough(m_rough),
    m_mean(m_mean),m_std(m_std),diff_particles_perc(diff_particles_perc),g1(g1),g2(g2),w_g(w_g),albedo_flake(albedo_flake), albedo_diff(albedo_diff){
        
//        this->m_rough = m_rough;
//        printf("m_rough = %f\n", m_rough);
//        printf("Initialization works");        
    }

    PBRT_CPU_GPU
    SGGXPhaseFunctionMixedRefl(Float m_mean,Float m_std): m_mean(m_mean),m_std(m_std) {}

Transform RotateX(Float theta) const{
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    SquareMatrix<4> m(1,        0,         0, 0,
                      0, cosTheta, -sinTheta, 0,
                      0, sinTheta,  cosTheta, 0,
                      0,        0,         0, 1);
    return Transform(m, Transpose(m));
}

Transform RotateY(Float theta) const{
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    SquareMatrix<4> m( cosTheta, 0, sinTheta, 0,
                              0, 1,        0, 0,
                      -sinTheta, 0, cosTheta, 0,
                              0, 0,        0, 1);
    return Transform(m, Transpose(m));
}
template <typename T>
inline T lin_interp(T a,T b,T w) const{
    if (w <=1 && w >= 0){
        return w *a + (1-w) * b ;
    }else{
        std::cout<<"Err weight is not in [0,1], returning first value only" <<std::endl;
        return a;
    }
    
}
PBRT_CPU_GPU
sggx::Ellipsoid SampleMicroflakeOrientation(Point2f u,Float m_rough) const{
        Vector3f initial_vec(0,0,1);
        Float sampled_angle = GaussianSample(u,m_mean,m_std);
        // std::cout<<"Sampled angle: "<< sampled_angle << " m_mean " << m_mean << " u " << u << " m_std "<< m_std << std::endl;

        Transform m_transform = RotateY(sampled_angle);
        Transform m_transform2 = RotateX(sampled_angle);
        // Transform m_transform3 = Rotatez(sampled_angle);
        auto t_transform = m_transform2;
        //better not concatenate the 2 transformations
        // auto t_transform = m_transform2*m_transform;
        Vector3f rotate_vec = t_transform(initial_vec);
        return sggx::Ellipsoid::fromSurface(rotate_vec, m_rough);
    }
inline Float GaussianSample(Point2f u, Float mu = 0, Float m_std = 1) const {
    if (m_std == 0.0){
        return mu;
    }
    else{
        // std::cout<<"Warning Gaussian Sample is used"<<std::endl;
        //catmull-box method to transofrm a constant dirstribution to a normal distribution
        Float z1 = std::sqrt(-2 * std::log(u.x)) * std::cos(2*Pi*u.y);
        return z1 *m_std + mu;
        //z2 = std::sqrtf(-2 * std::log(u.x)) * std::sin(2*PI*u.y);
        //account for a non-normalized std
    }
//        return 1 / std::sqrt(2 * Pi * m_rough * m_rough) *
//               std::exp(-Sqr(x - mu) / (2 * m_rough * m_rough));
};

    Float diffuserPDF(Vector3f wo, Vector3f wi) const { 
        // std::cout<<"g1 "<<g1 << " g2 "<< g2<< std::endl;
        Float g1_r = HenyeyGreenstein(Dot(wo, wi), g1);
        Float g2_r = HenyeyGreenstein(Dot(wo, wi), g2);
        return lin_interp(g1_r,g2_r,w_g);
    }
    pstd::optional<PhaseFunctionSample> diffuserSample(Vector3f wo, Point2f u1,Float u2) const { 
        Float pdf(0.0f);
        _MediumFlags p_flag = _MediumFlags::Diffuser;
        if (u2 < w_g){
            Vector3f wi = SampleHenyeyGreenstein(wo, g1, u1, &pdf);
            return PhaseFunctionSample{pdf, wi, pdf,p_flag};
        }else{
            // std::cout<<"g2 should not be selected"<<std::endl;
            Vector3f wi = SampleHenyeyGreenstein(wo, g2, u1, &pdf);
            return PhaseFunctionSample{pdf, wi, pdf,p_flag};

        }
    }

        SampledSpectrum albedoW(Vector3f wo, Vector3f wi, Point2f u,Point2f u2) const { 
            // std::cout<<"p Specular"<<std::endl;
        // std::cout<<" wo "<< wo << " wi "<< wi << " u "<< u << std::endl;
        // std::cout<< "m_rough "<< m_rough << std::endl;
        sggx::Ellipsoid S = SampleMicroflakeOrientation(u,m_rough);

        Float microflake = sggx::evalPhaseSpecular(wo, wi, S);
        // Float microflake = sggx::evalPhaseDiffuse(wo, wi, S,u2.x,u2.y);
        Float diffuser = diffuserPDF(wo,wi);
        Float pa = sggx::projArea(wi, S);
        // std::cout<< "diff_particles_perc "<< diff_particles_perc << std::endl;
        // std::cout<< "albedo_flake "<< albedo_flake << std::endl;
        // std::cout<< "albedo_diff "<< albedo_flake << std::endl;
        // auto res = lin_interp(diffuser*albedo_diff,albedo_flake*microflake*pa,diff_particles_perc);
        auto res = (diff_particles_perc)*diffuser*albedo_diff + (1-diff_particles_perc)*(albedo_flake*microflake*pa);
        
        // auto res = lin_interp(diffuser,microflake*pa,diff_particles_perc);
        // auto res = lin_interp(diffuser*pa,microflake,diff_particles_perc);
        //normalize over the area
        auto normalization_factor = diff_particles_perc + (1-diff_particles_perc) * pa;
        return res / normalization_factor;
        
        // lin_interp(1.0f,pa,diff_particles_perc);
        // std::cout << "res "<< res << " microflake "<< microflake << " diff "<<diffuser<< std::endl;
        //return res;
    }

    
    PBRT_CPU_GPU
    Float p(Vector3f wo, Vector3f wi, Point2f u,Point2f u2) const { 
            // std::cout<<"p Specular"<<std::endl;
        // std::cout<<" wo "<< wo << " wi "<< wi << " u "<< u << std::endl;
        // std::cout<< "m_rough "<< m_rough << std::endl;
        sggx::Ellipsoid S = SampleMicroflakeOrientation(u,m_rough);

        Float microflake = sggx::evalPhaseSpecular(wo, wi, S);
        // Float microflake = sggx::evalPhaseDiffuse(wo, wi, S,u2.x,u2.y);

        Float diffuser = diffuserPDF(wo,wi);
        Float pa = sggx::projArea(wi, S);

        auto res = lin_interp(diffuser,microflake*pa,diff_particles_perc);
        // auto res = lin_interp(diffuser,microflake*pa,diff_particles_perc);
        // auto res = lin_interp(diffuser*pa,microflake,diff_particles_perc);
        //normalize over the area
        return res / lin_interp(1.0f,pa,diff_particles_perc);
        // std::cout << "res "<< res << " microflake "<< microflake << " diff "<<diffuser<< std::endl;
        //return res;
    }

    PBRT_CPU_GPU
    /*Float p(Vector3f wo, Vector3f wi, Point2f v, Point2f u) const { 
            // std::cout<<"p Diffuse"<<std::endl;
            sggx::Ellipsoid S = SampleMicroflakeOrientation(u,m_rough);
            return sggx::evalPhaseDiffuse(wo, wi, S, v.x, v.y);
    }*/

    PBRT_CPU_GPU
    pstd::optional<PhaseFunctionSample> Sample_p(Vector3f wo,Point2f u1,Point2f u2,Float u3,Point2f u4) const {
        // std::cout<<"Sample Specular"<<std::endl;
        Float pdf = 0.0;
        
        // std::cout<<"diff_particles_perc "<< diff_particles_perc << std::endl;
    
        if (u3 > diff_particles_perc){
        // std::cout<<"Should not happen"<<std::endl;

        // std::cout<<"u3 "<< u3 << std::endl;
        sggx::Ellipsoid S = SampleMicroflakeOrientation(u1,m_rough);
        Vector3f wi = sggx::samplePhaseSpecular(wo, S, u2.x, u2.y); 
        // Vector3f wi = sggx::samplePhaseDiffuse(wo, S, u2.x, u2.y,u4.x, u4.y); 

        pdf = sggx::evalPhaseSpecular(wo, wi, S);
        _MediumFlags p_flag = _MediumFlags::MicroFlake;
        return PhaseFunctionSample{pdf, wi, pdf,p_flag};

        }else{
        //  std::cout<<"Only happen"<<std::endl;
        
        //Vector3f wi = SampleHenyeyGreenstein(wo, g1, u1, &pdf);
        //return PhaseFunctionSample{pdf, wi, pdf};
        return diffuserSample(wo,u1,u2.x);

        }
    }

    PBRT_CPU_GPU
    Float projectedArea(Vector3f wo,Point2f v) const {

        sggx::Ellipsoid S = SampleMicroflakeOrientation(v,m_rough);
        auto res = sggx::projArea(wo,S);

        auto normalization_factor = diff_particles_perc + (1-diff_particles_perc) * res;
        // std::cout<<"res "<< res << " w " << wo<< std::endl;
        return  normalization_factor;
        // return PhaseFunctionSample{pdf, wi, pdf};
    }

    PBRT_CPU_GPU
    pstd::optional<PhaseFunctionSample> Sample_p(Vector3f wo, Point2f u1, Point2f u2,Point2f u3) const {
        // std::cout<<"Sample Diffuse"<<std::endl;

        sggx::Ellipsoid S = SampleMicroflakeOrientation(u3,m_rough);
        Vector3f wi = sggx::samplePhaseDiffuse(wo, S, u1.x, u1.y, u2.x, u2.y); 
        Float pdf = sggx::evalPhaseDiffuse(wo, wi, S, u1.x, u1.y);

        return PhaseFunctionSample{pdf, wi, pdf};
    }

    PBRT_CPU_GPU
    // Float PDF(Vector3f wo, Vector3f wi,Point2f v) const { return p(wo, wi,v); }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, Point2f u,Point2f v) const { return p(wo, wi, u,v); }

    static const char *Name() { return "SGGX Phase Function"; }

    std::string ToString() const;

  private:
    // SGGXPhaseFunctionMixedRefl Private Members
    Float g1,g2,w_g;
    Float m_rough;
    //sggx::Ellipsoid S;
    Float m_std;
    Float m_mean;
    RNG rng;
    Float max_iter;
    Float diff_particles_perc;
    SampledSpectrum albedo_flake,albedo_diff;
};


class SGGXPhaseFunctionMixedRefl2 {
  public:
    // SGGXFunction Public Methods
    SGGXPhaseFunctionMixedRefl2() = default;
    
    PBRT_CPU_GPU
    SGGXPhaseFunctionMixedRefl2(Float m_rough,Float m_mean,Float std_or,Float std_size,Float diff_particles_perc,Float g):m_rough(m_rough),
    m_mean(m_mean),std_or(std_or),std_size(std_size),diff_particles_perc(diff_particles_perc),g(g){
        
//        this->m_rough = m_rough;
//        printf("m_rough = %f\n", m_rough);
//        printf("Initialization works");        
    }

    // PBRT_CPU_GPU
    // SGGXPhaseFunctionMixedRefl2(Float m_mean,Float m_std): m_mean(m_mean),m_std(m_std) {}

Transform RotateX(Float theta) const{
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    SquareMatrix<4> m(1,        0,         0, 0,
                      0, cosTheta, -sinTheta, 0,
                      0, sinTheta,  cosTheta, 0,
                      0,        0,         0, 1);
    return Transform(m, Transpose(m));
}

Transform RotateY(Float theta) const{
    Float sinTheta = std::sin(Radians(theta));
    Float cosTheta = std::cos(Radians(theta));
    SquareMatrix<4> m( cosTheta, 0, sinTheta, 0,
                              0, 1,        0, 0,
                      -sinTheta, 0, cosTheta, 0,
                              0, 0,        0, 1);
    return Transform(m, Transpose(m));
}
PBRT_CPU_GPU
sggx::Ellipsoid SampleMicroflakeOrientation(Point2f u,Point2f u2,Float m_rough) const{
    
        Vector3f initial_vec(0,0,1);
        auto [sampled_angle,sampled_size] = GaussianSample(u,m_mean,m_rough,std_or,std_size);

        Transform m_transform2 = RotateX(sampled_angle);
        Vector3f rotate_vec = m_transform2(initial_vec);
        return sggx::Ellipsoid::fromSurface(rotate_vec, sampled_size);
    }



inline std::pair<Float,Float> GaussianSample(Point2f u, Float mu_or = 0,Float mu_size = 0, Float m_std_or = 1,Float m_std_size = 1) const {
    
        //catmull-box method to transofrm a constant dirstribution to a normal distribution
        Float z1 = std::sqrt(-2 * std::log(u.x)) * std::cos(2*Pi*u.y);
        Float z2 = std::sqrt(-2 * std::log(u.x)) * std::sin(2*Pi*u.y);
        z1 = z1 *m_std_or + mu_or;
        z2 = z2 *m_std_size + mu_size;

        if (m_std_or == 0.0){
            z1 = mu_or;
        }
        if (m_std_size == 0.0){
            z2 = mu_size;
        }
        return {z1,z2};
};


    PBRT_CPU_GPU


template <typename T>
std::vector<T> linspace(T a, T b, size_t N)const {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

    PBRT_CPU_GPU
    Float p(Vector3f wo, Vector3f wi, Point2f u1, Point2f u2) const { 
            Float res = 0;
            Float total_w = 0;
            Float total_w_or = 0;
            if (std_size == 0.0f && std_or == 0.0f){
                Vector3f initial_vec(0,0,1);   
                Transform m_transform = RotateX(m_mean);
                Vector3f rotate_vec = m_transform(initial_vec);
                sggx::Ellipsoid S = sggx::Ellipsoid::fromSurface(rotate_vec, m_rough);
                // std::cout<<"diff "<< diff_particles_perc << std::endl;
                // return (evalPhaseSpecular(wi,wo,S));
                return diff_particles_perc * (HenyeyGreenstein(Dot(wo, wi), g)) + (1-diff_particles_perc)*(evalPhaseSpecular(wi,wo,S));
                // return sggx::evalPhaseType(diff_particles_perc,wo,wi,S,u1.x,u1.y,u2.x,g);
                
            }
            else{
             return 0.0f;
            }
            /*else{
            Vector3f initial_vec(0,0,1);  
            

            auto [r_or,r_size] =GaussianSample(u1,m_mean,m_rough,std_or,std_size);
            r_size = std::max(0.01f,r_size);

            Float p_g_or(1.0f);
       
            if (std_size != 0 && std_or != 0 ){
            Transform m_transform = RotateX(r_or);
            Vector3f rotated_vec = m_transform(initial_vec);
            sggx::Ellipsoid S = sggx::Ellipsoid::fromSurface(rotated_vec, r_size);
            Float  evalPS = sggx::evalPhaseType(diff_particles_perc,wo,wi,S,u1.x,u1.y,u2.x,g);

            res = evalPS;
            }
            else{
                if (std_size != 0){
                Transform m_transform = RotateX(m_mean);
                Vector3f  rotated_vec = m_transform(initial_vec);
                sggx::Ellipsoid S = sggx::Ellipsoid::fromSurface(rotated_vec, r_size);
                Float  evalPS = sggx::evalPhaseType(diff_particles_perc,wo,wi,S,u1.x,u1.y,u2.x,g);
                res = evalPS;
                }
                else{
                Transform m_transform = RotateX(r_or);
                Vector3f  rotated_vec = m_transform(initial_vec);
                sggx::Ellipsoid S = sggx::Ellipsoid::fromSurface(rotated_vec, m_rough);
                Float  evalPS = sggx::evalPhaseType(diff_particles_perc,wo,wi,S,u1.x,u1.y,u2.x,g);

                res = evalPS ;
                }

            }
            }
            // return res;
            if(IsNaN(res)){
                std::cout<<" nan found "<<std::endl;
            }
            
            return res ;
     */

    }

    

    PBRT_CPU_GPU
    Float projectedArea(Vector3f wo,Point2f u1,Point2f u2) const {

        sggx::Ellipsoid S = SampleMicroflakeOrientation(u1,u2,m_rough);
        auto res = sggx::projArea(wo,S);
    
        //std::cout<<"res "<< res << " w " << wo<< std::endl;
        return  res;
        // return PhaseFunctionSample{pdf, wi, pdf};
    }

    PBRT_CPU_GPU
    pstd::optional<PhaseFunctionSample> Sample_p(Vector3f wo, Point2f u1, Point2f u2,Point2f u3,Float u4) const {
        //randomly sample the orientation
        sggx::Ellipsoid S = SampleMicroflakeOrientation(u1,u2,m_rough);
        //randomly sample the type of reflectance
        Vector3f wi;
        // Vector3f wi = sggx::samplePhaseDiffuse(wo, S, u2.x, u2.y, u3.x, u3.y); 
        Float pdf = 0.0f;
        if (u4 < diff_particles_perc){
        // wi = sggx::samplePhaseDiffuse(wo, S, u2.x, u2.y, u3.x, u3.y);             
        // pdf = sggx::evalPhaseDiffuse(wo, wi, S, u1.x, u1.y);
        // std::cout<<"should not happen"<<std::endl;
        Vector3f wi = SampleHenyeyGreenstein(wo, g, u1, &pdf);
        return PhaseFunctionSample{pdf, wi, pdf};
        
        }else{
            
        wi = sggx:: samplePhaseSpecular(wi, S, u3.x, u3.y);
        pdf = sggx::evalPhaseSpecular(wi,wo,S);
        
        // Float pdf = sggx::evalPhaseDiffuse(wo, wi, S, u1.x, u1.y);
        
        
        // Vector3f wi = sggx::samplePhaseSpecular(wo, S, u3.x, u3.y);
        // Float pdf = evalPhaseSpecular(wo, wi, S); 
        //         return sggx::evalPhaseType(diff_particles_perc,wo,wi,S,u1.x,u1.y,u2.x);
        
        // Float pdf = sggx::evalPhaseDiffuse(wo, wi, S, u1.x, u1.y);
        }
        return PhaseFunctionSample{pdf, wi, pdf};
    }

    // PBRT_CPU_GPU
    // Float PDF(Vector3f wo, Vector3f wi,Point2f v) const { return p(wo, wi,v); }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, Point2f u,Point2f v) const { return p(wo, wi, u,v); }

    static const char *Name() { return "SGGX Phase Function"; }

    std::string ToString() const;

  private:
    // SGGXPhaseFunctionMixedRefl2 Private Members
    Float g;
    Float m_rough;
    //sggx::Ellipsoid S;
    Float std_or;
    Float std_size;
    Float m_mean;
    Float diff_particles_perc;
    RNG rng;
    Float max_iter;
};



// MediumProperties Definition
struct MediumProperties {
    SampledSpectrum sigma_a, sigma_s;
    PhaseFunction phase;
    SampledSpectrum Le;
};

// HomogeneousMajorantIterator Definition
class HomogeneousMajorantIterator {
  public:
    // HomogeneousMajorantIterator Public Methods
    PBRT_CPU_GPU
    HomogeneousMajorantIterator() : called(true) {}
    PBRT_CPU_GPU
    HomogeneousMajorantIterator(Float tMin, Float tMax, SampledSpectrum sigma_maj)
        : seg{tMin, tMax, sigma_maj}, called(false) {}

    PBRT_CPU_GPU
    pstd::optional<RayMajorantSegment> Next() {
        if (called)
            return {};
        called = true;
        return seg;
    }

    std::string ToString() const;

  private:
    RayMajorantSegment seg;
    bool called;
};

// MajorantGrid Definition
struct MajorantGrid {
    // MajorantGrid Public Methods
    MajorantGrid() = default;
    MajorantGrid(Bounds3f bounds, Point3i res, Allocator alloc)
        : bounds(bounds), voxels(res.x * res.y * res.z, alloc), res(res) {}

    PBRT_CPU_GPU
    Float Lookup(int x, int y, int z) const {
        DCHECK(x >= 0 && x < res.x && y >= 0 && y < res.y && z >= 0 && z < res.z);
        return voxels[x + res.x * (y + res.y * z)];
    }
    PBRT_CPU_GPU
    void Set(int x, int y, int z, Float v) {
        DCHECK(x >= 0 && x < res.x && y >= 0 && y < res.y && z >= 0 && z < res.z);
        voxels[x + res.x * (y + res.y * z)] = v;
    }

    PBRT_CPU_GPU
    Bounds3f VoxelBounds(int x, int y, int z) const {
        Point3f p0(Float(x) / res.x, Float(y) / res.y, Float(z) / res.z);
        Point3f p1(Float(x + 1) / res.x, Float(y + 1) / res.y, Float(z + 1) / res.z);
        return Bounds3f(p0, p1);
    }

    // MajorantGrid Public Members
    Bounds3f bounds;
    pstd::vector<Float> voxels;
    Point3i res;
};

// DDAMajorantIterator Definition
class DDAMajorantIterator {
  public:
    // DDAMajorantIterator Public Methods
    DDAMajorantIterator() = default;
    PBRT_CPU_GPU
    DDAMajorantIterator(Ray ray, Float tMin, Float tMax, const MajorantGrid *grid,
                        SampledSpectrum sigma_t)
        : tMin(tMin), tMax(tMax), grid(grid), sigma_t(sigma_t) {
        // Set up 3D DDA for ray through the majorant grid
        Vector3f diag = grid->bounds.Diagonal();
        Ray rayGrid(Point3f(grid->bounds.Offset(ray.o)),
                    Vector3f(ray.d.x / diag.x, ray.d.y / diag.y, ray.d.z / diag.z));
        Point3f gridIntersect = rayGrid(tMin);
        for (int axis = 0; axis < 3; ++axis) {
            // Initialize ray stepping parameters for _axis_
            // Compute current voxel for axis and handle negative zero direction
            voxel[axis] =
                Clamp(gridIntersect[axis] * grid->res[axis], 0, grid->res[axis] - 1);
            deltaT[axis] = 1 / (std::abs(rayGrid.d[axis]) * grid->res[axis]);
            if (rayGrid.d[axis] == -0.f)
                rayGrid.d[axis] = 0.f;

            if (rayGrid.d[axis] >= 0) {
                // Handle ray with positive direction for voxel stepping
                Float nextVoxelPos = Float(voxel[axis] + 1) / grid->res[axis];
                nextCrossingT[axis] =
                    tMin + (nextVoxelPos - gridIntersect[axis]) / rayGrid.d[axis];
                step[axis] = 1;
                voxelLimit[axis] = grid->res[axis];

            } else {
                // Handle ray with negative direction for voxel stepping
                Float nextVoxelPos = Float(voxel[axis]) / grid->res[axis];
                nextCrossingT[axis] =
                    tMin + (nextVoxelPos - gridIntersect[axis]) / rayGrid.d[axis];
                step[axis] = -1;
                voxelLimit[axis] = -1;
            }
        }
    }

    PBRT_CPU_GPU
    pstd::optional<RayMajorantSegment> Next() {
        if (tMin >= tMax)
            return {};
        // Find _stepAxis_ for stepping to next voxel and exit point _tVoxelExit_
        int bits = ((nextCrossingT[0] < nextCrossingT[1]) << 2) +
                   ((nextCrossingT[0] < nextCrossingT[2]) << 1) +
                   ((nextCrossingT[1] < nextCrossingT[2]));
        const int cmpToAxis[8] = {2, 1, 2, 1, 2, 2, 0, 0};
        int stepAxis = cmpToAxis[bits];
        Float tVoxelExit = std::min(tMax, nextCrossingT[stepAxis]);

        // Get _maxDensity_ for current voxel and initialize _RayMajorantSegment_, _seg_
        SampledSpectrum sigma_maj = sigma_t * grid->Lookup(voxel[0], voxel[1], voxel[2]);
        RayMajorantSegment seg{tMin, tVoxelExit, sigma_maj};

        // Advance to next voxel in maximum density grid
        tMin = tVoxelExit;
        if (nextCrossingT[stepAxis] > tMax)
            tMin = tMax;
        voxel[stepAxis] += step[stepAxis];
        if (voxel[stepAxis] == voxelLimit[stepAxis])
            tMin = tMax;
        nextCrossingT[stepAxis] += deltaT[stepAxis];

        return seg;
    }

    std::string ToString() const;

  private:
    // DDAMajorantIterator Private Members
    SampledSpectrum sigma_t;
    Float tMin = Infinity, tMax = -Infinity;
    const MajorantGrid *grid;
    Float nextCrossingT[3], deltaT[3];
    int step[3], voxelLimit[3], voxel[3];
};

// HomogeneousMedium Definition
class HomogeneousMedium {
  public:
    // HomogeneousMedium Public Type Definitions
    using MajorantIterator = HomogeneousMajorantIterator;

    // HomogeneousMedium Public Methods
    HomogeneousMedium(Spectrum sigma_a, Spectrum sigma_s, Float sigmaScale, Spectrum Le,
                      Float LeScale, Float g, Allocator alloc)
        : sigma_a_spec(sigma_a, alloc),
          sigma_s_spec(sigma_s, alloc),
          Le_spec(Le, alloc),
          phase(g) {
        sigma_a_spec.Scale(sigmaScale);
        sigma_s_spec.Scale(sigmaScale);
        Le_spec.Scale(LeScale);
    }

    static HomogeneousMedium *Create(const ParameterDictionary &parameters,
                                     const FileLoc *loc, Allocator alloc);

    PBRT_CPU_GPU
    bool IsEmissive() const { return Le_spec.MaxValue() > 0; }

    PBRT_CPU_GPU
    MediumProperties SamplePoint(Point3f p, const SampledWavelengths &lambda) const {
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);
        SampledSpectrum Le = Le_spec.Sample(lambda);
        return MediumProperties{sigma_a, sigma_s, &phase, Le};
    }

    PBRT_CPU_GPU
    HomogeneousMajorantIterator SampleRay(Ray ray, Float tMax,
                                          const SampledWavelengths &lambda) const {
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);
        return HomogeneousMajorantIterator(0, tMax, sigma_a + sigma_s);
    }

    std::string ToString() const;

  private:
    // HomogeneousMedium Private Data
    DenselySampledSpectrum sigma_a_spec, sigma_s_spec, Le_spec;
    HGPhaseFunction phase;
};


class HomogeneousMediumSGGX {
  public:
    // HomogeneousMedium Public Type Definitions
    using MajorantIterator = HomogeneousMajorantIterator;
    std::string ToString() const;
    // HomogeneousMedium Public Methods
    HomogeneousMediumSGGX(Spectrum sigma_a, Spectrum sigma_s, Float sigmaScale, Spectrum Le,
                      Float LeScale, Float g, Vector3f dir,Allocator alloc)
        : sigma_a_spec(sigma_a, alloc),
          sigma_s_spec(sigma_s, alloc),
          Le_spec(Le, alloc),
          phase(g,dir)
          {
        sigma_a_spec.Scale(sigmaScale);
        sigma_s_spec.Scale(sigmaScale);
        Le_spec.Scale(LeScale);
        
    }

    static HomogeneousMediumSGGX *Create(const ParameterDictionary &parameters,
                                     const FileLoc *loc, Allocator alloc);

    PBRT_CPU_GPU
    bool IsEmissive() const { return Le_spec.MaxValue() > 0; }

    PBRT_CPU_GPU
    MediumProperties SamplePoint(Point3f p, const SampledWavelengths &lambda) const {
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);
        SampledSpectrum Le = Le_spec.Sample(lambda);
        return MediumProperties{sigma_a, sigma_s, &phase, Le};
    }

    PBRT_CPU_GPU
    HomogeneousMajorantIterator SampleRay(Ray ray, Float tMax,
                                          const SampledWavelengths &lambda) const {
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);
        return HomogeneousMajorantIterator(0, tMax, sigma_a + sigma_s);
    }

    

  private:
    // HomogeneousMedium Private Data
    DenselySampledSpectrum sigma_a_spec, sigma_s_spec, Le_spec;
    //HGPhaseFunction phase;
    SGGXPhaseFunction phase;
};


// GridMedium Definition
class GridMedium {
  public:
    // GridMedium Public Type Definitions
    using MajorantIterator = DDAMajorantIterator;

    // GridMedium Public Methods
    GridMedium(const Bounds3f &bounds, const Transform &renderFromMedium,
               Spectrum sigma_a, Spectrum sigma_s, Float sigmaScale, Float g,
               SampledGrid<Float> density, pstd::optional<SampledGrid<Float>> temperature,
               Float temperatureScale, Float temperatureOffset,
               Spectrum Le, SampledGrid<Float> LeScale, Allocator alloc);

    static GridMedium *Create(const ParameterDictionary &parameters,
                              const Transform &renderFromMedium, const FileLoc *loc,
                              Allocator alloc);

    std::string ToString() const;

    PBRT_CPU_GPU
    bool IsEmissive() const { return isEmissive; }

    PBRT_CPU_GPU
    MediumProperties SamplePoint(Point3f p, const SampledWavelengths &lambda) const {
        // Sample spectra for grid medium $\sigmaa$ and $\sigmas$
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);

        // Scale scattering coefficients by medium density at _p_
        p = renderFromMedium.ApplyInverse(p);
        p = Point3f(bounds.Offset(p));
        Float d = densityGrid.Lookup(p);
        sigma_a *= d;
        sigma_s *= d;

        // Compute grid emission _Le_ at _p_
        SampledSpectrum Le(0.f);
        if (isEmissive) {
            Float scale = LeScale.Lookup(p);
            if (scale > 0) {
                // Compute emitted radiance using _temperatureGrid_ or _Le_spec_
                if (temperatureGrid) {
                    Float temp = temperatureGrid->Lookup(p);
                    // Added after book publication: optionally offset and scale
                    // temperature based on user-supplied parameters. (Match
                    // NanoVDBMedium functionality.)
                    temp = (temp - temperatureOffset) * temperatureScale;
                    if (temp > 100.f)
                        Le = scale * BlackbodySpectrum(temp).Sample(lambda);
                } else
                    Le = scale * Le_spec.Sample(lambda);
            }
        }

        return MediumProperties{sigma_a, sigma_s, &phase, Le};
    }

    PBRT_CPU_GPU
    DDAMajorantIterator SampleRay(Ray ray, Float raytMax,
                                  const SampledWavelengths &lambda) const {
        // Transform ray to medium's space and compute bounds overlap
        ray = renderFromMedium.ApplyInverse(ray, &raytMax);
        Float tMin, tMax;
        if (!bounds.IntersectP(ray.o, ray.d, raytMax, &tMin, &tMax))
            return {};
        DCHECK_LE(tMax, raytMax);

        // Sample spectra for grid medium $\sigmaa$ and $\sigmas$
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);

        SampledSpectrum sigma_t = sigma_a + sigma_s;
        return DDAMajorantIterator(ray, tMin, tMax, &majorantGrid, sigma_t);
    }

  private:
    // GridMedium Private Members
    Bounds3f bounds;
    Transform renderFromMedium;
    DenselySampledSpectrum sigma_a_spec, sigma_s_spec;
    SampledGrid<Float> densityGrid;
    HGPhaseFunction phase;
    pstd::optional<SampledGrid<Float>> temperatureGrid;
    DenselySampledSpectrum Le_spec;
    SampledGrid<Float> LeScale;
    bool isEmissive;
    Float temperatureScale, temperatureOffset;
    MajorantGrid majorantGrid;
};

// RGBGridMedium Definition
class RGBGridMedium {
  public:
    // RGBGridMedium Public Type Definitions
    using MajorantIterator = DDAMajorantIterator;

    // RGBGridMedium Public Methods
    RGBGridMedium(const Bounds3f &bounds, const Transform &renderFromMedium, Float g,
                  pstd::optional<SampledGrid<RGBUnboundedSpectrum>> sigma_a,
                  pstd::optional<SampledGrid<RGBUnboundedSpectrum>> sigma_s,
                  Float sigmaScale, pstd::optional<SampledGrid<RGBIlluminantSpectrum>> Le,
                  Float LeScale, Allocator alloc);

    static RGBGridMedium *Create(const ParameterDictionary &parameters,
                                 const Transform &renderFromMedium, const FileLoc *loc,
                                 Allocator alloc);

    std::string ToString() const;

    PBRT_CPU_GPU
    bool IsEmissive() const { return LeGrid && LeScale > 0; }

    PBRT_CPU_GPU
    MediumProperties SamplePoint(Point3f p, const SampledWavelengths &lambda) const {
        p = renderFromMedium.ApplyInverse(p);
        p = Point3f(bounds.Offset(p));
        // Compute $\sigmaa$ and $\sigmas$ for _RGBGridMedium_
        auto convert = [=] PBRT_CPU_GPU(RGBUnboundedSpectrum s) {
            return s.Sample(lambda);
        };
        SampledSpectrum sigma_a =
            sigmaScale *
            (sigma_aGrid ? sigma_aGrid->Lookup(p, convert) : SampledSpectrum(1.f));
        SampledSpectrum sigma_s =
            sigmaScale *
            (sigma_sGrid ? sigma_sGrid->Lookup(p, convert) : SampledSpectrum(1.f));

        // Find emitted radiance _Le_ for _RGBGridMedium_
        SampledSpectrum Le(0.f);
        if (LeGrid && LeScale > 0) {
            auto convert = [=] PBRT_CPU_GPU(RGBIlluminantSpectrum s) {
                return s.Sample(lambda);
            };
            Le = LeScale * LeGrid->Lookup(p, convert);
        }

        return MediumProperties{sigma_a, sigma_s, &phase, Le};
    }

    PBRT_CPU_GPU
    DDAMajorantIterator SampleRay(Ray ray, Float raytMax,
                                  const SampledWavelengths &lambda) const {
        // Transform ray to medium's space and compute bounds overlap
        ray = renderFromMedium.ApplyInverse(ray, &raytMax);
        Float tMin, tMax;
        if (!bounds.IntersectP(ray.o, ray.d, raytMax, &tMin, &tMax))
            return {};
        DCHECK_LE(tMax, raytMax);

        SampledSpectrum sigma_t(1);
        return DDAMajorantIterator(ray, tMin, tMax, &majorantGrid, sigma_t);
    }

  private:
    // RGBGridMedium Private Members
    Bounds3f bounds;
    Transform renderFromMedium;
    pstd::optional<SampledGrid<RGBIlluminantSpectrum>> LeGrid;
    Float LeScale;
    HGPhaseFunction phase;
    pstd::optional<SampledGrid<RGBUnboundedSpectrum>> sigma_aGrid, sigma_sGrid;
    Float sigmaScale;
    MajorantGrid majorantGrid;
};

// CloudMedium Definition
class CloudMedium {
  public:
    // CloudMedium Public Type Definitions
    using MajorantIterator = HomogeneousMajorantIterator;

    // CloudMedium Public Methods
    static CloudMedium *Create(const ParameterDictionary &parameters,
                               const Transform &renderFromMedium, const FileLoc *loc,
                               Allocator alloc);

    std::string ToString() const {
        return StringPrintf("[ CloudMedium bounds: %s renderFromMedium: %s phase: %s "
                            "sigma_a_spec: %s sigma_s_spec: %s density: %f wispiness: %f "
                            "frequency: %f ]",
                            bounds, renderFromMedium, phase, sigma_a_spec, sigma_s_spec,
                            density, wispiness, frequency);
    }

    CloudMedium(const Bounds3f &bounds, const Transform &renderFromMedium,
                Spectrum sigma_a, Spectrum sigma_s, Float g, Float density,
                Float wispiness, Float frequency, Allocator alloc)
        : bounds(bounds),
          renderFromMedium(renderFromMedium),
          sigma_a_spec(sigma_a, alloc),
          sigma_s_spec(sigma_s, alloc),
          phase(g),
          density(density),
          wispiness(wispiness),
          frequency(frequency) {}

    PBRT_CPU_GPU
    bool IsEmissive() const { return false; }

    PBRT_CPU_GPU
    MediumProperties SamplePoint(Point3f p, const SampledWavelengths &lambda) const {
        // Compute sampled spectra for cloud $\sigmaa$ and $\sigmas$ at _p_
        Float density = Density(renderFromMedium.ApplyInverse(p));
        SampledSpectrum sigma_a = density * sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = density * sigma_s_spec.Sample(lambda);

        return MediumProperties{sigma_a, sigma_s, &phase, SampledSpectrum(0.f)};
    }

    PBRT_CPU_GPU
    HomogeneousMajorantIterator SampleRay(Ray ray, Float raytMax,
                                          const SampledWavelengths &lambda) const {
        // Transform ray to medium's space and compute bounds overlap
        ray = renderFromMedium.ApplyInverse(ray, &raytMax);
        Float tMin, tMax;
        if (!bounds.IntersectP(ray.o, ray.d, raytMax, &tMin, &tMax))
            return {};
        DCHECK_LE(tMax, raytMax);

        // Compute $\sigmat$ bound for cloud medium and initialize majorant iterator
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);
        SampledSpectrum sigma_t = sigma_a + sigma_s;
        return HomogeneousMajorantIterator(tMin, tMax, sigma_t);
    }

  private:
    // CloudMedium Private Methods
    PBRT_CPU_GPU
    Float Density(Point3f p) const {
        Point3f pp = frequency * p;
        if (wispiness > 0) {
            // Perturb cloud lookup point _pp_ using noise
            Float vomega = 0.05f * wispiness, vlambda = 10.f;
            for (int i = 0; i < 2; ++i) {
                pp += vomega * DNoise(vlambda * pp);
                vomega *= 0.5f;
                vlambda *= 1.99f;
            }
        }
        // Sum scales of noise to approximate cloud density
        Float d = 0;
        Float omega = 0.5f, lambda = 1.f;
        for (int i = 0; i < 5; ++i) {
            d += omega * Noise(lambda * pp);
            omega *= 0.5f;
            lambda *= 1.99f;
        }

        // Model decrease in density with altitude and return final cloud density
        d = Clamp((1 - p.y) * 4.5f * density * d, 0, 1);
        d += 2 * std::max<Float>(0, 0.5f - p.y);
        return Clamp(d, 0, 1);
    }

    // CloudMedium Private Members
    Bounds3f bounds;
    Transform renderFromMedium;
    HGPhaseFunction phase;
    DenselySampledSpectrum sigma_a_spec, sigma_s_spec;
    Float density, wispiness, frequency;
};

// NanoVDBMedium Definition
// NanoVDBBuffer Definition
class NanoVDBBuffer {
  public:
    static inline void ptrAssert(void *ptr, const char *msg, const char *file, int line,
                                 bool abort = true) {
        if (abort)
            LOG_FATAL("%p: %s (%s:%d)", ptr, msg, file, line);
        else
            LOG_ERROR("%p: %s (%s:%d)", ptr, msg, file, line);
    }

    NanoVDBBuffer() = default;
    NanoVDBBuffer(Allocator alloc) : alloc(alloc) {}
    NanoVDBBuffer(size_t size, Allocator alloc = {}) : alloc(alloc) { init(size); }
    NanoVDBBuffer(const NanoVDBBuffer &) = delete;
    NanoVDBBuffer(NanoVDBBuffer &&other) noexcept
        : alloc(std::move(other.alloc)),
          bytesAllocated(other.bytesAllocated),
          ptr(other.ptr) {
        other.bytesAllocated = 0;
        other.ptr = nullptr;
    }
    NanoVDBBuffer &operator=(const NanoVDBBuffer &) = delete;
    NanoVDBBuffer &operator=(NanoVDBBuffer &&other) noexcept {
        // Note, this isn't how std containers work, but it's expedient for
        // our purposes here...
        clear();
        // operator= was deleted? Fine.
        new (&alloc) Allocator(other.alloc.resource());
        bytesAllocated = other.bytesAllocated;
        ptr = other.ptr;
        other.bytesAllocated = 0;
        other.ptr = nullptr;
        return *this;
    }
    ~NanoVDBBuffer() { clear(); }

    void init(uint64_t size) {
        if (size == bytesAllocated)
            return;
        if (bytesAllocated > 0)
            clear();
        if (size == 0)
            return;
        bytesAllocated = size;
        ptr = (uint8_t *)alloc.allocate_bytes(bytesAllocated, 128);
    }

    const uint8_t *data() const { return ptr; }
    uint8_t *data() { return ptr; }
    uint64_t size() const { return bytesAllocated; }
    bool empty() const { return size() == 0; }

    void clear() {
        alloc.deallocate_bytes(ptr, bytesAllocated, 128);
        bytesAllocated = 0;
        ptr = nullptr;
    }

    static NanoVDBBuffer create(uint64_t size, const NanoVDBBuffer *context = nullptr) {
        return NanoVDBBuffer(size, context ? context->GetAllocator() : Allocator());
    }

    Allocator GetAllocator() const { return alloc; }

  private:
    Allocator alloc;
    size_t bytesAllocated = 0;
    uint8_t *ptr = nullptr;
};

class NanoVDBMedium {
  public:
    using MajorantIterator = DDAMajorantIterator;
    // NanoVDBMedium Public Methods
    static NanoVDBMedium *Create(const ParameterDictionary &parameters,
                                 const Transform &renderFromMedium, const FileLoc *loc,
                                 Allocator alloc);

    std::string ToString() const;

    NanoVDBMedium(const Transform &renderFromMedium, Spectrum sigma_a, Spectrum sigma_s,
                  Float sigmaScale, Float g, nanovdb::GridHandle<NanoVDBBuffer> dg,
                  nanovdb::GridHandle<NanoVDBBuffer> tg, Float LeScale,
                  Float temperatureOffset, Float temperatureScale, Allocator alloc);

    PBRT_CPU_GPU
    bool IsEmissive() const { return temperatureFloatGrid && LeScale > 0; }

    PBRT_CPU_GPU
    MediumProperties SamplePoint(Point3f p, const SampledWavelengths &lambda) const {
        
        // Sample spectra for grid $\sigmaa$ and $\sigmas$
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);

        // Scale scattering coefficients by medium density at _p_
        p = renderFromMedium.ApplyInverse(p);

        nanovdb::Vec3<float> pIndex =
            densityFloatGrid->worldToIndexF(nanovdb::Vec3<float>(p.x, p.y, p.z));
        using Sampler = nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
        Float d = Sampler(densityFloatGrid->tree())(pIndex);

        return MediumProperties{sigma_a * d, sigma_s * d, &phase, Le(p, lambda)};
    }

    PBRT_CPU_GPU
    DDAMajorantIterator SampleRay(Ray ray, Float raytMax,
                                  const SampledWavelengths &lambda) const {
        // Transform ray to medium's space and compute bounds overlap
        ray = renderFromMedium.ApplyInverse(ray, &raytMax);
        Float tMin, tMax;
        if (!bounds.IntersectP(ray.o, ray.d, raytMax, &tMin, &tMax))
            return {};
        DCHECK_LE(tMax, raytMax);

        // Sample spectra for grid $\sigmaa$ and $\sigmas$
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);

        SampledSpectrum sigma_t = sigma_a + sigma_s;
        return DDAMajorantIterator(ray, tMin, tMax, &majorantGrid, sigma_t);
    }

  private:
    // NanoVDBMedium Private Methods
    PBRT_CPU_GPU
    SampledSpectrum Le(Point3f p, const SampledWavelengths &lambda) const {
        if (!temperatureFloatGrid)
            return SampledSpectrum(0.f);
        nanovdb::Vec3<float> pIndex =
            temperatureFloatGrid->worldToIndexF(nanovdb::Vec3<float>(p.x, p.y, p.z));
        using Sampler = nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
        Float temp = Sampler(temperatureFloatGrid->tree())(pIndex);
        temp = (temp - temperatureOffset) * temperatureScale;
        if (temp <= 100.f)
            return SampledSpectrum(0.f);
        return LeScale * BlackbodySpectrum(temp).Sample(lambda);
    }

    // NanoVDBMedium Private Members
    Bounds3f bounds;
    Transform renderFromMedium;
    DenselySampledSpectrum sigma_a_spec, sigma_s_spec;
    HGPhaseFunction phase;
    MajorantGrid majorantGrid;
    nanovdb::GridHandle<NanoVDBBuffer> densityGrid;
    nanovdb::GridHandle<NanoVDBBuffer> temperatureGrid;
    const nanovdb::FloatGrid *densityFloatGrid = nullptr;
    const nanovdb::FloatGrid *temperatureFloatGrid = nullptr;
    Float LeScale, temperatureOffset, temperatureScale;
};

class SGGXMedium {
  public:
    using MajorantIterator = DDAMajorantIterator;
    // SGGXMedium Public Methods
    static SGGXMedium *Create(const ParameterDictionary &parameters,
                                 const Transform &renderFromMedium, const FileLoc *loc,
                                 Allocator alloc);

    std::string ToString() const;

    SGGXMedium(const Transform &renderFromMedium, Spectrum sigma_a, Spectrum sigma_s,
                  Float sigmaScale, Float g, nanovdb::GridHandle<NanoVDBBuffer> dg,
                  nanovdb::GridHandle<NanoVDBBuffer> tg, Float LeScale,
                  Float temperatureOffset, Float temperatureScale, Allocator alloc);

    PBRT_CPU_GPU
    bool IsEmissive() const { return temperatureFloatGrid && LeScale > 0; }

    PBRT_CPU_GPU
    MediumProperties SamplePoint(Point3f p, const SampledWavelengths &lambda) const {
        
        //printf("We are using SGGX NanoVDBMedium\n");

        // Sample spectra for grid $\sigmaa$ and $\sigmas$
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);

        // Scale scattering coefficients by medium density at _p_
        p = renderFromMedium.ApplyInverse(p);

        nanovdb::Vec3<float> pIndex =
            densityFloatGrid->worldToIndexF(nanovdb::Vec3<float>(p.x, p.y, p.z));
        using Sampler = nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
        Float d = Sampler(densityFloatGrid->tree())(pIndex);

        return MediumProperties{sigma_a * d, sigma_s * d, &phase, Le(p, lambda)};
    }

    PBRT_CPU_GPU
    DDAMajorantIterator SampleRay(Ray ray, Float raytMax,
                                  const SampledWavelengths &lambda) const {
        // Transform ray to medium's space and compute bounds overlap
        ray = renderFromMedium.ApplyInverse(ray, &raytMax);
        Float tMin, tMax;
        if (!bounds.IntersectP(ray.o, ray.d, raytMax, &tMin, &tMax))
            return {};
        DCHECK_LE(tMax, raytMax);

        // Sample spectra for grid $\sigmaa$ and $\sigmas$
        SampledSpectrum sigma_a = sigma_a_spec.Sample(lambda);
        SampledSpectrum sigma_s = sigma_s_spec.Sample(lambda);

        SampledSpectrum sigma_t = sigma_a + sigma_s;
        return DDAMajorantIterator(ray, tMin, tMax, &majorantGrid, sigma_t);
    }

  private:
    // SGGXMedium Private Methods
    PBRT_CPU_GPU
    SampledSpectrum Le(Point3f p, const SampledWavelengths &lambda) const {
        if (!temperatureFloatGrid)
            return SampledSpectrum(0.f);
        nanovdb::Vec3<float> pIndex =
            temperatureFloatGrid->worldToIndexF(nanovdb::Vec3<float>(p.x, p.y, p.z));
        using Sampler = nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;
        Float temp = Sampler(temperatureFloatGrid->tree())(pIndex);
        temp = (temp - temperatureOffset) * temperatureScale;
        if (temp <= 100.f)
            return SampledSpectrum(0.f);
        return LeScale * BlackbodySpectrum(temp).Sample(lambda);
    }

    // SGGXMedium Private Members
    Bounds3f bounds;
    Transform renderFromMedium;
    DenselySampledSpectrum sigma_a_spec, sigma_s_spec;
    //HGPhaseFunction phase;
    SGGXPhaseFunction phase;
    MajorantGrid majorantGrid;
    nanovdb::GridHandle<NanoVDBBuffer> densityGrid;
    nanovdb::GridHandle<NanoVDBBuffer> temperatureGrid;
    const nanovdb::FloatGrid *densityFloatGrid = nullptr;
    const nanovdb::FloatGrid *temperatureFloatGrid = nullptr;
    Float LeScale, temperatureOffset, temperatureScale;
};




inline Float PhaseFunction::p(Vector3f wo, Vector3f wi) const {
    auto p = [&](auto ptr) { return ptr->p(wo, wi); };
    return Dispatch(p);
}

inline pstd::optional<PhaseFunctionSample> PhaseFunction::Sample_p(Vector3f wo,
                                                                   Point2f u) const {
    auto sample = [&](auto ptr) { return ptr->Sample_p(wo, u); };
    return Dispatch(sample);
}

inline Float PhaseFunction::PDF(Vector3f wo, Vector3f wi) const {
    auto pdf = [&](auto ptr) { return ptr->PDF(wo, wi); };
    return Dispatch(pdf);
}

inline pstd::optional<RayMajorantSegment> RayMajorantIterator::Next() {
    auto next = [](auto ptr) { return ptr->Next(); };
    return Dispatch(next);
}

inline MediumProperties Medium::SamplePoint(Point3f p,
                                            const SampledWavelengths &lambda) const {
    auto sample = [&](auto ptr) { return ptr->SamplePoint(p, lambda); };
    return Dispatch(sample);
}

// Medium Sampling Function Definitions
inline RayMajorantIterator Medium::SampleRay(Ray ray, Float tMax,
                                             const SampledWavelengths &lambda,
                                             ScratchBuffer &buf) const {
    // Explicit capture to work around MSVC weirdness; it doesn't see |buf| otherwise...
    auto sample = [ray, tMax, lambda, &buf](auto medium) {
        // Return _RayMajorantIterator_ for medium's majorant iterator
        using ConcreteMedium = typename std::remove_reference_t<decltype(*medium)>;
        using Iter = typename ConcreteMedium::MajorantIterator;
        Iter *iter = (Iter *)buf.Alloc(sizeof(Iter), alignof(Iter));
        *iter = medium->SampleRay(ray, tMax, lambda);
        return RayMajorantIterator(iter);
    };
    return DispatchCPU(sample);
}

template <typename F>
PBRT_CPU_GPU SampledSpectrum SampleT_maj(Ray ray, Float tMax, Float u, RNG &rng,
                                         const SampledWavelengths &lambda, F callback) {
    auto sample = [&](auto medium) {
        using M = typename std::remove_reference_t<decltype(*medium)>;
        return SampleT_maj<M>(ray, tMax, u, rng, lambda, callback);
    };
    return ray.medium.Dispatch(sample);
}

template <typename ConcreteMedium, typename F>
PBRT_CPU_GPU SampledSpectrum SampleT_maj(Ray ray, Float tMax, Float u, RNG &rng,
                                         const SampledWavelengths &lambda, F callback) {
    // Normalize ray direction and update _tMax_ accordingly
    tMax *= Length(ray.d);
    ray.d = Normalize(ray.d);

    // Initialize _MajorantIterator_ for ray majorant sampling
    ConcreteMedium *medium = ray.medium.Cast<ConcreteMedium>();
    typename ConcreteMedium::MajorantIterator iter = medium->SampleRay(ray, tMax, lambda);

    // Generate ray majorant samples until termination
    SampledSpectrum T_maj(1.f);
    bool done = false;
    while (!done) {
        // Get next majorant segment from iterator and sample it
        pstd::optional<RayMajorantSegment> seg = iter.Next();
        if (!seg)
            return T_maj;
        // Handle zero-valued majorant for current segment
        if (seg->sigma_maj[0] == 0) {
            Float dt = seg->tMax - seg->tMin;
            // Handle infinite _dt_ for ray majorant segment
            if (IsInf(dt))
                dt = std::numeric_limits<Float>::max();

            T_maj *= FastExp(-dt * seg->sigma_maj);
            continue;
        }

        // Generate samples along current majorant segment
        Float tMin = seg->tMin;
        while (true) {
            // Try to generate sample along current majorant segment
            Float t = tMin + SampleExponential(u, seg->sigma_maj[0]);
            PBRT_DBG("Sampled t = %f from tMin %f u %f sigma_maj[0] %f\n", t, tMin, u,
                     seg->sigma_maj[0]);
            u = rng.Uniform<Float>();
            if (t < seg->tMax) {
                // Call callback function for sample within segment
                PBRT_DBG("t < seg->tMax\n");
                T_maj *= FastExp(-(t - tMin) * seg->sigma_maj);
                MediumProperties mp = medium->SamplePoint(ray(t), lambda);
                if (!callback(ray(t), mp, seg->sigma_maj, T_maj)) {
                    // Returning out of doubly-nested while loop is not as good perf. wise
                    // on the GPU vs using "done" here.
                    done = true;
                    break;
                }
                T_maj = SampledSpectrum(1.f);
                tMin = t;

            } else {
                // Handle sample past end of majorant segment
                Float dt = seg->tMax - tMin;
                // Handle infinite _dt_ for ray majorant segment
                if (IsInf(dt))
                    dt = std::numeric_limits<Float>::max();

                T_maj *= FastExp(-dt * seg->sigma_maj);
                PBRT_DBG("Past end, added dt %f * maj[0] %f\n", dt, seg->sigma_maj[0]);
                break;
            }
        }
    }
    return SampledSpectrum(1.f);
}

}  // namespace pbrt

#endif  // PBRT_MEDIA_H
