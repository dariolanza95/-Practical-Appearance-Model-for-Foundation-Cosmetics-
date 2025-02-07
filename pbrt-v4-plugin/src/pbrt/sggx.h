// Adaption of Tizian implementation: Practical Multiple-Scattering Sheen Using Linearly Transformed Cosines (https://github.com/tizian/ltc-sheen/tree/master)
// Original source code: The SGGX Microflake Distribution (Supplemental Material) https://drive.google.com/file/d/0BzvWIdpUpRx_djVyMG9jMnltdTg/view?resourcekey=0-VTvjBPesVjrNy4SH2ShqDw

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_SGGX_H
#define PBRT_MATERIALS_SGGX_H

// materials/sggx.h*
#include "pbrt.h"
#include <stdio.h>
#include <iostream>
namespace pbrt {

/*
    From the reference implementaiton of
    "The SGGX Microflake Distribution" by Heitz et al. 2015.
*/

namespace sggx {

struct Ellipsoid {
    Float xx, yy, zz, xy, xz, yz;

    Ellipsoid() {}

    //The SGGX Microflake Distribution: Initializing the SGGX Representation (General Approach, Page 4, Equation 9)

    //The SGGX Microflake Distribution: Initializing the SGGX Representation (Surface-like Distributions, Page 5)
    static Ellipsoid fromSurface(const Vector3f &t, const Float sigma)  {
        Float x = t.x;
        Float y = t.y;
        Float z = t.z;
        Ellipsoid S;
        S.xx = x*x + sigma*sigma*(y*y + z*z);
        S.yy = y*y + sigma*sigma*(x*x + z*z);
        S.zz = z*z + sigma*sigma*(x*x + y*y);
        S.xy = x*y - sigma*sigma*x*y;
        S.xz = x*z - sigma*sigma*x*z;
        S.yz = y*z - sigma*sigma*y*z;
        return S;
    }
    
    //The SGGX Microflake Distribution: Initializing the SGGX Representation (Fiber-like Distributions, Page 5)
    static Ellipsoid fromFiber(const Vector3f &t, Float sigma) {
        Float x = t.x;
        Float y = t.y;
        Float z = t.z;
        Ellipsoid S;
        S.xx = sigma*sigma*x*x + y*y + z*z;
        S.yy = sigma*sigma*y*y + x*x + z*z;
        S.zz = sigma*sigma*z*z + x*x + y*y;
        S.xy = sigma*sigma*x*y - x*y;
        S.xz = sigma*sigma*x*z - x*z;
        S.yz = sigma*sigma*y*z - y*z;
        return S;
    }
};

//The SGGX Microflake Distribution (Supplemental Material): Page 13
inline Float D(const Vector3f &wm, const Ellipsoid &S) {
    Float detS = S.xx*S.yy*S.zz -
                 S.xx*S.yz*S.yz -
                 S.yy*S.xz*S.xz -
                 S.zz*S.xy*S.xy +
                 2.f*S.xy*S.xz*S.yz;

    Float den = wm.x*wm.x*(S.yy*S.zz-S.yz*S.yz) +
                wm.y*wm.y*(S.xx*S.zz-S.xz*S.xz) +
                wm.z*wm.z*(S.xx*S.yy-S.xy*S.xy)
              + 2.f*(wm.x*wm.y*(S.xz*S.yz-S.zz*S.xy) +
                     wm.x*wm.z*(S.xy*S.yz-S.yy*S.xz) +
                     wm.y*wm.z*(S.xy*S.xz-S.xx*S.yz));

    return std::pow(std::abs(detS), 1.5f) / (Pi*den*den);
}

//The SGGX Microflake Distribution (Supplemental Material): Page 13
inline Float projArea(const Vector3f &wi, const Ellipsoid &S) {
    Float sigma2 = wi.x*wi.x*S.xx +
                   wi.y*wi.y*S.yy +
                   wi.z*wi.z*S.zz
                 + 2.f*(wi.x*wi.y*S.xy +
                        wi.x*wi.z*S.xz +
                        wi.y*wi.z*S.yz);
    // return sigma2 > 0.0f? std::sqrt(sigma2):0.0f;
    return std::sqrt(std::max(0.f, sigma2));
}

inline Float visibleD(const Vector3f &wi, const Vector3f &wm, const Ellipsoid &S) {
    return std::max(0.f, std::min(1.f, Dot(wi, wm))) * D(wm, S) / projArea(wi, S);
}

//The SGGX Microflake Distribution (Supplemental Material): Page 14
inline void buildOrthonormalBasis(Vector3f &omega_1, Vector3f &omega_2, Vector3f omega_3)
{
    if (omega_3.z < -0.9999999f)
    {
        omega_1 = Vector3f(0.0f, -1.0f, 0.0f);
        omega_2 = Vector3f(-1.0f, 0.0f, 0.0f);
    }
    else {
        const Float a = 1.0f / (1.0f + omega_3.z);
        const Float b = -omega_3.x*omega_3.y*a;
        omega_1 = Vector3f(1.0f - omega_3.x*omega_3.x*a, b, -omega_3.x);
        omega_2 = Vector3f(b, 1.0f - omega_3.y*omega_3.y*a, -omega_3.y);
    }
}


//The SGGX Microflake Distribution (Supplemental Material): Page 14
inline Vector3f sampleVisibleD(const Vector3f &wi, const Ellipsoid &S, Float u1, Float u2) {
    // Generate sample (u, v, w)
    Float r   = std::sqrt(u1),
          phi = 2.f * Pi * u2;

    Float u = r * std::sin(phi),
          v = r * std::cos(phi),
          w = std::sqrt(1.f - u*u - v*v);

    // Build orthonormal basis
    Vector3f wk, wj;
    //CoordinateSystem(wi, &wk, &wj);
    buildOrthonormalBasis(wk, wj,wi);
    // Project S in this basis
    Float S_kk = wk.x*wk.x*S.xx + wk.y*wk.y*S.yy + wk.z*wk.z*S.zz
               + 2.f * (wk.x*wk.y*S.xy + wk.x*wk.z*S.xz + wk.y*wk.z*S.yz);
    Float S_jj = wj.x*wj.x*S.xx + wj.y*wj.y*S.yy + wj.z*wj.z*S.zz
               + 2.f * (wj.x*wj.y*S.xy + wj.x*wj.z*S.xz + wj.y*wj.z*S.yz);
    Float S_ii = wi.x*wi.x*S.xx + wi.y*wi.y*S.yy + wi.z*wi.z*S.zz
               + 2.f * (wi.x*wi.y*S.xy + wi.x*wi.z*S.xz + wi.y*wi.z*S.yz);
    Float S_kj = wk.x*wj.x*S.xx + wk.y*wj.y*S.yy + wk.z*wj.z*S.zz
               + (wk.x*wj.y + wk.y*wj.x)*S.xy
               + (wk.x*wj.z + wk.z*wj.x)*S.xz
               + (wk.y*wj.z + wk.z*wj.y)*S.yz;
    Float S_ki = wk.x*wi.x*S.xx + wk.y*wi.y*S.yy + wk.z*wi.z*S.zz
               + (wk.x*wi.y + wk.y*wi.x)*S.xy
               + (wk.x*wi.z + wk.z*wi.x)*S.xz
               + (wk.y*wi.z + wk.z*wi.y)*S.yz;
    Float S_ji = wj.x*wi.x*S.xx + wj.y*wi.y*S.yy + wj.z*wi.z*S.zz
               + (wj.x*wi.y + wj.y*wi.x)*S.xy
               + (wj.x*wi.z + wj.z*wi.x)*S.xz
               + (wj.y*wi.z + wj.z*wi.y)*S.yz;

    // Compute normal
    Float sqrtDetSkji = std::sqrt(std::abs(S_kk*S_jj*S_ii - S_kj*S_kj*S_ii - S_ki*S_ki*S_jj - S_ji*S_ji*S_kk + 2.f*S_kj*S_ki*S_ji));
    Float invSqrtSii = 1.f / std::sqrt(S_ii);
    Float tmp = std::sqrt(S_jj*S_ii - S_ji*S_ji);
    Vector3f Mk(sqrtDetSkji / tmp, 0.f, 0.f);
    Vector3f Mj(-invSqrtSii*(S_ki*S_ji - S_kj*S_ii) / tmp, invSqrtSii*tmp, 0);
    Vector3f Mi(invSqrtSii*S_ki, invSqrtSii*S_ji, invSqrtSii*S_ii);
    Vector3f wm_kji = Normalize(u*Mk + v*Mj + w*Mi);

    // Rotate back to world basis
    return wm_kji.x * wk + wm_kji.y * wj + wm_kji.z * wi;
}




//The SGGX Microflake Distribution (Supplemental Material): Page 15
inline Float evalPhaseSpecular(const Vector3f &wi, const Vector3f &wo, const Ellipsoid &S) {
    // std::cout<<"EvalSpecular"<<std::endl;
    Vector3f wh = Normalize(wi + wo);
    //if one index is nan also the others will be
    if(IsNaN(wh[0])){
        // std::cout<<"skipping"<<std::endl;
        return 0.f;
    }
    // SameHemisphere(wi, wo)
    //Vector3f wh = Normalize(wo);
    Float s = projArea(wi, S);
    Float rez = 0.25f * D(wh, S) / s;

    // std::cout<<"eval "<< s <<  " rez "<< rez << " wh " << wh << " wi "<< wi << " wo " << wo << std::endl;
    return s == 0.f ? 0.f : 0.25f * D(wh, S) / s;
}

//The SGGX Microflake Distribution (Supplemental Material): Page 15
inline Vector3f samplePhaseSpecular(const Vector3f &wi, const Ellipsoid &S, Float u1, Float u2) {
    // std::cout<<"SampleSpecular"<<std::endl;

    Vector3f wm = sampleVisibleD(wi, S, u1, u2);
    return Normalize(2.f*Dot(wm, wi)*wm - wi);
}


//The SGGX Microflake Distribution (Supplemental Material): Page 16
inline Float evalPhaseDiffuse(const Vector3f &wi, const Vector3f &wo, const Ellipsoid &S, Float u1, Float u2)
{
    // std::cout<<"EvalDiff"<<std::endl;

    Vector3f wm = sampleVisibleD(wi, S, u1, u2);// sample VNDF

    return 1.0f / Pi * std::max(0.0f, Dot(wo, wm));// eval diffuse
}



//The SGGX Microflake Distribution (Supplemental Material): Page 16
inline Vector3f samplePhaseDiffuse(const Vector3f &wi, const Ellipsoid &S, Float u1, Float u2, Float u3, Float u4)
{
    // std::cout<<"SampleDiff"<<std::endl;

    Vector3f wi_t(wi);
    Vector3f wm = sampleVisibleD(wi_t, S, u1, u2);// sample VNDF

    // sample diffuse reflection
    Vector3f w1, w2;
    buildOrthonormalBasis(w1, w2, wm);
    Float r1 = 2.0f*u3 - 1.0f;
    Float r2 = 2.0f*u4 - 1.0f;

    // concentric map code from
    // http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
    Float phi, r;

    if (r1 == 0 && r2 == 0) {
        r = phi = 0;
    } else if (r1*r1 > r2*r2) {
        r = r1;
        phi = (Pi/4.0f) * (r2/r1);
    } else {
        r = r2;
        phi = (Pi/2.0f) - (r1/r2) * (Pi/4.0f);
    }

    Float x = r*cosf(phi);
    Float y = r*sinf(phi);
    Float z = sqrtf(1.0f - x*x - y*y);
    Vector3f wo = x*w1 + y*w2 + z*wm;

    return wo;
}
inline Float evalPhaseType(const Float diff_particles_perc, const Vector3f &wi, const Vector3f &wo, const Ellipsoid &S,Float u1, Float u2,Float u3,Float g = 0.0f){
    if (u3 < diff_particles_perc){
        return HenyeyGreenstein(Dot(wo, wi), g); 
        // return evalPhaseDiffuse(wi,wo,S,u1,u2);
    }else{        
        return evalPhaseSpecular(wi,wo,S);

    }
}
inline Vector3f samplePhaseType(const Float diff_particles_perc,const Vector3f &wi, const Ellipsoid &S, Float u1, Float u2, Float u3, Float u4,Float u5){
    if (u5 < diff_particles_perc){
        return samplePhaseDiffuse(wi, S, u1,u2, u3, u4);
        // return evalPhaseDiffuse(wi,wo,S,u1,u2);


    }else{
        return samplePhaseSpecular(wi, S, u1, u2);

    }
}

}  // namespace sggx



}  // namespace pbrt

#endif  // PBRT_MATERIALS_SGGX_H
