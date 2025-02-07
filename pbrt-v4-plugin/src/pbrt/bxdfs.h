// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_BXDFS_H
#define PBRT_BXDFS_H

#include <pbrt/pbrt.h>

#include <pbrt/base/bxdf.h>
#include <pbrt/interaction.h>
#include <pbrt/media.h>
#include <pbrt/options.h>
#include <pbrt/util/math.h>
#include <pbrt/util/memory.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/scattering.h>
#include <pbrt/util/spectrum.h>
#include <pbrt/util/taggedptr.h>
#include <pbrt/util/vecmath.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>



namespace pbrt {
struct CosmeticStructBxdf{
    
    CosmeticStructBxdf() = default;
          
    //  CosmeticStructBxdf() {};
    PBRT_CPU_GPU
     CosmeticStructBxdf(Float sigma_t,Float thickness,Float diff_particles_perc,Float platelets_roughness,Float std_size, Float platelets_orientation, Float std_or, Float g1, Float g2,Float w_g,const SampledSpectrum &albedo_flakes,const SampledSpectrum &albedo_diff,bool use_specular):     
          sigma_t(sigma_t),
          thickness(thickness),
          diff_conc(diff_particles_perc),
          platelets_orientation(platelets_orientation),
          platelets_roughness(platelets_roughness),
          std_or(std_or),
          std_size(std_size),
          g1(g1),
          g2(g2),
          w_g(w_g),
          albedo_flakes(albedo_flakes),
          albedo_diff(albedo_diff),
          use_specular(use_specular) {
            // std::cout<<"platelets_orientation " << platelets_orientation << "sigma_t "<< sigma_t<< std::endl;
            // std::cout<<"platelets_orientation " << platelets_orientation << "sigma_t "<< sigma_t<< std::endl;
          }
      Float sigma_t,diff_conc, thickness; 
      Float platelets_roughness,std_size,platelets_orientation,std_or;
      Float g1, g2, w_g;
      bool use_specular;
      SampledSpectrum albedo_flakes, albedo_diff;

};



// DiffuseBxDF Definition


class DiffuseBxDF {
  public:
    // DiffuseBxDF Public Methods
    DiffuseBxDF() = default;
    PBRT_CPU_GPU
    DiffuseBxDF(SampledSpectrum R) : R(R) {}

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        if (!SameHemisphere(wo, wi))
            return SampledSpectrum(0.f);
        return R * InvPi;
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(
        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            return {};
        // Sample cosine-weighted hemisphere to compute _wi_ and _pdf_
        Vector3f wi = SampleCosineHemisphere(u);
        if (wo.z < 0)
            wi.z *= -1;
        Float pdf = CosineHemispherePDF(AbsCosTheta(wi));

        return BSDFSample(R * InvPi, wi, pdf, BxDFFlags::DiffuseReflection);
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        if (!(sampleFlags & BxDFReflTransFlags::Reflection) || !SameHemisphere(wo, wi))
            return 0;
        return CosineHemispherePDF(AbsCosTheta(wi));
    }

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "DiffuseBxDF"; }

    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() {}

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        return R ? BxDFFlags::DiffuseReflection : BxDFFlags::Unset;
    }

  private:
    SampledSpectrum R;
};

// DiffuseTransmissionBxDF Definition
class DiffuseTransmissionBxDF {
  public:
    // DiffuseTransmissionBxDF Public Methods
    DiffuseTransmissionBxDF() = default;
    PBRT_CPU_GPU
    DiffuseTransmissionBxDF(SampledSpectrum R, SampledSpectrum T) : R(R), T(T) {}

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        return SameHemisphere(wo, wi) ? (R * InvPi) : (T * InvPi);
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(
        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        // Compute reflection and transmission probabilities for diffuse BSDF
        Float pr = R.MaxComponentValue(), pt = T.MaxComponentValue();
        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            pr = 0;
        if (!(sampleFlags & BxDFReflTransFlags::Transmission))
            pt = 0;
        if (pr == 0 && pt == 0)
            return {};

        // Randomly sample diffuse BSDF reflection or transmission
        if (uc < pr / (pr + pt)) {
            // Sample diffuse BSDF reflection
            Vector3f wi = SampleCosineHemisphere(u);
            if (wo.z < 0)
                wi.z *= -1;
            Float pdf = CosineHemispherePDF(AbsCosTheta(wi)) * pr / (pr + pt);
            return BSDFSample(f(wo, wi, mode), wi, pdf, BxDFFlags::DiffuseReflection);

        } else {
            // Sample diffuse BSDF transmission
            Vector3f wi = SampleCosineHemisphere(u);
            if (wo.z > 0)
                wi.z *= -1;
            Float pdf = CosineHemispherePDF(AbsCosTheta(wi)) * pt / (pr + pt);
            return BSDFSample(f(wo, wi, mode), wi, pdf, BxDFFlags::DiffuseTransmission);
        }
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        // Compute reflection and transmission probabilities for diffuse BSDF
        Float pr = R.MaxComponentValue(), pt = T.MaxComponentValue();
        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            pr = 0;
        if (!(sampleFlags & BxDFReflTransFlags::Transmission))
            pt = 0;
        if (pr == 0 && pt == 0)
            return {};

        if (SameHemisphere(wo, wi))
            return pr / (pr + pt) * CosineHemispherePDF(AbsCosTheta(wi));
        else
            return pt / (pr + pt) * CosineHemispherePDF(AbsCosTheta(wi));
    }

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "DiffuseTransmissionBxDF"; }

    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() {}

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        return ((R ? BxDFFlags::DiffuseReflection : BxDFFlags::Unset) |
                (T ? BxDFFlags::DiffuseTransmission : BxDFFlags::Unset));
    }

  private:
    // DiffuseTransmissionBxDF Private Members
    SampledSpectrum R, T;
};
// NormalizedFresnelBxDF Definition
class NormalizedFresnelBxDF {
  public:
    // NormalizedFresnelBxDF Public Methods
    NormalizedFresnelBxDF() = default;
    PBRT_CPU_GPU
    NormalizedFresnelBxDF(Float eta) : eta(eta) {}

    PBRT_CPU_GPU
    BSDFSample Sample_f(Vector3f wo, Float uc, Point2f u, TransportMode mode,
                        BxDFReflTransFlags sampleFlags) const {
        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            return {};

        // Cosine-sample the hemisphere, flipping the direction if necessary
        Vector3f wi = SampleCosineHemisphere(u);
        if (wo.z < 0)
            wi.z *= -1;
        return BSDFSample(f(wo, wi, mode), wi, PDF(wo, wi, mode, sampleFlags),
                          BxDFFlags::DiffuseReflection);
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags) const {
        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            return 0;
        return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
    }

    PBRT_CPU_GPU
    void Regularize() {}

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "NormalizedFresnelBxDF"; }

    std::string ToString() const;

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        return BxDFFlags(BxDFFlags::Reflection | BxDFFlags::Diffuse);
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        if (!SameHemisphere(wo, wi))
            return SampledSpectrum(0.f);
        // Compute $\Sw$ factor for BSSRDF value
        Float c = 1 - 2 * FresnelMoment1(1 / eta);
        SampledSpectrum f((1 - FrDielectric(CosTheta(wi), eta)) / (c * Pi));

        // Update BSSRDF transmission term to account for adjoint light transport
        if (mode == TransportMode::Radiance)
            f *= Sqr(eta);

        return f;
    }

  private:
    Float eta;
};


// DielectricBxDF Definition
class DielectricBxDF {
  public:
    // DielectricBxDF Public Methods
    DielectricBxDF() = default;
    PBRT_CPU_GPU
    DielectricBxDF(Float eta, TrowbridgeReitzDistribution mfDistrib)
        : eta(eta), mfDistrib(mfDistrib) {}

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        BxDFFlags flags = (eta == 1) ? BxDFFlags::Transmission
                                     : (BxDFFlags::Reflection | BxDFFlags::Transmission);
        return flags |
               (mfDistrib.EffectivelySmooth() ? BxDFFlags::Specular : BxDFFlags::Glossy);
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(
        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const;

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const;
    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const;

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "DielectricBxDF"; }

    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() { mfDistrib.Regularize(); }

  private:
    // DielectricBxDF Private Members
    Float eta;
    TrowbridgeReitzDistribution mfDistrib;
};

// ThinDielectricBxDF Definition
class ThinDielectricBxDF {
  public:
    // ThinDielectricBxDF Public Methods
    ThinDielectricBxDF() = default;
    PBRT_CPU_GPU
    ThinDielectricBxDF(Float eta) : eta(eta) {}

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        return SampledSpectrum(0);
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(Vector3f wo, Float uc, Point2f u,
                                        TransportMode mode,
                                        BxDFReflTransFlags sampleFlags) const {
        Float R = FrDielectric(AbsCosTheta(wo), eta), T = 1 - R;
        // Compute _R_ and _T_ accounting for scattering between interfaces
        if (R < 1) {
            R += Sqr(T) * R / (1 - Sqr(R));
            T = 1 - R;
        }

        // Compute probabilities _pr_ and _pt_ for sampling reflection and transmission
        Float pr = R, pt = T;
        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            pr = 0;
        if (!(sampleFlags & BxDFReflTransFlags::Transmission))
            pt = 0;
        if (pr == 0 && pt == 0)
            return {};

        if (uc < pr / (pr + pt)) {
            // Sample perfect specular dielectric BRDF
            Vector3f wi(-wo.x, -wo.y, wo.z);
            SampledSpectrum fr(R / AbsCosTheta(wi));
            return BSDFSample(fr, wi, pr / (pr + pt), BxDFFlags::SpecularReflection);

        } else {
            // Sample perfect specular transmission at thin dielectric interface
            Vector3f wi = -wo;
            SampledSpectrum ft(T / AbsCosTheta(wi));
            return BSDFSample(ft, wi, pt / (pr + pt), BxDFFlags::SpecularTransmission);
        }
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags) const {
        return 0;
    }

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "ThinDielectricBxDF"; }

    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() { /* TODO */
    }

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        return (BxDFFlags::Reflection | BxDFFlags::Transmission | BxDFFlags::Specular);
    }

  private:
    Float eta;
};

// ConductorBxDF Definition
class ConductorBxDF {
  public:
    // ConductorBxDF Public Methods
    ConductorBxDF() = default;
    PBRT_CPU_GPU
    ConductorBxDF(const TrowbridgeReitzDistribution &mfDistrib, SampledSpectrum eta,
                  SampledSpectrum k)
        : mfDistrib(mfDistrib), eta(eta), k(k) {}

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        return mfDistrib.EffectivelySmooth() ? BxDFFlags::SpecularReflection
                                             : BxDFFlags::GlossyReflection;
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(
        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            return {};
        if (mfDistrib.EffectivelySmooth()) {
            // Sample perfect specular conductor BRDF
            Vector3f wi(-wo.x, -wo.y, wo.z);
            SampledSpectrum f = FrComplex(AbsCosTheta(wi), eta, k) / AbsCosTheta(wi);
            return BSDFSample(f, wi, 1, BxDFFlags::SpecularReflection);
        }
        // Sample rough conductor BRDF
        // Sample microfacet normal $\wm$ and reflected direction $\wi$
        if (wo.z == 0)
            return {};
        Vector3f wm = mfDistrib.Sample_wm(wo, u);
        Vector3f wi = Reflect(wo, wm);
        if (!SameHemisphere(wo, wi))
            return {};

        // Compute PDF of _wi_ for microfacet reflection
        Float pdf = mfDistrib.PDF(wo, wm) / (4 * AbsDot(wo, wm));

        Float cosTheta_o = AbsCosTheta(wo), cosTheta_i = AbsCosTheta(wi);
        if (cosTheta_i == 0 || cosTheta_o == 0)
            return {};
        // Evaluate Fresnel factor _F_ for conductor BRDF
        SampledSpectrum F = FrComplex(AbsDot(wo, wm), eta, k);

        SampledSpectrum f =
            mfDistrib.D(wm) * F * mfDistrib.G(wo, wi) / (4 * cosTheta_i * cosTheta_o);
        return BSDFSample(f, wi, pdf, BxDFFlags::GlossyReflection);
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        if (!SameHemisphere(wo, wi))
            return {};
        if (mfDistrib.EffectivelySmooth())
            return {};
        // Evaluate rough conductor BRDF
        // Compute cosines and $\wm$ for conductor BRDF
        Float cosTheta_o = AbsCosTheta(wo), cosTheta_i = AbsCosTheta(wi);
        if (cosTheta_i == 0 || cosTheta_o == 0)
            return {};
        Vector3f wm = wi + wo;
        if (LengthSquared(wm) == 0)
            return {};
        wm = Normalize(wm);

        // Evaluate Fresnel factor _F_ for conductor BRDF
        SampledSpectrum F = FrComplex(AbsDot(wo, wm), eta, k);

        return mfDistrib.D(wm) * F * mfDistrib.G(wo, wi) / (4 * cosTheta_i * cosTheta_o);
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags) const {
        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            return 0;
        if (!SameHemisphere(wo, wi))
            return 0;
        if (mfDistrib.EffectivelySmooth())
            return 0;
        // Evaluate sampling PDF of rough conductor BRDF
        Vector3f wm = wo + wi;
        CHECK_RARE(1e-5f, LengthSquared(wm) == 0);
        if (LengthSquared(wm) == 0)
            return 0;
        wm = FaceForward(Normalize(wm), Normal3f(0, 0, 1));
        return mfDistrib.PDF(wo, wm) / (4 * AbsDot(wo, wm));
    }

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "ConductorBxDF"; }
    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() { mfDistrib.Regularize(); }

  private:
    // ConductorBxDF Private Members
    TrowbridgeReitzDistribution mfDistrib;
    SampledSpectrum eta, k;
};

// TopOrBottomBxDF Definition
template <typename TopBxDF, typename BottomBxDF>
class TopOrBottomBxDF {
  public:
    // TopOrBottomBxDF Public Methods
    TopOrBottomBxDF() = default;
    PBRT_CPU_GPU
    TopOrBottomBxDF &operator=(const TopBxDF *t) {
        top = t;
        bottom = nullptr;
        return *this;
    }
    PBRT_CPU_GPU
    TopOrBottomBxDF &operator=(const BottomBxDF *b) {
        bottom = b;
        top = nullptr;
        return *this;
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        return top ? top->f(wo, wi, mode) : bottom->f(wo, wi, mode);
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(
        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        return top ? top->Sample_f(wo, uc, u, mode, sampleFlags)
                   : bottom->Sample_f(wo, uc, u, mode, sampleFlags);
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        return top ? top->PDF(wo, wi, mode, sampleFlags)
                   : bottom->PDF(wo, wi, mode, sampleFlags);
    }

    PBRT_CPU_GPU
    BxDFFlags Flags() const { return top ? top->Flags() : bottom->Flags(); }

  private:
    const TopBxDF *top = nullptr;
    const BottomBxDF *bottom = nullptr;
};

// LayeredBxDF Definition
template <typename TopBxDF, typename BottomBxDF, bool twoSided>
class LayeredBxDF {
  public:
    // LayeredBxDF Public Methods
    LayeredBxDF() = default;
    PBRT_CPU_GPU
    LayeredBxDF(TopBxDF top, BottomBxDF bottom, Float thickness,
                const SampledSpectrum &albedo, Float g, int maxDepth, int nSamples)
        : top(top),
          bottom(bottom),
          thickness(std::max(thickness, std::numeric_limits<Float>::min())),
          g(g),
          albedo(albedo),
          maxDepth(maxDepth),
          nSamples(nSamples) {}

    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() {
        top.Regularize();
        bottom.Regularize();
    }

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        BxDFFlags topFlags = top.Flags(), bottomFlags = bottom.Flags();
        CHECK(IsTransmissive(topFlags) ||
              IsTransmissive(bottomFlags));  // otherwise, why bother?

        BxDFFlags flags = BxDFFlags::Reflection;
        if (IsSpecular(topFlags))
            flags = flags | BxDFFlags::Specular;

        if (IsDiffuse(topFlags) || IsDiffuse(bottomFlags) || albedo)
            flags = flags | BxDFFlags::Diffuse;
        else if (IsGlossy(topFlags) || IsGlossy(bottomFlags))
            flags = flags | BxDFFlags::Glossy;

        if (IsTransmissive(topFlags) && IsTransmissive(bottomFlags))
            flags = flags | BxDFFlags::Transmission;

        return flags;
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        SampledSpectrum f(0.);
        // Estimate _LayeredBxDF_ value _f_ using random sampling
        // Set _wo_ and _wi_ for layered BSDF evaluation
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }

        // Determine entrance interface for layered BSDF
        TopOrBottomBxDF<TopBxDF, BottomBxDF> enterInterface;
        bool enteredTop = twoSided || wo.z > 0;
        if (enteredTop)
            enterInterface = &top;
        else
            enterInterface = &bottom;

        // Determine exit interface and exit $z$ for layered BSDF
        TopOrBottomBxDF<TopBxDF, BottomBxDF> exitInterface, nonExitInterface;
        if (SameHemisphere(wo, wi) ^ enteredTop) {
            exitInterface = &bottom;
            nonExitInterface = &top;
        } else {
            exitInterface = &top;
            nonExitInterface = &bottom;
        }
        Float exitZ = (SameHemisphere(wo, wi) ^ enteredTop) ? 0 : thickness;

        // Account for reflection at the entrance interface
        if (SameHemisphere(wo, wi))
            f = nSamples * enterInterface.f(wo, wi, mode);

        // Declare _RNG_ for layered BSDF evaluation
   
        // std::cout<< GetOptions().seed<<std::endl;
        // SetOption().seed = 42;
        // std::cout<< GetOptions().seed<<std::endl;

        RNG rng(Hash(GetOptions().seed, wo), Hash(wi));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };
        for (int s = 0; s < nSamples; ++s) {
            // Sample random walk through layers to estimate BSDF value
            // Sample transmission direction through entrance interface
            Float uc = r();
            pstd::optional<BSDFSample> wos = enterInterface.Sample_f(
                wo, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Transmission);
            if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0)
                continue;

            // Sample BSDF for virtual light from _wi_
            uc = r();
            pstd::optional<BSDFSample> wis = exitInterface.Sample_f(
                wi, uc, Point2f(r(), r()), !mode, BxDFReflTransFlags::Transmission);
            if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0)
                continue;
            // Declare state for random walk through BSDF layers
            SampledSpectrum beta = wos->f * AbsCosTheta(wos->wi) / wos->pdf;
            Float z = enteredTop ? thickness : 0;
            Vector3f w = wos->wi;
            HGPhaseFunction phase(g);  // Original code: Henyey-Greenstein phase function
            // phase.SetMicroflakes(30.0f,1.0f);//Fiber-like distributions

            // Diffuse phase function
            Point2f u1{r(), r()}, u2{r(), r()};
            // bool is_specular_phase = true;

            for (int depth = 0; depth < maxDepth; ++depth) {
                // Sample next event for layered BSDF evaluation random walk
                PBRT_DBG("beta: %f %f %f %f, w: %f %f %f, f: %f %f %f %f\n", beta[0],
                         beta[1], beta[2], beta[3], w.x, w.y, w.z, f[0], f[1], f[2],
                         f[3]);
                // Possibly terminate layered BSDF random walk with Russian roulette
                if (depth > 3 && beta.MaxComponentValue() < 0.25f) {
                    Float q = std::max<Float>(0, 1 - beta.MaxComponentValue());
                    if (r() < q)
                        break;
                    beta /= 1 - q;
                    PBRT_DBG("After RR with q = %f, beta: %f %f %f %f\n", q, beta[0],
                             beta[1], beta[2], beta[3]);
                }

                // Account for media between layers and possibly scatter
                if (!albedo) {
                    // Advance to next layer boundary and update _beta_ for transmittance
                    z = (z == thickness) ? 0 : thickness;
                    beta *= Tr(thickness, w);

                } else {
                    // Sample medium scattering for layered BSDF evaluation
                    Float sigma_t = 1;
                    Float dz = SampleExponential(r(), sigma_t / std::abs(w.z));
                    Float zp = w.z > 0 ? (z + dz) : (z - dz);
                    DCHECK_RARE(1e-5, z == zp);
                    if (z == zp)
                        continue;
                    if (0 < zp && zp < thickness) {
                        // Handle scattering event in layered BSDF medium
                        // Account for scattering through _exitInterface_ using _wis_
                        Float wt = 1;
                        if (!IsSpecular(exitInterface.Flags()))
                            wt = PowerHeuristic(1, wis->pdf, 1, phase.PDF(-w, -wis->wi));
                        f += beta * albedo * phase.p(-w, -wis->wi) * wt *
                             Tr(zp - exitZ, wis->wi) * wis->f / wis->pdf;

                        // Sample phase function and update layered path state
                        Point2f u{r(), r()};
                        pstd::optional<PhaseFunctionSample> ps = phase.Sample_p(-w, Point2f(r(), r()) );
                        if (!ps || ps->pdf == 0 || ps->wi.z == 0)
                            continue;
                        beta *= albedo * ps->p / ps->pdf;
                        w = ps->wi;
                        z = zp;

                        // Possibly account for scattering through _exitInterface_
                        if (((z < exitZ && w.z > 0) || (z > exitZ && w.z < 0)) &&
                            !IsSpecular(exitInterface.Flags())) {
                            // Account for scattering through _exitInterface_
                            SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                            if (fExit) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                Float wt = PowerHeuristic(1, ps->pdf, 1, exitPDF);
                                f += beta * Tr(zp - exitZ, ps->wi) * fExit * wt;
                            }
                        }

                        continue;
                    }
                    z = Clamp(zp, 0, thickness);
                }

                // Account for scattering at appropriate interface
                if (z == exitZ) {
                    // Account for reflection at _exitInterface_
                    Float uc = r();
                    pstd::optional<BSDFSample> bs = exitInterface.Sample_f(
                        -w, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        break;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;
                    w = bs->wi;

                } else {
                    // Account for scattering at _nonExitInterface_
                    if (!IsSpecular(nonExitInterface.Flags())) {
                        // Add NEE contribution along presampled _wis_ direction
                        Float wt = 1;
                        if (!IsSpecular(exitInterface.Flags()))
                            wt = PowerHeuristic(1, wis->pdf, 1,
                                                nonExitInterface.PDF(-w, -wis->wi, mode));
                        f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                             AbsCosTheta(wis->wi) * wt * Tr(thickness, wis->wi) * wis->f /
                             wis->pdf;
                    }
                    // Sample new direction using BSDF at _nonExitInterface_
                    Float uc = r();
                    Point2f u(r(), r());
                    pstd::optional<BSDFSample> bs = nonExitInterface.Sample_f(
                        -w, uc, u, mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        break;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;
                    w = bs->wi;

                    if (!IsSpecular(exitInterface.Flags())) {
                        // Add NEE contribution along direction from BSDF sample
                        SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                        if (fExit) {
                            Float wt = 1;
                            if (!IsSpecular(nonExitInterface.Flags())) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                wt = PowerHeuristic(1, bs->pdf, 1, exitPDF);
                            }
                            f += beta * Tr(thickness, bs->wi) * fExit * wt;
                        }
                    }
                }
            }
        }

        return f / nSamples;
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(
        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // Set _wo_ for layered BSDF sampling
        bool flipWi = false;
        if (twoSided && wo.z < 0) {
            wo = -wo;
            flipWi = true;
        }

        // Sample BSDF at entrance interface to get initial direction _w_
        bool enteredTop = twoSided || wo.z > 0;
        pstd::optional<BSDFSample> bs =
            enteredTop ? top.Sample_f(wo, uc, u, mode) : bottom.Sample_f(wo, uc, u, mode);
        if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
            return {};
        if (bs->IsReflection()) {
            if (flipWi)
                bs->wi = -bs->wi;
            bs->pdfIsProportional = true;
            return bs;
        }
        Vector3f w = bs->wi;
        bool specularPath = bs->IsSpecular();

        // Declare _RNG_ for layered BSDF sampling
        RNG rng(Hash(GetOptions().seed, wo), Hash(uc, u));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Declare common variables for layered BSDF sampling
        SampledSpectrum f = bs->f * AbsCosTheta(bs->wi);
        Float pdf = bs->pdf;
        Float z = enteredTop ? thickness : 0;
        HGPhaseFunction phase(g);  // Original code: Henyey-Greenstein phase function

        // Diffuse phase function
        // Point2f u1{r(), r()}, u2{r(), r()};

        for (int depth = 0; depth < maxDepth; ++depth) {
            // Follow random walk through layers to sample layered BSDF
            // Possibly terminate layered BSDF sampling with Russian Roulette
            Float rrBeta = f.MaxComponentValue() / pdf;
            if (depth > 3 && rrBeta < 0.25f) {
                Float q = std::max<Float>(0, 1 - rrBeta);
                if (r() < q)
                    return {};
                pdf *= 1 - q;
            }
            if (w.z == 0)
                return {};

            if (albedo) {
                // Sample potential scattering event in layered medium
                Float sigma_t = 1;
                Float dz = SampleExponential(r(), sigma_t / AbsCosTheta(w));
                Float zp = w.z > 0 ? (z + dz) : (z - dz);
                CHECK_RARE(1e-5, zp == z);
                if (zp == z)
                    return {};
                if (0 < zp && zp < thickness) {
                    // Update path state for valid scattering event between interfaces

                    pstd::optional<PhaseFunctionSample> ps =
                        phase.Sample_p(-w,  Point2f(r(), r()) );
                    if (!ps || ps->pdf == 0 || ps->wi.z == 0)
                        return {};
                    f *= albedo * ps->p;
                    pdf *= ps->pdf;
                    specularPath = false;
                    w = ps->wi;
                    z = zp;

                    continue;
                }
                z = Clamp(zp, 0, thickness);
                if (z == 0)
                    DCHECK_LT(w.z, 0);
                else
                    DCHECK_GT(w.z, 0);

            } else {
                // Advance to the other layer interface
                z = (z == thickness) ? 0 : thickness;
                f *= Tr(thickness, w);
            }
            // Initialize _interface_ for current interface surface
#ifdef interface  // That's enough out of you, Windows.
#undef interface
#endif
            TopOrBottomBxDF<TopBxDF, BottomBxDF> interface;
            if (z == 0)
                interface = &bottom;
            else
                interface = &top;

            // Sample interface BSDF to determine new path direction
            Float uc = r();
            Point2f u(r(), r());
            pstd::optional<BSDFSample> bs = interface.Sample_f(-w, uc, u, mode);
            if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                return {};
            f *= bs->f;
            pdf *= bs->pdf;
            specularPath &= bs->IsSpecular();
            w = bs->wi;

            // Return _BSDFSample_ if path has left the layers
            if (bs->IsTransmission()) {
                BxDFFlags flags = SameHemisphere(wo, w) ? BxDFFlags::Reflection
                                                        : BxDFFlags::Transmission;
                flags |= specularPath ? BxDFFlags::Specular : BxDFFlags::Glossy;
                if (flipWi)
                    w = -w;
                return BSDFSample(f, w, pdf, flags, 1.f, true);
            }

            // Scale _f_ by cosine term after scattering at the interface
            f *= AbsCosTheta(bs->wi);
        }
        return {};
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // Set _wo_ and _wi_ for layered BSDF evaluation
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }

        // Declare _RNG_ for layered PDF evaluation
        RNG rng(Hash(GetOptions().seed, wi), Hash(wo));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Update _pdfSum_ for reflection at the entrance layer
        bool enteredTop = twoSided || wo.z > 0;
        Float pdfSum = 0;
        if (SameHemisphere(wo, wi)) {
            auto reflFlag = BxDFReflTransFlags::Reflection;
            pdfSum += enteredTop ? nSamples * top.PDF(wo, wi, mode, reflFlag)
                                 : nSamples * bottom.PDF(wo, wi, mode, reflFlag);
        }

        for (int s = 0; s < nSamples; ++s) {
            // Evaluate layered BSDF PDF sample
            if (SameHemisphere(wo, wi)) {
                // Evaluate TRT term for PDF estimate
                TopOrBottomBxDF<TopBxDF, BottomBxDF> rInterface, tInterface;
                if (enteredTop) {
                    rInterface = &bottom;
                    tInterface = &top;
                } else {
                    rInterface = &top;
                    tInterface = &bottom;
                }
                // Sample _tInterface_ to get direction into the layers
                auto trans = BxDFReflTransFlags::Transmission;
                pstd::optional<BSDFSample> wos, wis;
                wos = tInterface.Sample_f(wo, r(), {r(), r()}, mode, trans);
                wis = tInterface.Sample_f(wi, r(), {r(), r()}, !mode, trans);

                // Update _pdfSum_ accounting for TRT scattering events
                if (wos && wos->f && wos->pdf > 0 && wis && wis->f && wis->pdf > 0) {
                    if (!IsNonSpecular(tInterface.Flags()))
                        pdfSum += rInterface.PDF(-wos->wi, -wis->wi, mode);
                    else {
                        // Use multiple importance sampling to estimate PDF product
                        pstd::optional<BSDFSample> rs =
                            rInterface.Sample_f(-wos->wi, r(), {r(), r()}, mode);
                        if (rs && rs->f && rs->pdf > 0) {
                            if (!IsNonSpecular(rInterface.Flags()))
                                pdfSum += tInterface.PDF(-rs->wi, wi, mode);
                            else {
                                // Compute MIS-weighted estimate of Equation
                                // (\ref{eq:pdf-triple-canceled-one})
                                Float rPDF = rInterface.PDF(-wos->wi, -wis->wi, mode);
                                Float wt = PowerHeuristic(1, wis->pdf, 1, rPDF);
                                pdfSum += wt * rPDF;

                                Float tPDF = tInterface.PDF(-rs->wi, wi, mode);
                                wt = PowerHeuristic(1, rs->pdf, 1, tPDF);
                                pdfSum += wt * tPDF;
                            }
                        }
                    }
                }

            } else {
                // Evaluate TT term for PDF estimate
                TopOrBottomBxDF<TopBxDF, BottomBxDF> toInterface, tiInterface;
                if (enteredTop) {
                    toInterface = &top;
                    tiInterface = &bottom;
                } else {
                    toInterface = &bottom;
                    tiInterface = &top;
                }

                Float uc = r();
                Point2f u(r(), r());
                pstd::optional<BSDFSample> wos = toInterface.Sample_f(wo, uc, u, mode);
                if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0 ||
                    wos->IsReflection())
                    continue;

                uc = r();
                u = Point2f(r(), r());
                pstd::optional<BSDFSample> wis = tiInterface.Sample_f(wi, uc, u, !mode);
                if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0 ||
                    wis->IsReflection())
                    continue;

                if (IsSpecular(toInterface.Flags()))
                    pdfSum += tiInterface.PDF(-wos->wi, wi, mode);
                else if (IsSpecular(tiInterface.Flags()))
                    pdfSum += toInterface.PDF(wo, -wis->wi, mode);
                else
                    pdfSum += (toInterface.PDF(wo, -wis->wi, mode) +
                               tiInterface.PDF(-wos->wi, wi, mode)) /
                              2;
            }
        }
        // Return mixture of PDF estimate and constant PDF
        return Lerp(0.9f, 1 / (4 * Pi), pdfSum / nSamples);
    }

  private:
    // LayeredBxDF Private Methods
    PBRT_CPU_GPU
    static Float Tr(Float dz, Vector3f w) { return FastExp(-std::abs(dz / w.z)); }

    // LayeredBxDF Private Members
    TopBxDF top;
    BottomBxDF bottom;
    Float thickness, g;
    SampledSpectrum albedo;
    int maxDepth, nSamples;
};

template <typename TopBxDF, typename BottomBxDF, bool twoSided>
class MyLayeredBxDF {
  public:
    // MyLayeredBxDF Public Methods
    MyLayeredBxDF() = default;
    PBRT_CPU_GPU
    MyLayeredBxDF(TopBxDF top, BottomBxDF bottom, Float thickness, Vector3f dir,
                  Float sigma, Float mean, Float std_or, Float std_size, bool use_refl,
                  const SampledSpectrum &albedo, Float g, int maxDepth, int nSamples)
        : top(top),
          bottom(bottom),
          thickness(std::max(thickness, std::numeric_limits<Float>::min())),
          dir(dir),
          sigma(sigma),
          mean(mean),
          std_or(std_or),
          std_size(std_size),
          use_specular(use_refl),
          g(g),
          albedo(albedo),
          maxDepth(maxDepth),
          nSamples(nSamples) {
        // std::cout<<"Mean is " << mean << std::endl;
        // std::cout<<"std SIZE is " << std_size << std::endl;
    }

    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() {
        top.Regularize();
        bottom.Regularize();
    }

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        BxDFFlags topFlags = top.Flags(), bottomFlags = bottom.Flags();
        CHECK(IsTransmissive(topFlags) ||
              IsTransmissive(bottomFlags));  // otherwise, why bother?

        BxDFFlags flags = BxDFFlags::Reflection;
        if (IsSpecular(topFlags))
            flags = flags | BxDFFlags::Specular;

        if (IsDiffuse(topFlags) || IsDiffuse(bottomFlags) || albedo)
            flags = flags | BxDFFlags::Diffuse;
        else if (IsGlossy(topFlags) || IsGlossy(bottomFlags))
            flags = flags | BxDFFlags::Glossy;

        if (IsTransmissive(topFlags) && IsTransmissive(bottomFlags))
            flags = flags | BxDFFlags::Transmission;

        return flags;
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        SampledSpectrum f(0.);
        // Next line is for debugging purposes only, D.L.
        //  wo = Vector3f(0,0,1);
        //  wi = Vector3f(0,0,1);
        //  Estimate _MyLayeredBxDF_ value _f_ using random sampling
        //  Set _wo_ and _wi_ for layered BSDF evaluation
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }

        // Determine entrance interface for layered BSDF
        TopOrBottomBxDF<TopBxDF, BottomBxDF> enterInterface;
        bool enteredTop = twoSided || wo.z > 0;
        if (enteredTop)
            enterInterface = &top;
        else
            enterInterface = &bottom;

        // Determine exit interface and exit $z$ for layered BSDF
        TopOrBottomBxDF<TopBxDF, BottomBxDF> exitInterface, nonExitInterface;
        if (SameHemisphere(wo, wi) ^ enteredTop) {
            exitInterface = &bottom;
            nonExitInterface = &top;
        } else {
            exitInterface = &top;
            nonExitInterface = &bottom;
        }
        Float exitZ = (SameHemisphere(wo, wi) ^ enteredTop) ? 0 : thickness;

        // Account for reflection at the entrance interface
        if (SameHemisphere(wo, wi))
            f = nSamples * enterInterface.f(wo, wi, mode);

        // Declare _RNG_ for layered BSDF evaluation
        RNG rng(Hash(GetOptions().seed, wo), Hash(wi));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        for (int s = 0; s < nSamples; ++s) {
            // Sample random walk through layers to estimate BSDF value
            // Sample transmission direction through entrance interface
            Float uc = r();
            pstd::optional<BSDFSample> wos = enterInterface.Sample_f(
                wo, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Transmission);
            if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0)
                continue;
            // Sample BSDF for virtual light from _wi_
            uc = r();
            pstd::optional<BSDFSample> wis = exitInterface.Sample_f(
                wi, uc, Point2f(r(), r()), !mode, BxDFReflTransFlags::Transmission);
            if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0)
                continue;

            // Declare state for random walk through BSDF layers
            SampledSpectrum beta = wos->f * AbsCosTheta(wos->wi) / wos->pdf;
            Float z = enteredTop ? thickness : 0;
            // Debugging D.L.
            //  wis->wi = -wi;
            //  wos->wi = -wo;

            Vector3f w = wos->wi;

            // HGPhaseFunction phase(g);//Original code: Henyey-Greenstein phase function
            SGGXPhaseFunctionNormalDistribution phase(sigma, mean, std_or);
            // SGGXPhaseFunctionMixedRefl phase(sigma,mean,std_or,std_size);
            
            // phase.SetMicroflakes(dir, sigma);//Fiber-like distributions
            // phase.SetMicroflakes(30.f, 0.0f);//Fiber-like distributions

            // Diffuse phase function
            Point2f u1{r(), r()}, u2{r(), r()};
            bool is_specular_phase = use_specular;

            for (int depth = 0; depth < maxDepth; ++depth) {
                // std::cout<<"depth "<< depth << std::endl;
                // Sample next event for layered BSDF evaluation random walk
                PBRT_DBG("beta: %f %f %f %f, w: %f %f %f, f: %f %f %f %f\n", beta[0],
                         beta[1], beta[2], beta[3], w.x, w.y, w.z, f[0], f[1], f[2],
                         f[3]);
                // Possibly terminate layered BSDF random walk with Russian roulette
                //  if (depth > 3 && beta.MaxComponentValue() < 0.25f) {
                //  Float q = std::max<Float>(0, 1 - beta.MaxComponentValue());
                //  if (r() < q)
                //  break;
                // beta /= 1 - q;
                // PBRT_DBG("After RR with q = %f, beta: %f %f %f %f\n", q, beta[0],
                //  beta[1], beta[2], beta[3]);
                // }

                // Account for media between layers and possibly scatter
                if (!albedo) {
                    // Advance to next layer boundary and update _beta_ for transmittance
                    z = (z == thickness) ? 0 : thickness;
                    beta *= Tr(thickness, w);

                } else {
                    // Sample medium scattering for layered BSDF evaluation
                    Float sigma_t = 1;
                    // Take into account the projected Area of particles along direction w
                    Float pa = phase.projectedArea(-w, Point2f(r(), r()));
                    if (pa <= 0) {
                        std::cout << "Projected area is NOT positive" << std::endl;
                    }
                    sigma_t *= pa;
                    
                    // std::cout<<" pa "<< pa << std::endl;
                    // Float distSampled = (sigmaT == 0.f) ? Infinity: -std::log(1.f -
                    // rng.Uniform<Float>() ) / sigmaT;
                    Float dz = SampleExponential(r(), sigma_t / std::abs(w.z));

                    Float zp = w.z > 0 ? (z + dz) : (z - dz);
                    // std::cout<<"zp " << zp << std::endl;
                    DCHECK_RARE(1e-5, z == zp);
                    if (z == zp)
                        continue;
                    if (0 < zp && zp < thickness) {
                        // Handle scattering event in layered BSDF medium
                        // Account for scattering through _exitInterface_ using _wis_
                        Float wt = 1;
                        pstd::optional<PhaseFunctionSample> ps;
                        // std::cout<<"-wis->wi "<< -wis->wi <<" -w "<< -w <<std::endl;
                        // Specular phase function
                        if (is_specular_phase) {
                            if (!IsSpecular(exitInterface.Flags()))
                                wt = PowerHeuristic(
                                    1, wis->pdf, 1,
                                    phase.PDF(-w, -wis->wi, Point2f(r(), r())));
                            Float m_pdf = phase.PDF(-w, -wis->wi, Point2f(r(), r()));

                            // std::cout<<"mpdf "<< m_pdf<<std::endl;
                            if (m_pdf <= 0.0f) {
                                // std::cout<<" m_pdf is zero " << m_pdf << " -w " << -w
                                // << " -wis->wi " << -wis->wi<< " Good luck! "
                                // <<std::endl;
                            }
                            auto sigmaTAlongWi =
                                phase.projectedArea(-wis->wi, Point2f(r(), r()));
                            f += beta * albedo *
                                 phase.p(-w, -wis->wi, Point2f(r(), r())) * wt *
                                 Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /
                                 wis->pdf;

                            /*auto tdebug = beta * albedo * phase.p(-w,
                               -wis->wi,Point2f(r(),r())) * wt * Tr(zp - exitZ, wis->wi) *
                               wis->f / wis->pdf;*/

                            // Sample phase function and update layered path state
                            // ps = phase.Sample_p(-w, Point2f(r(),
                            // r()),Point2f(r(),r()));
                            ps = phase.Sample_p(-w, Point2f(r(), r()), Point2f(r(), r()));

                            // std::cout<<"tdebug " << tdebug << std::endl;
                            // std::cout<<"ps->wi " << ps->wi << std::endl;
                        } else {
                            // Diffuse phase function
                            if (!IsSpecular(exitInterface.Flags()))
                                wt = PowerHeuristic(
                                    1, wis->pdf, 1,
                                    phase.PDF(-w, -wis->wi, u1, Point2f(r(), r())));
                            auto sigmaTAlongWi =
                                phase.projectedArea(-wis->wi, Point2f(r(), r()));
                            f += beta * albedo *
                                 phase.p(-w, -wis->wi, Point2f(r(), r()),
                                         Point2f(r(), r())) *
                                 wt * Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /
                                 wis->pdf;
                            // f += beta * albedo * phase.p(-w, -wis->wi,
                            // u1,Point2f(r(),r())) * wt *
                            //      Tr(zp - exitZ, wis->wi) * wis->f / wis->pdf;

                            // Sample phase function and update layered path state
                            ps = phase.Sample_p(-w, Point2f(r(), r()), Point2f(r(), r()),
                                                Point2f(r(), r()));
                        }

                        if (!ps || ps->pdf == 0 || ps->wi.z == 0)
                            continue;

                        beta *= albedo * ps->p / ps->pdf;
                        w = ps->wi;
                        z = zp;

                        // Possibly account for scattering through _exitInterface_
                        if (((z < exitZ && w.z > 0) || (z > exitZ && w.z < 0)) &&
                            !IsSpecular(exitInterface.Flags())) {
                            // Account for scattering through _exitInterface_
                            SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                            if (fExit) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                Float wt = PowerHeuristic(1, ps->pdf, 1, exitPDF);
                                auto sigmaTAlongWi =
                                    phase.projectedArea(-ps->wi, Point2f(r(), r()));

                                f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi) *
                                     fExit * wt;

                                // f += beta * Tr(zp - exitZ, ps->wi) * fExit * wt;
                            }
                        }

                        continue;
                    }
                    z = Clamp(zp, 0, thickness);
                }

                // Account for scattering at appropriate interface
                if (z == exitZ) {
                    // Account for reflection at _exitInterface_
                    Float uc = r();
                    pstd::optional<BSDFSample> bs = exitInterface.Sample_f(
                        -w, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        break;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;
                    w = bs->wi;

                } else {
                    // Account for scattering at _nonExitInterface_
                    if (!IsSpecular(nonExitInterface.Flags())) {
                        // Add NEE contribution along presampled _wis_ direction
                        Float wt = 1;
                        if (!IsSpecular(exitInterface.Flags()))
                            wt = PowerHeuristic(1, wis->pdf, 1,
                                                nonExitInterface.PDF(-w, -wis->wi, mode));

                        // f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                        //  AbsCosTheta(wis->wi) * wt * Tr(thickness, wis->wi) * wis->f /
                        //  wis->pdf;
                        auto sigmaTAlongWi =
                            phase.projectedArea(-wis->wi, Point2f(r(), r()));
                        f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                             AbsCosTheta(wis->wi) * wt *
                             Tr(thickness, wis->wi, sigmaTAlongWi) * wis->f / wis->pdf;
                    }
                    // Sample new direction using BSDF at _nonExitInterface_
                    Float uc = r();
                    Point2f u(r(), r());
                    pstd::optional<BSDFSample> bs = nonExitInterface.Sample_f(
                        -w, uc, u, mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        break;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;
                    w = bs->wi;

                    if (!IsSpecular(exitInterface.Flags())) {
                        // Add NEE contribution along direction from BSDF sample
                        SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                        if (fExit) {
                            Float wt = 1;
                            if (!IsSpecular(nonExitInterface.Flags())) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                wt = PowerHeuristic(1, bs->pdf, 1, exitPDF);
                            }
                            auto sigmaTAlongWi =
                                phase.projectedArea(-bs->wi, Point2f(r(), r()));

                            f += beta * Tr(thickness, bs->wi, sigmaTAlongWi) * fExit * wt;
                        }
                    }
                }
            }
        }

        return f / nSamples;
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(
        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // Set _wo_ for layered BSDF sampling
        bool flipWi = false;
        if (twoSided && wo.z < 0) {
            wo = -wo;
            flipWi = true;
        }

        // Sample BSDF at entrance interface to get initial direction _w_
        bool enteredTop = twoSided || wo.z > 0;
        pstd::optional<BSDFSample> bs =
            enteredTop ? top.Sample_f(wo, uc, u, mode) : bottom.Sample_f(wo, uc, u, mode);
        

        if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
            return {};
        if (bs->IsReflection()) {
            if (flipWi)
                bs->wi = -bs->wi;
            bs->pdfIsProportional = true;
            return bs;
        }
        Vector3f w = bs->wi;
        bool specularPath = bs->IsSpecular();

        // Declare _RNG_ for layered BSDF sampling
        RNG rng(Hash(GetOptions().seed, wo), Hash(uc, u));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Declare common variables for layered BSDF sampling
        SampledSpectrum f = bs->f * AbsCosTheta(bs->wi);
        Float pdf = bs->pdf;
        Float z = enteredTop ? thickness : 0;
        // HGPhaseFunction phase(g);//Original code: Henyey-Greenstein phase function
        //  SGGXPhaseFunctionMixedRefl phase(sigma,mean,std_or,std_size);
        SGGXPhaseFunctionNormalDistribution phase(sigma, mean, std_or);

        // phase.SetMicroflakes(dir, sigma);//Fiber-like distributions

        // Diffuse phase function
        Point2f u1{r(), r()}, u2{r(), r()};
        bool is_specular_phase = use_specular;

        for (int depth = 0; depth < maxDepth; ++depth) {
            // Follow random walk through layers to sample layered BSDF
            // Possibly terminate layered BSDF sampling with Russian Roulette
            Float rrBeta = f.MaxComponentValue() / pdf;
            if (depth > 3 && rrBeta < 0.25f) {
            Float q = std::max<Float>(0, 1 - rrBeta);
            if (r() < q)
            return {};
            pdf *= 1 - q;
            }
            if (w.z == 0)
                return {};

            if (albedo) {
                // Sample potential scattering event in layered medium
                Float sigma_t = 1;
                Float pa = phase.projectedArea(w, Point2f(r(), r()));
                sigma_t *= pa;

                Float dz = SampleExponential(r(), sigma_t / AbsCosTheta(w));
                Float zp = w.z > 0 ? (z + dz) : (z - dz);
                CHECK_RARE(1e-5, zp == z);
                if (zp == z)
                    return {};
                if (0 < zp && zp < thickness) {
                    // Update path state for valid scattering event between interfaces

                    pstd::optional<PhaseFunctionSample> ps;

                    if (is_specular_phase) {
                        ps =
                            phase.Sample_p(-w, Point2f(r(), r()),
                                           Point2f(r(), r()));  // Specular phase function
                    } else {
                        ps = phase.Sample_p(-w, Point2f(r(), r()), Point2f(r(), r()),
                                            Point2f(r(), r()));  // Diffuse phase function
                    }

                    if (!ps || ps->pdf == 0 || ps->wi.z == 0)
                        return {};
                    f *= albedo * ps->p;
                    pdf *= ps->pdf;
                    specularPath = false;
                    w = ps->wi;
                    z = zp;

                    continue;
                }
                z = Clamp(zp, 0, thickness);
                if (z == 0)
                    DCHECK_LT(w.z, 0);
                else
                    DCHECK_GT(w.z, 0);

            } else {
                // Advance to the other layer interface
                z = (z == thickness) ? 0 : thickness;
                f *= Tr(thickness, w);
            }
            // Initialize _interface_ for current interface surface
#ifdef interface  // That's enough out of you, Windows.
#undef interface
#endif
            TopOrBottomBxDF<TopBxDF, BottomBxDF> interface;
            if (z == 0)
                interface = &bottom;
            else
                interface = &top;

            // Sample interface BSDF to determine new path direction
            Float uc = r();
            Point2f u(r(), r());
            pstd::optional<BSDFSample> bs = interface.Sample_f(-w, uc, u, mode);
            if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                return {};
            f *= bs->f;
            pdf *= bs->pdf;
            specularPath &= bs->IsSpecular();
            w = bs->wi;

            // Return _BSDFSample_ if path has left the layers
            if (bs->IsTransmission()) {
                BxDFFlags flags = SameHemisphere(wo, w) ? BxDFFlags::Reflection
                                                        : BxDFFlags::Transmission;
                flags |= specularPath ? BxDFFlags::Specular : BxDFFlags::Glossy;
                if (flipWi)
                    w = -w;

                return BSDFSample(f, w, pdf, flags, 1.f, true);
            }

            // Scale _f_ by cosine term after scattering at the interface
            f *= AbsCosTheta(bs->wi);
        }
        return {};
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // Set _wo_ and _wi_ for layered BSDF evaluation
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }

        // Declare _RNG_ for layered PDF evaluation
        RNG rng(Hash(GetOptions().seed, wi), Hash(wo));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Update _pdfSum_ for reflection at the entrance layer
        bool enteredTop = twoSided || wo.z > 0;
        Float pdfSum = 0;
        if (SameHemisphere(wo, wi)) {
            auto reflFlag = BxDFReflTransFlags::Reflection;
            pdfSum += enteredTop ? nSamples * top.PDF(wo, wi, mode, reflFlag)
                                 : nSamples * bottom.PDF(wo, wi, mode, reflFlag);
        }

        for (int s = 0; s < nSamples; ++s) {
            // Evaluate layered BSDF PDF sample
            if (SameHemisphere(wo, wi)) {
                // Evaluate TRT term for PDF estimate
                TopOrBottomBxDF<TopBxDF, BottomBxDF> rInterface, tInterface;
                if (enteredTop) {
                    rInterface = &bottom;
                    tInterface = &top;
                } else {
                    rInterface = &top;
                    tInterface = &bottom;
                }
                // Sample _tInterface_ to get direction into the layers
                auto trans = BxDFReflTransFlags::Transmission;
                pstd::optional<BSDFSample> wos, wis;
                wos = tInterface.Sample_f(wo, r(), {r(), r()}, mode, trans);
                wis = tInterface.Sample_f(wi, r(), {r(), r()}, !mode, trans);

                // Update _pdfSum_ accounting for TRT scattering events
                if (wos && wos->f && wos->pdf > 0 && wis && wis->f && wis->pdf > 0) {
                    if (!IsNonSpecular(tInterface.Flags()))
                        pdfSum += rInterface.PDF(-wos->wi, -wis->wi, mode);
                    else {
                        // Use multiple importance sampling to estimate PDF product
                        pstd::optional<BSDFSample> rs =
                            rInterface.Sample_f(-wos->wi, r(), {r(), r()}, mode);
                        if (rs && rs->f && rs->pdf > 0) {
                            if (!IsNonSpecular(rInterface.Flags()))
                                pdfSum += tInterface.PDF(-rs->wi, wi, mode);
                            else {
                                // Compute MIS-weighted estimate of Equation
                                // (\ref{eq:pdf-triple-canceled-one})
                                Float rPDF = rInterface.PDF(-wos->wi, -wis->wi, mode);
                                Float wt = PowerHeuristic(1, wis->pdf, 1, rPDF);
                                pdfSum += wt * rPDF;

                                Float tPDF = tInterface.PDF(-rs->wi, wi, mode);
                                wt = PowerHeuristic(1, rs->pdf, 1, tPDF);
                                pdfSum += wt * tPDF;
                            }
                        }
                    }
                }

            } else {
                // Evaluate TT term for PDF estimate
                TopOrBottomBxDF<TopBxDF, BottomBxDF> toInterface, tiInterface;
                if (enteredTop) {
                    toInterface = &top;
                    tiInterface = &bottom;
                } else {
                    toInterface = &bottom;
                    tiInterface = &top;
                }

                Float uc = r();
                Point2f u(r(), r());
                pstd::optional<BSDFSample> wos = toInterface.Sample_f(wo, uc, u, mode);
                if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0 ||
                    wos->IsReflection())
                    continue;

                uc = r();
                u = Point2f(r(), r());
                pstd::optional<BSDFSample> wis = tiInterface.Sample_f(wi, uc, u, !mode);
                if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0 ||
                    wis->IsReflection())
                    continue;

                if (IsSpecular(toInterface.Flags()))
                    pdfSum += tiInterface.PDF(-wos->wi, wi, mode);
                else if (IsSpecular(tiInterface.Flags()))
                    pdfSum += toInterface.PDF(wo, -wis->wi, mode);
                else
                    pdfSum += (toInterface.PDF(wo, -wis->wi, mode) +
                               tiInterface.PDF(-wos->wi, wi, mode)) /
                              2;
            }
        }
        // Return mixture of PDF estimate and constant PDF
        return Lerp(0.9f, 1 / (4 * Pi), pdfSum / nSamples);
    }

  private:
    // MyLayeredBxDF Private Methods
    PBRT_CPU_GPU
    static Float Tr(Float dz, Vector3f w) {
        if (std::abs(dz) <= std::numeric_limits<Float>::min())
            return 1;
        return FastExp(-std::abs(dz / w.z));
    }
    static Float Tr(Float dz, Vector3f w, Float sigmaTAlongWi) {
        return FastExp(-std::abs(sigmaTAlongWi * dz / w.z));
    }

    // MyLayeredBxDF Private Members
    TopBxDF top;
    BottomBxDF bottom;
    Float thickness, g, sigma, mean, std_or, std_size;
    Vector3f dir;
    SampledSpectrum albedo;
    int maxDepth, nSamples;
    bool use_specular;
};

template <typename TopBxDF, typename BottomBxDF, bool twoSided>
class LayeredBxDFMixedParticles {
  public:
    // LayeredBxDFMixedParticles Public Methods
    LayeredBxDFMixedParticles() = default;
    PBRT_CPU_GPU
    LayeredBxDFMixedParticles(TopBxDF top, BottomBxDF bottom, Float thickness,
                              Vector3f dir, Float sigma, Float mean, Float std_or,
                              Float std_size, Float diff_conc, bool use_refl,
                              const SampledSpectrum &albedo_flakes,const SampledSpectrum &albedo_diff, Float g1,Float g2,Float w_g, int maxDepth,
                              int nSamples)
        : top(top),
          bottom(bottom),
          thickness(std::max(thickness, std::numeric_limits<Float>::min())),
          dir(dir),
          sigma(sigma),
          mean(mean),
          std_or(std_or),
          std_size(std_size),
          use_specular(use_refl),
          g1(g1),
          g2(g2),
          w_g(w_g),
          albedo_flakes(albedo_flakes),
          albedo_diff(albedo_diff),
          maxDepth(maxDepth),
          nSamples(nSamples),
          diff_conc(diff_conc),
          sigma_t(1.0f) {
            
        // std::cout<<"Mean is " << mean << std::endl;
        // std::cout<<"std SIZE is " << std_size << std::endl;
    }

    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() {
        top.Regularize();
        bottom.Regularize();
    }

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        //TODO - remove this        
        // SampledSpectrum albedo(0.0f);
        BxDFFlags topFlags = top.Flags(), bottomFlags = bottom.Flags();
        CHECK(IsTransmissive(topFlags) ||
              IsTransmissive(bottomFlags));  // otherwise, why bother?

        BxDFFlags flags = BxDFFlags::Reflection;
        if (IsSpecular(topFlags))
            flags = flags | BxDFFlags::Specular;

        if (IsDiffuse(topFlags) || IsDiffuse(bottomFlags) || albedo_diff || albedo_flakes)
            flags = flags | BxDFFlags::Diffuse;
        else if (IsGlossy(topFlags) || IsGlossy(bottomFlags))
            flags = flags | BxDFFlags::Glossy;

        if (IsTransmissive(topFlags) && IsTransmissive(bottomFlags))
            flags = flags | BxDFFlags::Transmission;

        return flags;
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        SampledSpectrum w_albedo2;
         w_albedo2 = diff_conc * albedo_diff + (1-diff_conc) * albedo_flakes;
                            // if (ps->p_flag == _MediumFlags::Diffuser){
                            //     albedo_particle = albedo_diff;
                            // }else{
                            //     albedo_particle = albedo_flakes;
                            // }
        SampledSpectrum f(0.);
        // Next line is for debugging purposes only, D.L.
        //  wo = Vector3f(0,0,1);
        //  wi = Vector3f(0,0,1);
        //  Estimate _LayeredBxDFMixedParticles_ value _f_ using random sampling
        //  Set _wo_ and _wi_ for layered BSDF evaluation
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }
        
        // Determine entrance interface for layered BSDF
        TopOrBottomBxDF<TopBxDF, BottomBxDF> enterInterface;
        bool enteredTop = twoSided || wo.z > 0;
        if (enteredTop)
            enterInterface = &top;
        else
            enterInterface = &bottom;
            

        // Determine exit interface and exit $z$ for layered BSDF
        TopOrBottomBxDF<TopBxDF, BottomBxDF> exitInterface, nonExitInterface;
        if (SameHemisphere(wo, wi) ^ enteredTop) {
            
            exitInterface = &bottom;
            nonExitInterface = &top;
        } else {
            exitInterface = &top;
            nonExitInterface = &bottom;
        }
        Float exitZ = (SameHemisphere(wo, wi) ^ enteredTop) ? 0 : thickness;

        // Account for reflection at the entrance interface
        if (SameHemisphere(wo, wi))
            f = nSamples * enterInterface.f(wo, wi, mode);

        // Declare _RNG_ for layered BSDF evaluation
        RNG rng(Hash(GetOptions().seed, wo), Hash(wi));

        
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };
        // std::cout << "wi "<< wi << std::endl;
        // std::cout << "wo "<< wo << std::endl;
        // auto t_seed = GetOptions().seed;
        // std::cout<< t_seed << std::endl;
        // std::cout<<"first r " << r() <<std::endl;
        for (int s = 0; s < nSamples; ++s) {
            // Sample random walk through layers to estimate BSDF value
            // Sample transmission direction through entrance interface
            Float uc = r();
            pstd::optional<BSDFSample> wos = enterInterface.Sample_f(
                wo, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Transmission);
            if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0)
                continue;

            // Sample BSDF for virtual light from _wi_
            uc = r();
            pstd::optional<BSDFSample> wis = exitInterface.Sample_f(
                wi, uc, Point2f(r(), r()), !mode, BxDFReflTransFlags::Transmission);
            if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0)
                continue;
            // std::cout<<"wos->wi" << wos->wi <<" wi " << wi << std::endl;
            // std::cout<<"wis->wi" << wis->wi <<" wo " << wo << std::endl;



            // Declare state for random walk through BSDF layers
            SampledSpectrum beta = wos->f * AbsCosTheta(wos->wi) / wos->pdf;
            Float z = enteredTop ? thickness : 0;
            // Debugging D.L.
            //  wis->wi = -wi;
            //  wos->wi = -wo;
            // std::cout<< " " << AbsCosTheta(wos->wi) << "  " << AbsCosTheta(wo) << std::endl;
            
            Vector3f w = wos->wi;

            // SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g);
            SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g,albedo_flakes,albedo_diff);
            
            // SGGXPhaseFunctionNormalDistribution phase(sigma,mean,std_or);

            // Diffuse phase function
            Point2f u1{r(), r()}, u2{r(), r()};
            bool is_specular_phase = use_specular;
            // SampledSpectrum albedo_particle;

            for (int depth = 0; depth < maxDepth; ++depth) {
                // Sample next event for layered BSDF evaluation random walk
                PBRT_DBG("beta: %f %f %f %f, w: %f %f %f, f: %f %f %f %f\n", beta[0],
                         beta[1], beta[2], beta[3], w.x, w.y, w.z, f[0], f[1], f[2],
                         f[3]);
                // Possibly terminate layered BSDF random walk with Russian roulette
                 if (depth > 128 && beta.MaxComponentValue() < 0.01f) {
                 Float q = std::max<Float>(0, 1 - beta.MaxComponentValue());
                 if (r() < q)
                 break;
                beta /= 1 - q;
                PBRT_DBG("After RR with q = %f, beta: %f %f %f %f\n", q, beta[0],
                 beta[1], beta[2], beta[3]);
                }
                // Account for media between layers and possibly scatter
                // if (!true) {
                    //line after was the usual one but in any case this if-else is not used
                if (!true) {
                    // Advance to next layer boundary and update _beta_ for transmittance
                    z = (z == thickness) ? 0 : thickness;
                    beta *= Tr(thickness, w);
                } else {
                    // Sample medium scattering for layered BSDF evaluation
                    // Take into account the projected Area of particles along direction w
                    Float pa = phase.projectedArea(-w, Point2f(r(), r()));
                    //TODO:
                    //change this put the weighting factor
                    // sigma_t *= pa; //Only for SGGX
                    // sigma_t = diff_conc * sigma_t + (1-diff_conc)*(sigma_t * pa) ;
                    Float reduced_st = diff_conc * sigma_t + (1-diff_conc)*(sigma_t * pa) ;
                    
                    if (reduced_st <= 0) {
                        std::cout << "Projected area is NOT positive" << std::endl;
                        continue;
                    }
                    // std::cout<<" pa "<< pa << std::endl;
                    // Float distSampled = (sigmaT == 0.f) ? Infinity: -std::log(1.f -
                    // rng.Uniform<Float>() ) / sigmaT;
                    Float dz = SampleExponential(r(), reduced_st / std::abs(w.z));

                    Float zp = w.z > 0 ? (z + dz) : (z - dz);
                    // std::cout<<"zp " << zp << std::endl;
                    DCHECK_RARE(1e-5, z == zp);
                    if (z == zp)
                        continue;
                    if (0 < zp && zp < thickness) {
                        // Handle scattering event in layered BSDF medium
                        // Account for scattering through _exitInterface_ using _wis_
                        Float wt = 1;
                        pstd::optional<PhaseFunctionSample> ps;
                        // std::cout<<"-wis->wi "<< -wis->wi <<" -w "<< -w <<std::endl;
                        // Specular phase function
                        // if (is_specular_phase) {
                            if (!IsSpecular(exitInterface.Flags()))
                                wt = PowerHeuristic(
                                    1, wis->pdf, 1,
                                    phase.PDF(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()) ) );
                            Float m_pdf = phase.PDF(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));

                            // std::cout<<"mpdf "<< m_pdf<<std::endl;
                            //if (m_pdf <= 0.0f) {
                                // std::cout<<" m_pdf is zero " << m_pdf << " -w " << -w
                                // << " -wis->wi " << -wis->wi<< " Good luck! "
                                // <<std::endl;
                            //}
                            // auto sigmaTAlongWi = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                            pa = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                            auto sigmaTAlongWi = diff_conc  + (1-diff_conc)*( pa) ;
                            
                            SampledSpectrum old_f = f;
                            Float pf_p = phase.p(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));
                            SampledSpectrum w_albedo = phase.albedoW(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));
                            // std::cout<<"w_albedo"<< w_albedo << std::endl;
                            f += beta * w_albedo *  wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            // f += beta * w_albedo * pf_p * wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            // f += beta * pf_p * wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            
                            SampledSpectrum albedo_particle;
                            ps = phase.Sample_p(-w, Point2f(r(), r()),Point2f(r(), r()),r(),Point2f(r(), r()));
                            if (ps->p_flag == _MediumFlags::Diffuser){
                                albedo_particle = albedo_diff;
                                // std::cout<<"albedo_diff "<< albedo_diff << std::endl;
                            }else{
                                albedo_particle = albedo_flakes;
                                // std::cout<<"albedo_flakes "<< albedo_flakes << std::endl;

                            }


                            if (IsNaN(f[0])) {
                                std::cout << " nan found in f()_s_05" << std::endl;
                                std::cout << "oldF " << old_f << " pf_p " << pf_p
                                          << std::endl;
                                std::cout << "beta " << beta << " albedo " << w_albedo
                                          << " wt " << wt << " wis->f " << wis->f
                                          << " wis->pdf" << wis->pdf << std::endl;
                                std::cout << "Tr "
                                          << Tr(zp - exitZ, wis->wi, sigmaTAlongWi)
                                          << std::endl;
                                // std::cout<<" m_pdf is zero " << m_pdf << " -w " << -w
                                // << " -wis->wi " << -wis->wi<< " Good luck! "
                                // <<std::endl;
                            }
     

                        if (!ps || ps->pdf == 0 || ps->wi.z == 0 || ps->p == 0.0f)
                            continue;
                        // original line
                        // beta *= albedo * ps->p / ps->pdf;
                        beta *= albedo_particle;
                        
                        //beta *=  w_albedo/pf_p;
                        // beta *=  w_albedo;
                        w = ps->wi;
                        z = zp;

                        // Possibly account for scattering through _exitInterface_
                        if (((z < exitZ && w.z > 0) || (z > exitZ && w.z < 0)) &&
                            !IsSpecular(exitInterface.Flags())) {
                            // Account for scattering through _exitInterface_
                            SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                            if (fExit) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                Float wt = PowerHeuristic(1, ps->pdf, 1, exitPDF);
                                // auto sigmaTAlongWi = phase.projectedArea(-ps->wi, Point2f(r(), r()));
                                pa = phase.projectedArea(-ps->wi, Point2f(r(), r()));
                                sigmaTAlongWi = diff_conc * 1 + (1-diff_conc)*(1 * pa) ;

                                f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi) *
                                     fExit * wt;

                                // f += beta * Tr(zp - exitZ, ps->wi) * fExit * wt;
                            }
                        }

                        continue;
                    }
                    z = Clamp(zp, 0, thickness);
                }

                // Account for scattering at appropriate interface
                if (z == exitZ) {
                    // Account for reflection at _exitInterface_
                    Float uc = r();
                    pstd::optional<BSDFSample> bs = exitInterface.Sample_f(
                        -w, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        break;
                    // std::cout << "Reflection back: "<< bs->wi << " w "<< w <<  std::endl;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;
                    w = bs->wi;

                } else {
                    // Account for scattering at _nonExitInterface_
                    if (!IsSpecular(nonExitInterface.Flags())) {
                        // Add NEE contribution along presampled _wis_ direction
                        Float wt = 1;
                        if (!IsSpecular(exitInterface.Flags()))
                            wt = PowerHeuristic(1.0f, wis->pdf, 1.0f,
                                                nonExitInterface.PDF(-w, -wis->wi, mode));

                        // f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                        //  AbsCosTheta(wis->wi) * wt * Tr(thickness, wis->wi) * wis->f /
                        //  wis->pdf;
                        // auto sigmaTAlongWi = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                        auto pa = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                        auto sigmaTAlongWi = diff_conc  + (1.0f-diff_conc)*(pa) ;
                        
                        f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                             AbsCosTheta(wis->wi) * wt *
                             Tr(thickness, wis->wi, sigmaTAlongWi) * wis->f / wis->pdf;
                    }
                    // Sample new direction using BSDF at _nonExitInterface_
                    Float uc = r();
                    Point2f u(r(), r());
                    pstd::optional<BSDFSample> bs = nonExitInterface.Sample_f(
                        -w, uc, u, mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        break;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;

                    w = bs->wi;

                    if (!IsSpecular(exitInterface.Flags())) {
                        // Add NEE contribution along direction from BSDF sample
                        SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                        if (fExit) {
                            Float wt = 1;
                            if (!IsSpecular(nonExitInterface.Flags())) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                wt = PowerHeuristic(1, bs->pdf, 1, exitPDF);
                            }
                            auto pa = phase.projectedArea(-bs->wi, Point2f(r(), r()));
                            // auto sigmaTAlongWi = phase.projectedArea(-bs->wi, Point2f(r(), r()));
                            auto sigmaTAlongWi = diff_conc * 1 + (1-diff_conc)*(1 * pa) ;

                            f += beta * Tr(thickness, bs->wi, sigmaTAlongWi) * fExit * wt;
                        }
                    }
                }
            }
        }
        if (f.HasNaNs()) {
            std::cout << " nan found in f()" << std::endl;
        }
        return f / nSamples;
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(
        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // Set _wo_ for layered BSDF sampling
        bool flipWi = false;
        if (twoSided && wo.z < 0) {
            wo = -wo;
            flipWi = true;
        }

        // Sample BSDF at entrance interface to get initial direction _w_
        bool enteredTop = twoSided || wo.z > 0;
        pstd::optional<BSDFSample> bs =
            enteredTop ? top.Sample_f(wo, uc, u, mode) : bottom.Sample_f(wo, uc, u, mode);
        if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
            return {};
        if (bs->IsReflection()) {
            if (flipWi)
                bs->wi = -bs->wi;
            bs->pdfIsProportional = true;
            return bs;
        }
        Vector3f w = bs->wi;
        bool specularPath = bs->IsSpecular();

        // Declare _RNG_ for layered BSDF sampling
        RNG rng(Hash(GetOptions().seed, wo), Hash(uc, u));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Declare common variables for layered BSDF sampling
        SampledSpectrum f = bs->f * AbsCosTheta(bs->wi);
                        if (IsInf(f[0])){
                    std::cout<< "it's INF line 1900"<< std::endl;
                }
        Float pdf = bs->pdf;
        Float z = enteredTop ? thickness : 0;
        // HGPhaseFunction phase(g);//Original code: Henyey-Greenstein phase function
        SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g,albedo_flakes,albedo_diff);
        bool is_specular_phase = use_specular;
        
        for (int depth = 0; depth < maxDepth; ++depth) {
            // Follow random walk through layers to sample layered BSDF
            // Possibly terminate layered BSDF sampling with Russian Roulette
            Float rrBeta = f.MaxComponentValue() / pdf;
            if (depth > 128 && rrBeta < 0.01f) {
            Float q = std::max<Float>(0, 1 - rrBeta);
            if (r() < q)
            return {};
            pdf *= 1 - q;
            }
            if (w.z == 0)
                return {};

            if (albedo_diff || albedo_flakes) {
                // Sample potential scattering event in layered medium
                Float pa = phase.projectedArea(w, Point2f(r(), r()));
                
                Float reduced_st = diff_conc * sigma_t + (1-diff_conc)*(sigma_t * pa) ;

                Float dz = SampleExponential(r(), reduced_st / AbsCosTheta(w));
                Float zp = w.z > 0 ? (z + dz) : (z - dz);
                CHECK_RARE(1e-5, zp == z);
                if (zp == z)
                    return {};
                if (0 < zp && zp < thickness) {
                    // Update path state for valid scattering event between interfaces

                    pstd::optional<PhaseFunctionSample> ps;
                    //get the albedo of the sampled particle
                    SampledSpectrum albedo_particle;
                    ps = phase.Sample_p(-w, Point2f(r(), r()),Point2f(r(), r()),r(),Point2f(r(), r())); 
                    if (ps->p_flag == _MediumFlags::Diffuser){
                                
                                albedo_particle = albedo_diff;
                            }else{
                                // std::cout<<"happens"<<std::endl;
                                albedo_particle = albedo_flakes * pa;
                    }
                    if (!ps || ps->pdf == 0 || ps->wi.z == 0){
                        // LOG_ERROR("early termination is not great",pPixel.x, pPixel.y, sampleIndex);
                        return {};

                    }
                    //auto old_f = f;
                    f *= albedo_particle;    
                    specularPath = false;
                    w = ps->wi;
                    z = zp;

                    continue;
                }
                z = Clamp(zp, 0, thickness);
                if (z == 0)
                    DCHECK_LT(w.z, 0);
                else
                    DCHECK_GT(w.z, 0);

            } else {
                // Advance to the other layer interface
                z = (z == thickness) ? 0 : thickness;
                f *= Tr(thickness, w);
            }
            // Initialize _interface_ for current interface surface
#ifdef interface  // That's enough out of you, Windows.
#undef interface
#endif
            TopOrBottomBxDF<TopBxDF, BottomBxDF> interface;
            if (z == 0)
                interface = &bottom;
            else
                interface = &top;

            // Sample interface BSDF to determine new path direction
            Float uc = r();
            Point2f u(r(), r());
            pstd::optional<BSDFSample> bs = interface.Sample_f(-w, uc, u, mode);
            if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0){
                        // LOG_ERROR("early termination is not great");
                return {};

            }
            f *= bs->f;
            if (IsInf(f[0])){
                    std::cout<< "it's INF line 1990"<< std::endl;
            }
            pdf *= bs->pdf;
            specularPath &= bs->IsSpecular();
            w = bs->wi;

            // Return _BSDFSample_ if path has left the layers
            if (bs->IsTransmission()) {
                BxDFFlags flags = SameHemisphere(wo, w) ? BxDFFlags::Reflection
                                                        : BxDFFlags::Transmission;
                flags |= specularPath ? BxDFFlags::Specular : BxDFFlags::Glossy;
                if (flipWi)
                    w = -w;
                if (f.HasNaNs()) {
                    std::cout << " nan found in sample()" << std::endl;
                }
                if (IsInf(f[0])){
                    std::cout<< "it's INF"<< std::endl;
                }
                return BSDFSample(f, w, pdf, flags, 1.f, true);
            }

            // Scale _f_ by cosine term after scattering at the interface
            f *= AbsCosTheta(bs->wi);
        }
        // LOG_ERROR("Reaching End");
        return {};
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // Set _wo_ and _wi_ for layered BSDF evaluation
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }

        // Declare _RNG_ for layered PDF evaluation
        RNG rng(Hash(GetOptions().seed, wi), Hash(wo));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Update _pdfSum_ for reflection at the entrance layer
        bool enteredTop = twoSided || wo.z > 0;
        Float pdfSum = 0;
        if (SameHemisphere(wo, wi)) {
            auto reflFlag = BxDFReflTransFlags::Reflection;
            pdfSum += enteredTop ? nSamples * top.PDF(wo, wi, mode, reflFlag)
                                 : nSamples * bottom.PDF(wo, wi, mode, reflFlag);
        }

        for (int s = 0; s < nSamples; ++s) {
            // Evaluate layered BSDF PDF sample
            if (SameHemisphere(wo, wi)) {
                // Evaluate TRT term for PDF estimate
                TopOrBottomBxDF<TopBxDF, BottomBxDF> rInterface, tInterface;
                if (enteredTop) {
                    rInterface = &bottom;
                    tInterface = &top;
                } else {
                    rInterface = &top;
                    tInterface = &bottom;
                }
                // Sample _tInterface_ to get direction into the layers
                auto trans = BxDFReflTransFlags::Transmission;
                pstd::optional<BSDFSample> wos, wis;
                wos = tInterface.Sample_f(wo, r(), {r(), r()}, mode, trans);
                wis = tInterface.Sample_f(wi, r(), {r(), r()}, !mode, trans);

                // Update _pdfSum_ accounting for TRT scattering events
                if (wos && wos->f && wos->pdf > 0 && wis && wis->f && wis->pdf > 0) {
                    if (!IsNonSpecular(tInterface.Flags()))
                        pdfSum += rInterface.PDF(-wos->wi, -wis->wi, mode);
                    else {
                        // Use multiple importance sampling to estimate PDF product
                        pstd::optional<BSDFSample> rs =
                            rInterface.Sample_f(-wos->wi, r(), {r(), r()}, mode);
                        if (rs && rs->f && rs->pdf > 0) {
                            if (!IsNonSpecular(rInterface.Flags()))
                                pdfSum += tInterface.PDF(-rs->wi, wi, mode);
                            else {
                                // Compute MIS-weighted estimate of Equation
                                // (\ref{eq:pdf-triple-canceled-one})
                                Float rPDF = rInterface.PDF(-wos->wi, -wis->wi, mode);
                                Float wt = PowerHeuristic(1, wis->pdf, 1, rPDF);
                                pdfSum += wt * rPDF;

                                Float tPDF = tInterface.PDF(-rs->wi, wi, mode);
                                wt = PowerHeuristic(1, rs->pdf, 1, tPDF);
                                pdfSum += wt * tPDF;
                            }
                        }
                    }
                }

            } else {
                // Evaluate TT term for PDF estimate
                TopOrBottomBxDF<TopBxDF, BottomBxDF> toInterface, tiInterface;

                if (enteredTop) {
                    toInterface = &top;
                    tiInterface = &bottom;
                } else {
                    toInterface = &bottom;
                    tiInterface = &top;
                }

                Float uc = r();
                Point2f u(r(), r());
                pstd::optional<BSDFSample> wos = toInterface.Sample_f(wo, uc, u, mode);
                if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0 ||
                    wos->IsReflection())
                    continue;

                uc = r();
                u = Point2f(r(), r());
                pstd::optional<BSDFSample> wis = tiInterface.Sample_f(wi, uc, u, !mode);
                if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0 ||
                    wis->IsReflection())
                    continue;

                if (IsSpecular(toInterface.Flags()))
                    pdfSum += tiInterface.PDF(-wos->wi, wi, mode);
                else if (IsSpecular(tiInterface.Flags()))
                    pdfSum += toInterface.PDF(wo, -wis->wi, mode);
                else
                    pdfSum += (toInterface.PDF(wo, -wis->wi, mode) +
                               tiInterface.PDF(-wos->wi, wi, mode)) /
                              2;
            }
        }
        // Return mixture of PDF estimate and constant PDF
        return Lerp(0.9f, 1 / (4 * Pi), pdfSum / nSamples);
    }

  private:
    // LayeredBxDFMixedParticles Private Methods
    PBRT_CPU_GPU
    static Float Tr(Float dz, Vector3f w) {
        if (std::abs(dz) <= std::numeric_limits<Float>::min())
            return 1;
        return FastExp(-std::abs(dz / w.z));
    }
    static Float Tr(Float dz, Vector3f w, Float sigmaTAlongWi) {
        return FastExp(-std::abs(sigmaTAlongWi * dz / w.z));
    }

    // LayeredBxDFMixedParticles Private Members
    TopBxDF top;
    BottomBxDF bottom;
    Float thickness, g1,g2,w_g, sigma, mean, std_or, std_size, diff_conc;
    Float sigma_t;

    Vector3f dir;
    SampledSpectrum albedo_diff,albedo_flakes;
    int maxDepth, nSamples;
    bool use_specular;
};

template <typename TopBxDF, typename BottomBxDF, bool twoSided>
class LayeredBxDFMixedParticlesTransmittance {
  public:
    // LayeredBxDFMixedParticlesTransmittance Public Methods
    LayeredBxDFMixedParticlesTransmittance() = default;
    PBRT_CPU_GPU
    LayeredBxDFMixedParticlesTransmittance(TopBxDF top, BottomBxDF bottom, Float thickness,
                              Vector3f dir, Float sigma, Float mean, Float std_or,
                              Float std_size, Float diff_conc, bool use_refl,
                              const SampledSpectrum &albedo_flakes,const SampledSpectrum &albedo_diff, Float g1,Float g2,Float w_g, int maxDepth,
                              int nSamples)
        : top(top),
          bottom(bottom),
          thickness(std::max(thickness, std::numeric_limits<Float>::min())),
          dir(dir),
          sigma(sigma),
          mean(mean),
          std_or(std_or),
          std_size(std_size),
          use_specular(use_refl),
          g1(g1),
          g2(g2),
          w_g(w_g),
          albedo_flakes(albedo_flakes),
          albedo_diff(albedo_diff),
          maxDepth(maxDepth),
          nSamples(nSamples),
          diff_conc(diff_conc),
          sigma_t(1.0f),
          eta(1.3f) {
            
        // std::cout<<"std SIZE is " << std_size << std::endl;
    }

    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() {
        top.Regularize();
        bottom.Regularize();
    }

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        //TODO - remove this        
        // SampledSpectrum albedo(0.0f);
        BxDFFlags topFlags = top.Flags(), bottomFlags = bottom.Flags();
        CHECK(IsTransmissive(topFlags) ||
              IsTransmissive(bottomFlags));  // otherwise, why bother?

        BxDFFlags flags = BxDFFlags::Reflection;
        if (IsSpecular(topFlags))
            flags = flags | BxDFFlags::Specular;

        if (IsDiffuse(topFlags) || IsDiffuse(bottomFlags) || albedo_diff || albedo_flakes)
            flags = flags | BxDFFlags::Diffuse;
        else if (IsGlossy(topFlags) || IsGlossy(bottomFlags))
            flags = flags | BxDFFlags::Glossy;

        if (IsTransmissive(topFlags) && IsTransmissive(bottomFlags))
            flags = flags | BxDFFlags::Transmission;

        return flags;
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        // The top interface is a fake one, since it is a dielectric with Ior = 1.0 and roughness 0.0
        // The bottom interface is an actual dielectric with ior > 1.0 and roughness > 0.0
        SampledSpectrum f(0.);
        // Next line is for debugging purposes only, D.L.
        //  Estimate _LayeredBxDFMixedParticlesTransmittance_ value _f_ using random sampling
        //  Set _wo_ and _wi_ for layered BSDF evaluation
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }
        // Determine entrance interface for layered BSDF
        TopBxDF enterInterface;

        bool enteredTop = twoSided || wo.z > 0;
        if (enteredTop)
            enterInterface = top;
        else
            enterInterface = bottom;

        TopBxDF exitInterface, nonExitInterface;
        if (SameHemisphere(wo, wi) ^ enteredTop) {
            exitInterface = bottom;
            nonExitInterface = top;
        } else {
            exitInterface = top;
            nonExitInterface = bottom;
        }
        Float exitZ = (SameHemisphere(wo, wi) ^ enteredTop) ? 0 : thickness;

        // Account for reflection at the entrance interface
        if (SameHemisphere(wo, wi))
        {   
            f = nSamples * enterInterface.f(wo, wi, mode);
        }        

        // Declare _RNG_ for layered BSDF evaluation
        RNG rng(Hash(GetOptions().seed, wo), Hash(wi));
        
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        for (int s = 0; s < nSamples; ++s) {
            // Sample random walk through layers to estimate BSDF value
            // Sample transmission direction through entrance interface
            Float uc = r();
            pstd::optional<BSDFSample> wos = enterInterface.Sample_f(
                wo, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Transmission);
            if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0)
                continue;

            // Sample BSDF for virtual light from _wi_
            uc = r();
            pstd::optional<BSDFSample> wis = exitInterface.Sample_f(
                wi, uc, Point2f(r(), r()), !mode, BxDFReflTransFlags::Transmission);
            if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0)
                continue;
         
            // Declare state for random walk through BSDF layers
            SampledSpectrum beta = wos->f * AbsCosTheta(wos->wi) / wos->pdf;
            Float z = enteredTop ? thickness : 0;
            // std::cout<<" z is "<< z << std::endl;
            // Debugging D.L.
            //  wis->wi = -wi;
            //  wos->wi = -wo;
            // std::cout<< " " << AbsCosTheta(wos->wi) << "  " << AbsCosTheta(wo) << std::endl;
            
            Vector3f w = wos->wi;

            // SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g);
            SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g,albedo_flakes,albedo_diff);
            
            // SGGXPhaseFunctionNormalDistribution phase(sigma,mean,std_or);

            // Diffuse phase function
            Point2f u1{r(), r()}, u2{r(), r()};
            bool is_specular_phase = use_specular;
            // SampledSpectrum albedo_particle;

            for (int depth = 0; depth < maxDepth; ++depth) {
                // Sample next event for layered BSDF evaluation random walk
                PBRT_DBG("beta: %f %f %f %f, w: %f %f %f, f: %f %f %f %f\n", beta[0],
                         beta[1], beta[2], beta[3], w.x, w.y, w.z, f[0], f[1], f[2],
                         f[3]);
                // Possibly terminate layered BSDF random walk with Russian roulette
                 if (depth > 128 && beta.MaxComponentValue() < 0.01f) {
                 Float q = std::max<Float>(0, 1 - beta.MaxComponentValue());
                 if (r() < q)
                 break;
                beta /= 1 - q;
                PBRT_DBG("After RR with q = %f, beta: %f %f %f %f\n", q, beta[0],
                 beta[1], beta[2], beta[3]);
                }
                // Account for media between layers and possibly scatter
                // if (!true) {
                    //line after was the usual one but in any case this if-else is not used
                if (!true) {
                    // Advance to next layer boundary and update _beta_ for transmittance
                    z = (z == thickness) ? 0 : thickness;
                    beta *= Tr(thickness, w);
                } else {
                    // Sample medium scattering for layered BSDF evaluation
                    // Take into account the projected Area of particles along direction w
                    Float pa = phase.projectedArea(-w, Point2f(r(), r()));
                    //TODO:
                    //change this put the weighting factor
                    // sigma_t *= pa; //Only for SGGX
                    // sigma_t = diff_conc * sigma_t + (1-diff_conc)*(sigma_t * pa) ;
                    Float reduced_st = diff_conc * sigma_t + (1-diff_conc)*(sigma_t * pa) ;
                    
                    if (reduced_st <= 0) {
                        std::cout << "Projected area is NOT positive" << std::endl;
                        continue;
                    }
                    // std::cout<<" pa "<< pa << std::endl;
                    // Float distSampled = (sigmaT == 0.f) ? Infinity: -std::log(1.f -
                    // rng.Uniform<Float>() ) / sigmaT;
                    Float dz = SampleExponential(r(), reduced_st / std::abs(w.z));

                    Float zp = w.z > 0 ? (z + dz) : (z - dz);
                    // std::cout<<"zp " << zp << std::endl;
                    DCHECK_RARE(1e-5, z == zp);
                    if (z == zp)
                        continue;
                    if (0 < zp && zp < thickness) {
                        // Handle scattering event in layered BSDF medium
                        // Account for scattering through _exitInterface_ using _wis_
                        Float wt = 1;
                        pstd::optional<PhaseFunctionSample> ps;
                        if (!IsSpecular(exitInterface.Flags()))
                                wt = PowerHeuristic(
                                    1, wis->pdf, 1,
                                    phase.PDF(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()) ) );
                            Float m_pdf = phase.PDF(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));

                            pa = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                            auto sigmaTAlongWi = diff_conc  + (1-diff_conc)*( pa) ;
                            
                            SampledSpectrum old_f = f;
                            Float pf_p = phase.p(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));
                            SampledSpectrum w_albedo = phase.albedoW(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));
                            // std::cout<<"w_albedo"<< w_albedo << std::endl;
                            f += beta * w_albedo *  wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            // f += beta * w_albedo * pf_p * wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            // f += beta * pf_p * wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            
                            SampledSpectrum albedo_particle;
                            ps = phase.Sample_p(-w, Point2f(r(), r()),Point2f(r(), r()),r(),Point2f(r(), r()));
                            if (ps->p_flag == _MediumFlags::Diffuser){
                                albedo_particle = albedo_diff;
                                // std::cout<<"albedo_diff "<< albedo_diff << std::endl;
                            }else{
                                albedo_particle = albedo_flakes;
                                // std::cout<<"albedo_flakes "<< albedo_flakes << std::endl;

                            }


                            if (IsNaN(f[0])) {
                                std::cout << " nan found in f()_s_05" << std::endl;
                                std::cout << "oldF " << old_f << " pf_p " << pf_p
                                          << std::endl;
                                std::cout << "beta " << beta << " albedo " << w_albedo
                                          << " wt " << wt << " wis->f " << wis->f
                                          << " wis->pdf" << wis->pdf << std::endl;
                                std::cout << "Tr "
                                          << Tr(zp - exitZ, wis->wi, sigmaTAlongWi)
                                          << std::endl;
                            }
   

                        if (!ps || ps->pdf == 0 || ps->wi.z == 0 || ps->p == 0.0f)
                            continue;
                        // original line in code was the following one --> replaced since ps->p and ps->pdf cancels out as they are the same
                        // beta *= albedo * ps->p / ps->pdf;
                        beta *= albedo_particle;
                        w = ps->wi;
                        z = zp;

                        // Possibly account for scattering through _exitInterface_
                        if (((z < exitZ && w.z > 0) || (z > exitZ && w.z < 0)) &&
                            !IsSpecular(exitInterface.Flags())) {
                            // Account for scattering through _exitInterface_
                            SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                            if (fExit) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                Float wt = PowerHeuristic(1, ps->pdf, 1, exitPDF);
                                // auto sigmaTAlongWi = phase.projectedArea(-ps->wi, Point2f(r(), r()));
                                pa = phase.projectedArea(-ps->wi, Point2f(r(), r()));
                                sigmaTAlongWi = diff_conc * 1 + (1-diff_conc)*(1 * pa) ;

                                f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi) *fExit * wt;

                                // f += beta * Tr(zp - exitZ, ps->wi) * fExit * wt;
                            }
                        }

                        continue;
                    }
                    z = Clamp(zp, 0, thickness);
                }

                // Account for scattering at appropriate interface
                if (z == exitZ) {
                    // Account for reflection at _exitInterface_
                    Float uc = r();
                    pstd::optional<BSDFSample> bs = exitInterface.Sample_f(
                        -w, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0){
                        break;
                    }
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;
                    w = bs->wi;

                } else {
                    // Account for scattering at _nonExitInterface_
                    if (!IsSpecular(nonExitInterface.Flags())) {
                        // Add NEE contribution along presampled _wis_ direction
                        Float wt = 1;
                        if (!IsSpecular(exitInterface.Flags()))
                            wt = PowerHeuristic(1.0f, wis->pdf, 1.0f,
                                                nonExitInterface.PDF(-w, -wis->wi, mode));

                        auto pa = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                        auto sigmaTAlongWi = diff_conc  + (1.0f-diff_conc)*(pa) ;
                        
                        f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                             AbsCosTheta(wis->wi) * wt *
                             Tr(thickness, wis->wi, sigmaTAlongWi) * wis->f / wis->pdf;
                    }
                    // Sample new direction using BSDF at _nonExitInterface_
                    Float uc = r();
                    Point2f u(r(), r());
                    pstd::optional<BSDFSample> bs = nonExitInterface.Sample_f(
                        -w, uc, u, mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        break;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;

                    w = bs->wi;

                    if (!IsSpecular(exitInterface.Flags())) {
                        // Add NEE contribution along direction from BSDF sample
                        SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                        if (fExit) {
                            Float wt = 1;
                            if (!IsSpecular(nonExitInterface.Flags())) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                wt = PowerHeuristic(1, bs->pdf, 1, exitPDF);
                            }
                            auto pa = phase.projectedArea(-bs->wi, Point2f(r(), r()));
                            // auto sigmaTAlongWi = phase.projectedArea(-bs->wi, Point2f(r(), r()));
                            auto sigmaTAlongWi = diff_conc * 1 + (1-diff_conc)*(1 * pa) ;

                            f += beta * Tr(thickness, bs->wi, sigmaTAlongWi) * fExit * wt;
                        }
                    }
                }
            }
        }
        if (f.HasNaNs()) {
            std::cout << " nan found in f()" << std::endl;
        }
        return f / nSamples;
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(
        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // return bottom.Sample_f(wo, uc, u, mode);

        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            return {};
        bool enteredTop = twoSided || wo.z > 0;

        // Set _wo_ for layered BSDF sampling
        bool flipWi = false;
        if (twoSided && wo.z < 0) {
            wo = -wo;
            flipWi = true;
        }

        // Sample BSDF at entrance interface to get initial direction _w_
        pstd::optional<BSDFSample> bs =
            enteredTop ? top.Sample_f(wo, uc, u, mode) : bottom.Sample_f(wo, uc, u, mode);

        if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
            return {};
        if (bs->IsReflection()) {
            // if(!enteredTop){
                // std::cout<<"SHOULD NOT HAPPEN"<<std::endl;
            // }
            if (flipWi)
                bs->wi = -bs->wi;
            bs->pdfIsProportional = true;
            return bs;
        }
        
        Vector3f w = bs->wi;
        // std::cout<<"bs->wi "<< w << " wo "<< wo << std::endl ;
        bool specularPath = bs->IsSpecular();

        // Declare _RNG_ for layered BSDF sampling
        RNG rng(Hash(GetOptions().seed, wo), Hash(uc, u));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Declare common variables for layered BSDF sampling
        SampledSpectrum f = bs->f * AbsCosTheta(bs->wi);
                        if (IsInf(f[0])){
                    std::cout<< "it's INF line 1900"<< std::endl;
                }
        Float pdf = bs->pdf;
        Float z = enteredTop ? thickness : 0;

        SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g,albedo_flakes,albedo_diff);
        bool is_specular_phase = use_specular;
        
        for (int depth = 0; depth < maxDepth; ++depth) {
            // Follow random walk through layers to sample layered BSDF
            // Possibly terminate layered BSDF sampling with Russian Roulette
            Float rrBeta = f.MaxComponentValue() / pdf;
            if (depth > 128 && rrBeta < 0.01f) {
            Float q = std::max<Float>(0, 1 - rrBeta);
            if (r() < q)
            return {};
            pdf *= 1 - q;
            }
            if (w.z == 0)
                return {};
            //if one of the two is selected than it makes sense to use this
            if (albedo_diff || albedo_flakes) {
                // Sample potential scattering event in layered medium
                Float pa = phase.projectedArea(w, Point2f(r(), r()));
                
                Float reduced_st = diff_conc * sigma_t + (1-diff_conc)*(sigma_t * pa) ;

                Float dz = SampleExponential(r(), reduced_st / AbsCosTheta(w));
                Float zp = w.z > 0 ? (z + dz) : (z - dz);
                CHECK_RARE(1e-5, zp == z);
                if (zp == z)
                    return {};
                //media interaction
                if (0 < zp && zp < thickness) {
                    pstd::optional<PhaseFunctionSample> ps;
                    //get the albedo of the sampled particle
                    SampledSpectrum albedo_particle;
                    ps = phase.Sample_p(-w, Point2f(r(), r()),Point2f(r(), r()),r(),Point2f(r(), r())); 
                    if (ps->p_flag == _MediumFlags::Diffuser){
                                
                                albedo_particle = albedo_diff;
                            }else{
                                // std::cout<<"happens"<<std::endl;
                                albedo_particle = albedo_flakes;
                    }
                    // if (!albedo_particle){
                        // LOG_ERROR("Albedo particle is not instantiated");
                    // } 
                    if (!ps || ps->pdf == 0 || ps->wi.z == 0){
                        // LOG_ERROR("early termination is not great",pPixel.x, pPixel.y, sampleIndex);
                        return {};

                    }
                    auto old_f = f;
                    if (ps->p_flag == _MediumFlags::Diffuser){
                        f *= albedo_particle;
                    }else{
                        // f *= albedo_particle;
                        f *= albedo_particle * pa;
                    }

                    
                    specularPath = false;
                    w = ps->wi;
                    z = zp;

                    continue;
                }
                z = Clamp(zp, 0, thickness);
                if (z == 0)
                    DCHECK_LT(w.z, 0);
                else
                    DCHECK_GT(w.z, 0);

            } else {
                // Advance to the other layer interface
                z = (z == thickness) ? 0 : thickness;
                f *= Tr(thickness, w);
            }
            // Initialize _interface_ for current interface surface
#ifdef interface  // That's enough out of you, Windows.
#undef interface
#endif
            BottomBxDF interface;
            
            if ( z == 0){
                interface = bottom;
            }
            else{
                interface = top;
            }

            // Sample interface BSDF to determine new path direction
            Float uc = r();
            Point2f u(r(), r());
            
            // std::cout<<"w "<< w << std::endl ;
            pstd::optional<BSDFSample> bs = interface.Sample_f(-w, uc, u, mode);
            if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0){
                        // LOG_ERROR("early termination is not great");
                return {};

            }
            
            f *= bs->f;
            if (IsInf(f[0])){
                    std::cout<< "it's INF line 1990"<< std::endl;
            }
            pdf *= bs->pdf;
            specularPath &= bs->IsSpecular();
            w = bs->wi;
            
            if (bs->IsTransmission() ) {
                BxDFFlags flags = SameHemisphere(wo, w) ? BxDFFlags::Reflection
                                                        : BxDFFlags::Transmission;
                flags |= specularPath ? BxDFFlags::Specular : BxDFFlags::Glossy;
                if (flipWi)
                    w = -w;
                if (f.HasNaNs()) {
                    std::cout << " nan found in sample()" << std::endl;
                }
                if (IsInf(f[0])){
                    std::cout<< "it's INF"<< std::endl;
                }
                // std::cout<<"exiting f "<< f << " w " << w <<  " depth "<< depth << std::endl; 
                if ( z == thickness){
                    //scatter towards the air
                return BSDFSample(f, w, pdf, flags,1.0f);
                }else{
                    //scatter inside the skin
                return BSDFSample(f, w, pdf, flags,1.3);

                }
                // Old Code
                // return BSDFSample(f, w, pdf, flags, 1.0f, true);
            }
            // Scale _f_ by cosine term after scattering at the interface
            f *= AbsCosTheta(bs->wi);
        }
        // LOG_ERROR("Reaching End");
        return {};
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // Set _wo_ and _wi_ for layered BSDF evaluation
        bool enteredTop = twoSided || wo.z > 0;

        //****************************************************
        //                   Debug code
        //****************************************************
        
        //  if(!enteredTop){
        //     // return top.PDF(wo, wi, mode, sampleFlags);
        //     wo = -wo;
        //     if (!(sampleFlags & BxDFReflTransFlags::Reflection))
        //         return 0;
        //     return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;

        //  }

        //****************************************************
        //                   END Debug code
        //****************************************************

        // return bottom.PDF(wo,wi,mode,sampleFlags);
        
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }
                
        // Declare _RNG_ for layered PDF evaluation
        RNG rng(Hash(GetOptions().seed, wi), Hash(wo));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Update _pdfSum_ for reflection at the entrance layer
        Float pdfSum = 0;
        if (SameHemisphere(wo, wi)) {
            auto reflFlag = BxDFReflTransFlags::Reflection;
            pdfSum += enteredTop ? nSamples * top.PDF(wo, wi, mode, reflFlag)
                                 : nSamples * bottom.PDF(wo, wi, mode, reflFlag);
        }

        for (int s = 0; s < nSamples; ++s) {
            // Evaluate layered BSDF PDF sample
            if (SameHemisphere(wo, wi)) {
                // Evaluate TRT term for PDF estimate
                // Old Code
                // TopOrBottomBxDF<TopBxDF, BottomBxDF> rInterface, tInterface;
                // const TopBxDF *rInterface, *tInterface;
                TopBxDF rInterface, tInterface;
                
                if (enteredTop) {
                    rInterface = bottom;
                    tInterface = top;
                } else {
                    rInterface = top;
                    tInterface = bottom;
                }
                // Sample _tInterface_ to get direction into the layers
                auto trans = BxDFReflTransFlags::Transmission;
                pstd::optional<BSDFSample> wos, wis;
                wos = tInterface.Sample_f(wo, r(), {r(), r()}, mode, trans);
                wis = tInterface.Sample_f(wi, r(), {r(), r()}, !mode, trans);

                // Update _pdfSum_ accounting for TRT scattering events
                if (wos && wos->f && wos->pdf > 0 && wis && wis->f && wis->pdf > 0) {
                    if (!IsNonSpecular(tInterface.Flags()))
                        pdfSum += rInterface.PDF(-wos->wi, -wis->wi, mode);
                    else {
                        // Use multiple importance sampling to estimate PDF product
                        pstd::optional<BSDFSample> rs =
                            rInterface.Sample_f(-wos->wi, r(), {r(), r()}, mode);
                        if (rs && rs->f && rs->pdf > 0) {
                            if (!IsNonSpecular(rInterface.Flags()))
                                pdfSum += tInterface.PDF(-rs->wi, wi, mode);
                            else {
                                // Compute MIS-weighted estimate of Equation
                                // (\ref{eq:pdf-triple-canceled-one})
                                Float rPDF = rInterface.PDF(-wos->wi, -wis->wi, mode);
                                Float wt = PowerHeuristic(1, wis->pdf, 1, rPDF);
                                pdfSum += wt * rPDF;

                                Float tPDF = tInterface.PDF(-rs->wi, wi, mode);
                                wt = PowerHeuristic(1, rs->pdf, 1, tPDF);
                                pdfSum += wt * tPDF;
                            }
                        }
                    }
                }

            } else {
                
                // Evaluate TT term for PDF estimate
                // TopOrBottomBxDF<TopBxDF, BottomBxDF> toInterface, tiInterface;
                TopBxDF toInterface, tiInterface;
                
                if (enteredTop) {
                    toInterface = top;
                    tiInterface = bottom;
                } else {
                    toInterface = bottom;
                    tiInterface = top;
                }

                Float uc = r();
                Point2f u(r(), r());
                
                pstd::optional<BSDFSample> wos = toInterface.Sample_f(wo, uc, u, mode);
                if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0 ||
                    wos->IsReflection())
                    continue;

                uc = r();
                u = Point2f(r(), r());
                pstd::optional<BSDFSample> wis = tiInterface.Sample_f(wi, uc, u, !mode);
                if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0 ||
                    wis->IsReflection())
                    continue;

                if (IsSpecular(toInterface.Flags()))
                    pdfSum += tiInterface.PDF(-wos->wi, wi, mode);
                else if (IsSpecular(tiInterface.Flags()))
                    pdfSum += toInterface.PDF(wo, -wis->wi, mode);
                else
                    pdfSum += (toInterface.PDF(wo, -wis->wi, mode) +
                               tiInterface.PDF(-wos->wi, wi, mode)) /
                              2;
            }
        }
        // Return mixture of PDF estimate and constant PDF
        return Lerp(0.9f, 1 / (4 * Pi), pdfSum / nSamples);
    }

  private:
    // LayeredBxDFMixedParticlesTransmittance Private Methods
    PBRT_CPU_GPU
    static Float Tr(Float dz, Vector3f w) {
        if (std::abs(dz) <= std::numeric_limits<Float>::min())
            return 1;
        return FastExp(-std::abs(dz / w.z));
    }
    static Float Tr(Float dz, Vector3f w, Float sigmaTAlongWi) {
        return FastExp(-std::abs(sigmaTAlongWi * dz / w.z));
    }

    // LayeredBxDFMixedParticlesTransmittance Private Members
    TopBxDF top;
    BottomBxDF bottom;
    Float thickness, g1,g2,w_g, sigma, mean, std_or, std_size, diff_conc;
    Float sigma_t,eta;

    Vector3f dir;
    
    SampledSpectrum albedo_diff,albedo_flakes;
    int maxDepth, nSamples;
    bool use_specular;
};


template <typename TopBxDF, typename BottomBxDF, bool twoSided>
class LayeredBxDFMixedParticlesNormalizedFresnelAdvanced {
  public:
    // LayeredBxDFMixedParticlesNormalizedFresnelAdvanced Public Methods
    LayeredBxDFMixedParticlesNormalizedFresnelAdvanced() = default;
    PBRT_CPU_GPU
    LayeredBxDFMixedParticlesNormalizedFresnelAdvanced(TopBxDF top, BottomBxDF bottom, CosmeticStructBxdf media1,CosmeticStructBxdf media2,  int maxDepth,int nSamples)
        : top(top),
          bottom(bottom),
          media1(media1),
          media2(media2),
          maxDepth(maxDepth),
          nSamples(nSamples),
          eta(1.3f) {
            
            medias[0] = media1;
            medias[1] = media2;

        // std::cout<<"std SIZE is " << std_size << std::endl;
    }



    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() {
        top.Regularize();
        bottom.Regularize();
    }

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        //TODO - remove this        
        // SampledSpectrum albedo(0.0f);
        BxDFFlags topFlags = top.Flags(), bottomFlags = bottom.Flags();
        CHECK(IsTransmissive(topFlags) ||
              IsTransmissive(bottomFlags));  // otherwise, why bother?

        BxDFFlags flags = BxDFFlags::Reflection;
        if (IsSpecular(topFlags))
            flags = flags | BxDFFlags::Specular;

        if (IsDiffuse(topFlags) || IsDiffuse(bottomFlags) || media1.albedo_diff || media1.albedo_flakes)
            flags = flags | BxDFFlags::Diffuse;
        else if (IsGlossy(topFlags) || IsGlossy(bottomFlags))
            flags = flags | BxDFFlags::Glossy;

        if (IsTransmissive(topFlags) && IsTransmissive(bottomFlags))
            flags = flags | BxDFFlags::Transmission;

        return flags;
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        //TODO: implement using hemisphere sampling (check the NormalizedFresnelBxDF class)
        // The top interface is a fake one, since it is a dielectric with Ior = 1.0 and roughness 0.0
        // The bottom interface is an actual dielectric with ior > 1.0 and roughness > 0.0
        SampledSpectrum f(0.0);

//        return f;
        return bottom.f(wo, wi, mode);

        // Next line is for debugging purposes only, D.L.
        //  Estimate _LayeredBxDFMixedParticlesTransmittance_ value _f_ using random sampling
        //  Set _wo_ and _wi_ for layered BSDF evaluation
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }
        // Determine entrance interface for layered BSDF
        // TopBxDF enterInterface;
        TopOrBottomBxDF<TopBxDF, BottomBxDF> enterInterface;
        bool enteredTop = false;
        CosmeticStructBxdf selected_media;
        enterInterface = &bottom;
        int exit_media_index;
        int last_index_media = 1;
        TopOrBottomBxDF<TopBxDF, BottomBxDF> exitInterface, nonExitInterface;
        // TopBxDF exitInterface, nonExitInterface;

        int media_index = 0;
        if ((medias[0].thickness < MachineEpsilon && medias[0].thickness > -MachineEpsilon) && (medias[1].thickness < MachineEpsilon && medias[1].thickness > -MachineEpsilon)){
            //just get the value of the bottomLayer material
            return bottom.f(wo, wi, mode);
        }
    
        // Float exitZ = (SameHemisphere(wo, wi) ^ enteredTop) ? 0 : thickness;
        RNG rng(Hash(GetOptions().seed, wo), Hash(wi));
        
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };
        // Account for reflection at the entrance interface
        /*
        if (SameHemisphere(wo, wi))
        {   
            f = nSamples * enterInterface.f(wo, wi, mode);
            selected_media = medias[0];
            if ( selected_media.thickness!=0){
            SGGXPhaseFunctionMixedRefl phase_bot(selected_media.platelets_roughness, selected_media.platelets_orientation, selected_media.std_or, selected_media.diff_conc,selected_media.g1,selected_media.g2,selected_media.w_g,selected_media.albedo_flakes,selected_media.albedo_diff);
            auto pa_bot = phase_bot.projectedArea(-wi, Point2f(r(), r()));
            auto sigmaTAlongWi_bot = medias[0].diff_conc  + (1-medias[0].diff_conc)*( pa_bot) ;
            f *=   Tr(selected_media.thickness, wi, sigmaTAlongWi_bot);
            }
            selected_media = medias[1];
            
            if ( selected_media.thickness!=0){
            SGGXPhaseFunctionMixedRefl phase_top(selected_media.platelets_roughness, selected_media.platelets_orientation, selected_media.std_or, selected_media.diff_conc,selected_media.g1,selected_media.g2,selected_media.w_g,selected_media.albedo_flakes,selected_media.albedo_diff);
            auto pa_top = phase_top.projectedArea(-wi, Point2f(r(), r()));
            auto sigmaTAlongWi_top = medias[1].diff_conc  + (1-selected_media.diff_conc)*( pa_top) ;
            f *=   Tr(selected_media.thickness, wi, sigmaTAlongWi_top);
            }
        }*/
        //return f/nSamples;
        
        selected_media = medias[0];
        Float thickness = selected_media.thickness;
        if (!SameHemisphere(wo, wi) ) {
            exitInterface = &bottom;
            nonExitInterface = &top;
            exit_media_index = 0;
            return {};            
        } else {
            exitInterface = &top;
            nonExitInterface = &bottom;
            exit_media_index = 1;
        }
 
        //Float exitZ =  selected_media.thickness;
        Float exitZ = (SameHemisphere(wo, wi) ^ enteredTop) ? 0 : thickness;

        // Declare _RNG_ for layered BSDF evaluation


        for (int s = 0; s < nSamples; ++s) {
            if (SameHemisphere(wo, wi))
            f +=  enterInterface.f(wo, wi, mode);

            // Sample random walk through layers to estimate BSDF value
            // Sample transmission direction through entrance interface
            Float uc = r();
            pstd::optional<BSDFSample> wos = enterInterface.Sample_f(
                wo, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Reflection);
            if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0)
                continue;

            // Sample BSDF for virtual light from _wi_
            uc = r();
            pstd::optional<BSDFSample> wis = exitInterface.Sample_f(
                wi, uc, Point2f(r(), r()), !mode, BxDFReflTransFlags::Transmission);
            if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0)
                continue;
         
            // Declare state for random walk through BSDF layers
            SampledSpectrum beta = wos->f * AbsCosTheta(wos->wi) / wos->pdf;
            Float z = enteredTop ? thickness : 0;            
            Vector3f w = wos->wi;
            // SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g);
            
            // SGGXPhaseFunctionNormalDistribution phase(sigma,mean,std_or);

            // Diffuse phase function
            Point2f u1{r(), r()}, u2{r(), r()};
            // SampledSpectrum albedo_particle;

            for (int depth = 0; depth < maxDepth; ++depth) {
                selected_media = medias[media_index];
                SGGXPhaseFunctionMixedRefl phase(selected_media.platelets_roughness, selected_media.platelets_orientation, selected_media.std_or, selected_media.diff_conc,selected_media.g1,selected_media.g2,selected_media.w_g,selected_media.albedo_flakes,selected_media.albedo_diff);
                thickness = selected_media.thickness;
//                if (thickness > 0)
//                std::cout<<"thickness "<< thickness << std::endl;
                // Sample next event for layered BSDF evaluation random walk
                PBRT_DBG("beta: %f %f %f %f, w: %f %f %f, f: %f %f %f %f\n", beta[0],
                         beta[1], beta[2], beta[3], w.x, w.y, w.z, f[0], f[1], f[2],
                         f[3]);
                // Possibly terminate layered BSDF random walk with Russian roulette
                if (depth > 128 && beta.MaxComponentValue() < 0.01f) {
                 Float q = std::max<Float>(0, 1 - beta.MaxComponentValue());
                 if (r() < q)
                 break;
                beta /= 1 - q;
                PBRT_DBG("After RR with q = %f, beta: %f %f %f %f\n", q, beta[0],
                 beta[1], beta[2], beta[3]);
                }

                    // Sample medium scattering for layered BSDF evaluation
                    // Take into account the projected Area of particles along direction w
                    Float pa = phase.projectedArea(-w, Point2f(r(), r()));
                    Float reduced_st = selected_media.diff_conc * selected_media.sigma_t + (1-selected_media.diff_conc)*(selected_media.sigma_t * pa) ;
                    
                    if (reduced_st <= 0) {
                        std::cout << "Projected area is NOT positive" << std::endl;
                        continue;
                    }
                    // std::cout<<" pa "<< pa << std::endl;
                    // Float distSampled = (sigmaT == 0.f) ? Infinity: -std::log(1.f -
                    // rng.Uniform<Float>() ) / sigmaT;
                    Float dz = SampleExponential(r(), reduced_st / std::abs(w.z));

                    Float zp = w.z > 0 ? (z + dz) : (z - dz);
                    // std::cout<<"zp " << zp << std::endl;
                    DCHECK_RARE(1e-5, z == zp);
                    if (z == zp)
                        continue;
                    if (0 < zp && zp < thickness && thickness > 0 ) {
                        // Handle scattering event in layered BSDF medium
                        // Account for scattering through _exitInterface_ using _wis_
                        
                        Float wt = 1;
                        pstd::optional<PhaseFunctionSample> ps;
                        if (!IsSpecular(exitInterface.Flags()))
                                wt = PowerHeuristic(1, wis->pdf, 1,phase.PDF(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()) ) );
                        Float m_pdf = phase.PDF(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));

                            pa = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                            
                            SampledSpectrum old_f = f;
                            Float pf_p = phase.p(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));
                            SampledSpectrum w_albedo = phase.albedoW(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));
                            // std::cout<<"w_albedo"<< w_albedo << std::endl;
                            auto sigmaTAlongWi = selected_media.diff_conc  + (1-selected_media.diff_conc)*( pa) ;
                            if (media_index == exit_media_index){
                                f += beta * w_albedo *  wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            }else{
                                auto sigmaTAlongWi_exit_media = medias[exit_media_index].diff_conc  + (1-medias[exit_media_index].diff_conc)*( pa) ;
                                //f += beta * w_albedo *  wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                                
                                f += beta * w_albedo *  wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * Tr(medias[exit_media_index].thickness, wis->wi, sigmaTAlongWi_exit_media) * wis->f /wis->pdf;
                            }
                            // f += beta * w_albedo * pf_p * wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            // f += beta * pf_p * wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            
                            SampledSpectrum albedo_particle;
                            ps = phase.Sample_p(-w, Point2f(r(), r()),Point2f(r(), r()),r(),Point2f(r(), r()));
                            if (ps->p_flag == _MediumFlags::Diffuser){
                                albedo_particle = selected_media.albedo_diff;
                                // std::cout<<"albedo_diff "<< albedo_diff << std::endl;
                            }else{
                                albedo_particle = selected_media.albedo_flakes;
                                // std::cout<<"albedo_flakes "<< albedo_flakes << std::endl;

                            }


                            if (IsNaN(f[0])) {
                                std::cout << " l3170" << std::endl;
                                
                                std::cout << " nan found in f()_s_05" << std::endl;
                                std::cout << "oldF " << old_f << " pf_p " << pf_p
                                          << std::endl;
                                std::cout << "beta " << beta << " albedo " << w_albedo
                                          << " wt " << wt << " wis->f " << wis->f
                                          << " wis->pdf" << wis->pdf << std::endl;
                                std::cout << "Tr "
                                          << Tr(zp - exitZ, wis->wi, sigmaTAlongWi)
                                          << std::endl;
                            }
   

                        if (!ps || ps->pdf == 0 || ps->wi.z == 0 || ps->p == 0.0f)
                            continue;
                        // original line in code was the following one --> replaced since ps->p and ps->pdf cancels out as they are the same
                        // beta *= albedo * ps->p / ps->pdf;
                        beta *= albedo_particle;
                        w = ps->wi;
                        z = zp;

                        // Possibly account for scattering through _exitInterface_
                        if (((z < exitZ && w.z > 0 ) || (z > exitZ && w.z < 0 )) ){                        
                       // if (((z < exitZ && w.z > 0 ) || (z > exitZ && w.z < 0 )) && !IsSpecular(exitInterface.Flags())) {
                            // Account for scattering through _exitInterface_
                            SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                            if (true) {
//                                std::cout<<"herz l3222"<<std::endl;
                                Float exitPDF = exitInterface.PDF(-w, wi, mode, BxDFReflTransFlags::Transmission);
                                Float wt = PowerHeuristic(1, ps->pdf, 1, exitPDF);
                                // auto sigmaTAlongWi = phase.projectedArea(-ps->wi, Point2f(r(), r()));
                                pa = phase.projectedArea(-ps->wi, Point2f(r(), r()));
                                if (w.z > 0){
                                    //in case the ray is pointing upward
                                    sigmaTAlongWi = selected_media.diff_conc + (1-selected_media.diff_conc)*(1 * pa) ;
                                    if (media_index < last_index_media && medias[media_index+1].thickness > 0){
                                        //in case we are in the bottom layer, we have to account for the top layer
                                        auto sigmaTAlongWi_top_media = medias[media_index+1].diff_conc + (1-medias[media_index+1].diff_conc)*(1 * pa) ;
                                        f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi)*Tr( medias[media_index+1].thickness, ps->wi, sigmaTAlongWi_top_media)  * wt;
                                    }else{
                                        f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi)  * wt;
                                    }
                                }
                                if (w.z < 0){
                                    //in case the ray is pointing upward
                                    sigmaTAlongWi = selected_media.diff_conc + (1-selected_media.diff_conc)*(1 * pa) ;
                                    if (media_index > 0 && medias[media_index-1].thickness>0){
                                        //in case we are in the top layer, we have to account for the bottom layer
                                        auto sigmaTAlongWi_bot_media = medias[media_index-1].diff_conc + (1-medias[media_index-1].diff_conc)*(1 * pa) ;

                                        f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi)*Tr(0 - medias[media_index-1].thickness, ps->wi, sigmaTAlongWi_bot_media)  * wt;
                                    }else{
                                        f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi)  * wt;
                                    }
                                }
                                
                                // f += beta * Tr(zp - exitZ, ps->wi) * fExit * wt;
                            }
                        }

                        continue;
                    }

                    z = Clamp(zp, 0, thickness);
                    if ( z == 0){
                        media_index--;
                        if (media_index>-1){
                            selected_media = medias[media_index];
                            z = selected_media.thickness;
                            if (exit_media_index == 0){
                                exitZ = 0;
                            }else{
                                exitZ = selected_media.thickness;
                            }
                            continue;
                        }
                    }else{
                        media_index++;
                        if (media_index<2){
                            selected_media = medias[media_index];
                            z = 0;
                            if (exit_media_index == 0){
                                exitZ = 0;
                            }else{
                                exitZ = selected_media.thickness;
                            }
                            continue;
                        }

                    }
                
                if (media_index == 2 || media_index == -1){
                    media_index = (media_index == 2) ? 1 : 0;
                }else{
                    std::cout<<"media index " << media_index << " thickness"<<thickness<<std::endl;
                }
                
                //TODO fix the next part
                // Account for scattering at appropriate interface
                if (z == exitZ  ) {
                    if (media_index != exit_media_index ){
                    //    if (!medias[exit_media_index].thickness == 0)
//                            std::cout<<"Media index << "<<  media_index << "is different than exit_media_index " << exit_media_index <<std::endl;
                        continue;
                        
                    }
                    Float uc = r();
                    pstd::optional<BSDFSample> bs = exitInterface.Sample_f(
                        -w, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Reflection);
                    
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0){
                        //no reflection is possible, then let's start a new iteration
                        break;
                    }
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;
                    w = bs->wi;

                } else {
                    // Account for scattering at _nonExitInterface_
                    if ( !IsSpecular(nonExitInterface.Flags()) ) {

                        //std::cout<<"IsNotSpecular"<<std::endl;
                        // Add NEE contribution along presampled _wis_ direction

                        Float wt = 1;
                        if (media_index == exit_media_index){
                            //std::cout<<"Media index == exit_media_index. ERROR?!?"<<std::endl;
                        continue;
                        }
//                        std::cout<<"reached l3324"<<std::endl;

                        if (!IsSpecular(exitInterface.Flags()))
                            wt = PowerHeuristic(1.0f, wis->pdf, 1.0f,
                                                nonExitInterface.PDF(-w, -wis->wi, mode));

                        auto pa = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                        
                        auto sigmaTAlongWi = medias[0].diff_conc  + (1.0f-medias[0].diff_conc)*(pa) ;
                        //auto sigmaTAlongWi = selected_media.diff_conc  + (1.0f-selected_media.diff_conc)*(pa) ;
                        auto sigmaTAlongWi_exit_media = medias[exit_media_index].diff_conc  + (1.0f-medias[exit_media_index].diff_conc)*(pa) ;

                        // f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                        //      AbsCosTheta(wis->wi) * wt *
                        //      Tr(thickness, wis->wi, sigmaTAlongWi) * wis->f / wis->pdf;
                        
                        f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                             AbsCosTheta(wis->wi) * wt *
                             Tr(medias[0].thickness, wis->wi, sigmaTAlongWi)* Tr(medias[1].thickness, wis->wi, sigmaTAlongWi_exit_media) * wis->f / wis->pdf;
                    }
                    // Sample new direction using BSDF at _nonExitInterface_
                    Float uc = r();
                    Point2f u(r(), r());
                    pstd::optional<BSDFSample> bs = nonExitInterface.Sample_f(
                        -w, uc, u, mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        break;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;

                    w = bs->wi;

                    /*if (!IsSpecular(exitInterface.Flags())) {
                        // Add NEE contribution along direction from BSDF sample
                        SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                        if (fExit) {
                            Float wt = 1;
                            if (!IsSpecular(nonExitInterface.Flags())) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                wt = PowerHeuristic(1, bs->pdf, 1, exitPDF);
                            }
                            auto pa = phase.projectedArea(-bs->wi, Point2f(r(), r()));
                            // auto sigmaTAlongWi = phase.projectedArea(-bs->wi, Point2f(r(), r()));
                            //TODO fix this
                            auto sigmaTAlongWi = selected_media.diff_conc * 1 + (1-selected_media.diff_conc)*(1 * pa) ;

                            f += beta * Tr(thickness, bs->wi, sigmaTAlongWi) * fExit * wt;
                        }
                    }*/
                }
            }
        }
        if (f.HasNaNs()) {
            std::cout << " nan found in f()" << std::endl;
        }
        return f / nSamples;
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(

      

        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now

        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            return {};
        
        // Set _wo_ for layered BSDF sampling
        bool flipWi = false;
        if (twoSided && wo.z < 0) {
            wo = -wo;
            flipWi = true;
        }
        bool enteredTop = false;
        // pstd::optional<BSDFSample> bs =            enteredTop ? top.Sample_f(wo, uc, u, mode) : bottom.Sample_f(wo, uc, u, mode);
        pstd::optional<BSDFSample> bs = bottom.Sample_f(wo, uc, u, mode,sampleFlags);
        return bs;
        // We don't need to select one of the two materials, since this material is used only for BSSRDF scattering, 
        // hence we will always sample the bot interface 
        //

        //this vector will always be toward the upward direction
        
//        Vector3f wi = SampleCosineHemisphere(u);
//        SampledSpectrum f = bottom.f(wo, wi, mode) * AbsCosTheta(wi);
//        Float pdf = AbsCosTheta(wi) * InvPi;

//        if (wo.z < 0)
//            wi.z *= -1;
        SampledSpectrum f = bs->f;
        Float pdf = bs->pdf;
        Vector3f w  = bs->wi;
        BxDFFlags m_flags = SameHemisphere(wo, w) ? BxDFFlags::Reflection: BxDFFlags::Transmission;

        //return BSDFSample(SampledSpectrum(0.0f), w, pdf,m_flags,1.3);
        
        if (medias[0].thickness == 0 && medias[1].thickness == 0){
            BxDFFlags flags = SameHemisphere(wo, w) ? BxDFFlags::Reflection: BxDFFlags::Transmission;
            return BSDFSample(f, w, pdf, flags,1.3);
            //return bs;
        }

//        pstd::optional<BSDFSample> bs;
//        return BSDFSample(f(wo, wi, mode), wi, PDF(wo, wi, mode, sampleFlags),BxDFFlags::DiffuseReflection);
        //return bs;
        //special case for when we apply no cosmetic at all.
        
        //if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
        //    return {};
        
        /*if (bs->IsReflection()) {
                    if (flipWi)
                bs->wi = -bs->wi;
            bs->pdfIsProportional = true;
            return bs;
        }*/
        
        //Vector3f w = bs->wi;
        // std::cout<<"bs->wi "<< w << " wo "<< wo << std::endl ;
        bool specularPath = false;//bs->IsSpecular();

        // Declare _RNG_ for layered BSDF sampling
        RNG rng(Hash(GetOptions().seed, wo), Hash(uc, u));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };
        
        // size_t media_index = 0;
        size_t last_index = 1;
               
        // Declare common variables for layered BSDF sampling
        
        if (IsInf(f[0])){
                    std::cout<< "it's INF line 1900"<< std::endl;
                }
        size_t media_index = enteredTop ? 1 : 0;

        Float thickness = medias[media_index].thickness;
        Float z = enteredTop ? thickness : 0;
        // std::cout<<"media1.diffprec" << medias[0].diff_conc << "media2.diffprce " << medias[1].diff_conc<< std::endl; 
        // phase.toString();
        // SGGXPhaseFunctionMixedRefl phase(media1.sigma, media1.mean, media1.std_or, media1.diff_conc,media1.g1,media1.g2,media1.w_g,media1.albedo_flakes,media1.albedo_diff);
        // std::cout<<"maxDepth "<< maxDepth << " thick "<< thickness << " sigmat " << sigma_t << std::endl;
        
        for (int depth = 0; depth < maxDepth; ++depth) {
            // std::cout<<" d "<< depth << std::endl;
            // Follow random walk through layers to sample layered BSDF
            // Possibly terminate layered BSDF sampling with Russian Roulette
            Float rrBeta = f.MaxComponentValue() / pdf;
            if (depth > 128 && rrBeta < 0.01f) {
            Float q = std::max<Float>(0, 1 - rrBeta);
            if (r() < q)
            return {};
            pdf *= 1 - q;
            }
            if (w.z == 0)
                return {};
            //if one of the two is selected than it makes sense to use this
            
            CosmeticStructBxdf selected_media = medias[media_index];
            thickness = selected_media.thickness;
            
            if (thickness == 0.0){
                if (w.z > 0){
                media_index++;
                if(media_index < 2){
                z = 0;
                continue;
                }else{
                    return bs;
                }
                //this media has no thickness go to the next one
                }
            }
            Float sigma_t = selected_media.sigma_t;

            SGGXPhaseFunctionMixedRefl phase(selected_media.platelets_roughness, selected_media.platelets_orientation, selected_media.std_or, selected_media.diff_conc,selected_media.g1,selected_media.g2,selected_media.w_g,selected_media.albedo_flakes,selected_media.albedo_diff);
                   
            if (selected_media.albedo_diff || selected_media.albedo_flakes) {
                // Sample potential scattering event in layered medium
                Float pa = phase.projectedArea(w, Point2f(r(), r()));
                
                Float reduced_st = selected_media.diff_conc * sigma_t + (1-selected_media.diff_conc)*(sigma_t * pa) ;

                Float dz = SampleExponential(r(), reduced_st / AbsCosTheta(w));
                Float zp = w.z > 0 ? (z + dz) : (z - dz);
                CHECK_RARE(1e-5, zp == z);
                if (zp == z)
                    return {};
                //media interaction
                if (0 < zp && zp < thickness) {
                    pstd::optional<PhaseFunctionSample> ps;
                    //get the albedo of the sampled particle
                    SampledSpectrum albedo_particle;
                    ps = phase.Sample_p(-w, Point2f(r(), r()),Point2f(r(), r()),r(),Point2f(r(), r())); 
                    if (ps->p_flag == _MediumFlags::Diffuser){
                                
                                albedo_particle = selected_media.albedo_diff;
                            }else{
                                // std::cout<<"happens"<<std::endl;
                                albedo_particle = selected_media.albedo_flakes;
                    }
                    // if (!albedo_particle){
                        // LOG_ERROR("Albedo particle is not instantiated");
                    // } 
                    if (!ps || ps->pdf == 0 || ps->wi.z == 0){
                        // LOG_ERROR("early termination is not great",pPixel.x, pPixel.y, sampleIndex);
                        return {};

                    }
                    auto old_f = f;
                    
                    if (ps->p_flag == _MediumFlags::Diffuser){
                        f *= albedo_particle;
                    }else{
                        // f *= albedo_particle;
                        
                        f *= albedo_particle * pa;
                    
                    }
                    if (IsInf(f[0])){
                    std::cout<< "it's INF"<< std::endl;
                    std::cout<<"albedo particle" << albedo_particle << " pa "<< pa << media_index << std::endl;
            }
                    
                    specularPath = false;
                    w = ps->wi;
                    z = zp;
                    continue;
                }
                z = Clamp(zp, 0, thickness);
                
                if (z == 0){
                    DCHECK_LT(w.z, 0);
                    media_index --;
                    //going down to 1 media
                    if (media_index > -1){
                    thickness = medias[media_index].thickness;
                    z = thickness;
                        
                    }

                }
                else{
                    media_index++;
                    DCHECK_GT(w.z, 0);
                    if (media_index < last_index+1){
                    //moving upward
                    thickness = medias[media_index].thickness;
                    z = 0;
                        
                    }

                }
                    if (IsInf(f[0])){
                    std::cout<< "it's INF line 3478"<< std::endl;
            }
                } else {
                    // Advance to the other layer interface
                    z = (z == thickness) ? 0 : thickness;
                    f *= Tr(thickness, w);
                }
            if (media_index > -1 || media_index < last_index+1){
                //sample is still inside one of the medias
                continue;
            }else{
                //clip to legal values
                // std::cout<<"media index "<< media_index<< std::endl;
                media_index = (media_index == -1) ? 0 : last_index;
            }
            // Initialize _interface_ for current interface surface
#ifdef interface  // That's enough out of you, Windows.
#undef interface
#endif
            TopOrBottomBxDF<TopBxDF, BottomBxDF> interface;
            //BottomBxDF interface;
            bool isBot = false;
            if ( z == 0){
                isBot = true;
                interface = &bottom;
            }
            else{
                interface = &top;
            }

            // Sample interface BSDF to determine new path direction
            Float uc = r();
            Point2f u(r(), r());
            
            // std::cout<<"w "<< w << std::endl ;
            pstd::optional<BSDFSample> bs = interface.Sample_f(-w, uc, u, mode);
            if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0){
                // std::cout<<"-w " << -w << " uc "<< uc << " u " << u << " mode " << mode << std::endl;
                return {};
            }

            f *= bs->f;

            pdf *= bs->pdf;
            specularPath &= bs->IsSpecular();
            w = bs->wi;
            
            if (bs->IsTransmission() ) {
                //BxDFFlags flags = SameHemisphere(wo, w) ? BxDFFlags::Reflection : BxDFFlags::Transmission;
                BxDFFlags flags =  BxDFFlags::Transmission;
                flags |= specularPath ? BxDFFlags::Specular : BxDFFlags::Glossy;
                if (flipWi)
                    w = -w;
                if (f.HasNaNs()) {
                    std::cout << " nan found in sample()" << std::endl;
                }
                if (IsInf(f[0])){
                    std::cout<< "it's INF"<< std::endl;
                }
                // std::cout<<"exiting f "<< f << " w " << w <<  " depth "<< depth << std::endl; 
                if ( z == thickness){
                    //scatter towards the air
                return BSDFSample(f, w, pdf, flags,1.0f);
                }else{
                    //scatter inside the skin
                return BSDFSample(f, w, pdf, flags,1.3);

                }
                // Old Code
                // return BSDFSample(f, w, pdf, flags, 1.3f, true);
            }
            // Scale _f_ by cosine term after scattering at the interface
            f *= AbsCosTheta(bs->wi);
        }
        //TODO check why this is happening
        LOG_ERROR("End Reached");
        return {};
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // Set _wo_ and _wi_ for layered BSDF evaluation
        bool enteredTop = twoSided || wo.z > 0;

        //****************************************************
        //                   Debug code
        //****************************************************
        
        //  if(!enteredTop){
        //     // return top.PDF(wo, wi, mode, sampleFlags);
        //     wo = -wo;
        //     if (!(sampleFlags & BxDFReflTransFlags::Reflection))
        //         return 0;
        //     return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;

        //  }

        //****************************************************
        //                   END Debug code
        //****************************************************

        // return bottom.PDF(wo,wi,mode,sampleFlags);
        
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }
                
        // Declare _RNG_ for layered PDF evaluation
        RNG rng(Hash(GetOptions().seed, wi), Hash(wo));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Update _pdfSum_ for reflection at the entrance layer
        Float pdfSum = 0;
        if (SameHemisphere(wo, wi)) {
            auto reflFlag = BxDFReflTransFlags::Reflection;
            pdfSum += enteredTop ? nSamples * top.PDF(wo, wi, mode, reflFlag)
                                 : nSamples * bottom.PDF(wo, wi, mode, reflFlag);
        }

        for (int s = 0; s < nSamples; ++s) {
            // Evaluate layered BSDF PDF sample
            if (SameHemisphere(wo, wi)) {
                // Evaluate TRT term for PDF estimate
                // Old Code
                TopOrBottomBxDF<TopBxDF, BottomBxDF> rInterface, tInterface;
                // const TopBxDF *rInterface, *tInterface;
                // TopBxDF rInterface, tInterface;
                
                if (enteredTop) {
                    rInterface = &bottom;
                    tInterface = &top;
                } else {
                    rInterface = &top;
                    tInterface = &bottom;
                }
                // Sample _tInterface_ to get direction into the layers
                auto trans = BxDFReflTransFlags::Transmission;
                pstd::optional<BSDFSample> wos, wis;
                wos = tInterface.Sample_f(wo, r(), {r(), r()}, mode, trans);
                wis = tInterface.Sample_f(wi, r(), {r(), r()}, !mode, trans);

                // Update _pdfSum_ accounting for TRT scattering events
                if (wos && wos->f && wos->pdf > 0 && wis && wis->f && wis->pdf > 0) {
                    if (!IsNonSpecular(tInterface.Flags()))
                        pdfSum += rInterface.PDF(-wos->wi, -wis->wi, mode);
                    else {
                        // Use multiple importance sampling to estimate PDF product
                        pstd::optional<BSDFSample> rs =
                            rInterface.Sample_f(-wos->wi, r(), {r(), r()}, mode);
                        if (rs && rs->f && rs->pdf > 0) {
                            if (!IsNonSpecular(rInterface.Flags()))
                                pdfSum += tInterface.PDF(-rs->wi, wi, mode);
                            else {
                                // Compute MIS-weighted estimate of Equation
                                // (\ref{eq:pdf-triple-canceled-one})
                                Float rPDF = rInterface.PDF(-wos->wi, -wis->wi, mode);
                                Float wt = PowerHeuristic(1, wis->pdf, 1, rPDF);
                                pdfSum += wt * rPDF;

                                Float tPDF = tInterface.PDF(-rs->wi, wi, mode);
                                wt = PowerHeuristic(1, rs->pdf, 1, tPDF);
                                pdfSum += wt * tPDF;
                            }
                        }
                    }
                }

            } else {
                
                // Evaluate TT term for PDF estimate
                TopOrBottomBxDF<TopBxDF, BottomBxDF> toInterface, tiInterface;
                // TopBxDF toInterface, tiInterface;
                
                if (enteredTop) {
                    toInterface = &top;
                    tiInterface = &bottom;
                } else {
                    toInterface = &bottom;
                    tiInterface = &top;
                }

                Float uc = r();
                Point2f u(r(), r());
                
                pstd::optional<BSDFSample> wos = toInterface.Sample_f(wo, uc, u, mode);
                if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0 ||
                    wos->IsReflection())
                    continue;

                uc = r();
                u = Point2f(r(), r());
                pstd::optional<BSDFSample> wis = tiInterface.Sample_f(wi, uc, u, !mode);
                if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0 ||
                    wis->IsReflection())
                    continue;

                if (IsSpecular(toInterface.Flags()))
                    pdfSum += tiInterface.PDF(-wos->wi, wi, mode);
                else if (IsSpecular(tiInterface.Flags()))
                    pdfSum += toInterface.PDF(wo, -wis->wi, mode);
                else
                    pdfSum += (toInterface.PDF(wo, -wis->wi, mode) +
                               tiInterface.PDF(-wos->wi, wi, mode)) /
                              2;
            }
        }
        // Return mixture of PDF estimate and constant PDF
        return Lerp(0.9f, 1 / (4 * Pi), pdfSum / nSamples);
    }

  private:
    // LayeredBxDFMixedParticlesNormalizedFresnelAdvanced Private Methods
    PBRT_CPU_GPU
    static Float Tr(Float dz, Vector3f w) {
        if (std::abs(dz) <= std::numeric_limits<Float>::min())
            return 1;
        return FastExp(-std::abs(dz / w.z));
    }
    static Float Tr(Float dz, Vector3f w, Float sigmaTAlongWi) {
        return FastExp(-std::abs(sigmaTAlongWi * dz / w.z));
    }

    // LayeredBxDFMixedParticlesNormalizedFresnelAdvanced Private Members
    TopBxDF top;
    BottomBxDF bottom;

    // Float thickness, g1,g2,w_g, sigma, mean, std_or, std_size, diff_conc;
    // Float sigma_t,eta;
    // Vector3f dir;
    // SampledSpectrum albedo_diff,albedo_flakes;
    Float eta;
    CosmeticStructBxdf media1,media2;
    CosmeticStructBxdf medias[2];
    int maxDepth, nSamples;
    bool use_specular;
};



template <typename TopBxDF, typename BottomBxDF, bool twoSided>
class LayeredBxDFMixedParticlesTransmittanceAdvanced {
  public:
    // LayeredBxDFMixedParticlesTransmittanceAdvanced Public Methods
    LayeredBxDFMixedParticlesTransmittanceAdvanced() = default;
    PBRT_CPU_GPU
    LayeredBxDFMixedParticlesTransmittanceAdvanced(TopBxDF top, BottomBxDF bottom, CosmeticStructBxdf media1,CosmeticStructBxdf media2,  int maxDepth,int nSamples)
        : top(top),
          bottom(bottom),
          media1(media1),
          media2(media2),
          maxDepth(maxDepth),
          nSamples(nSamples),
          eta(1.3f) {
            
            medias[0] = media1;
            medias[1] = media2;

        // std::cout<<"std SIZE is " << std_size << std::endl;
    }



    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() {
        top.Regularize();
        bottom.Regularize();
    }

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        //TODO - remove this        
        // SampledSpectrum albedo(0.0f);
        BxDFFlags topFlags = top.Flags(), bottomFlags = bottom.Flags();
        CHECK(IsTransmissive(topFlags) ||
              IsTransmissive(bottomFlags));  // otherwise, why bother?

        BxDFFlags flags = BxDFFlags::Reflection;
        if (IsSpecular(topFlags))
            flags = flags | BxDFFlags::Specular;

        if (IsDiffuse(topFlags) || IsDiffuse(bottomFlags) || media1.albedo_diff || media1.albedo_flakes)
            flags = flags | BxDFFlags::Diffuse;
        else if (IsGlossy(topFlags) || IsGlossy(bottomFlags))
            flags = flags | BxDFFlags::Glossy;

        if (IsTransmissive(topFlags) && IsTransmissive(bottomFlags))
            flags = flags | BxDFFlags::Transmission;

        return flags;
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        // The top interface is a fake one, since it is a dielectric with Ior = 1.0 and roughness 0.0
        // The bottom interface is an actual dielectric with ior > 1.0 and roughness > 0.0
        SampledSpectrum f(0.);
        // Next line is for debugging purposes only, D.L.
        //  Estimate _LayeredBxDFMixedParticlesTransmittance_ value _f_ using random sampling
        //  Set _wo_ and _wi_ for layered BSDF evaluation
        /*if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }*/
        // Determine entrance interface for layered BSDF
        //std::cout<<"t1 "<< medias[0].thickness << " t2 "<< medias[1].thickness << std::endl;
        TopBxDF enterInterface;
       // return bottom.f(wo, wi, mode);
        bool enteredTop = twoSided || wo.z > 0;
        CosmeticStructBxdf selected_media;
        int media_index = 0;
        if (enteredTop){
            enterInterface = top;
            media_index = 1;
        }
        else{
            media_index = 0;
            enterInterface = bottom;
        }
        selected_media = medias[media_index];

        int exit_media_index;
        int last_index_media = 1;
        TopBxDF exitInterface, nonExitInterface;
        
        if (SameHemisphere(wo, wi) ^ enteredTop) {
            exitInterface = bottom;
            nonExitInterface = top;
            exit_media_index = 0;
            
        } else {
            exitInterface = top;
            nonExitInterface = bottom;
            exit_media_index = 1;

        }
        // exitInterface = top;
        //    nonExitInterface = bottom;
        //    exit_media_index = 1;
        if ((medias[0].thickness < MachineEpsilon && medias[0].thickness > -MachineEpsilon) && (medias[1].thickness < MachineEpsilon && medias[1].thickness > -MachineEpsilon)){
            return bottom.f(wo, wi, mode);
        }

   
        Float thickness = selected_media.thickness;

        Float exitZ = (SameHemisphere(wo, wi) ^ enteredTop) ? 0 : thickness;
        // Float exitZ = (SameHemisphere(wo, wi) ^ enteredTop) ? 0 : thickness;

        // Account for reflection at the entrance interface
        if (SameHemisphere(wo, wi))
        {   
            f = nSamples * enterInterface.f(wo, wi, mode);
        }        
        //return f/nSamples;
        // Declare _RNG_ for layered BSDF evaluation
        RNG rng(Hash(GetOptions().seed, wo), Hash(wi));
        
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        for (int s = 0; s < nSamples; ++s) {
            
            // Sample random walk through layers to estimate BSDF value
            // Sample transmission direction through entrance interface
            Float uc = r();
            pstd::optional<BSDFSample> wos = enterInterface.Sample_f(
                wo, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Transmission);
            if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0)
                continue;

            // Sample BSDF for virtual light from _wi_
            uc = r();
            pstd::optional<BSDFSample> wis = exitInterface.Sample_f(
                wi, uc, Point2f(r(), r()), !mode, BxDFReflTransFlags::Transmission);
            if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0)
                continue;
         
            // Declare state for random walk through BSDF layers
            SampledSpectrum beta = wos->f * AbsCosTheta(wos->wi) / wos->pdf;
            Float z = enteredTop ? thickness : 0;
            // std::cout<<" z is "<< z << std::endl;
            // Debugging D.L.
            //  wis->wi = -wi;
            //  wos->wi = -wo;
            // std::cout<< " " << AbsCosTheta(wos->wi) << "  " << AbsCosTheta(wo) << std::endl;
            
            Vector3f w = wos->wi;

            // SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g);
            
            // SGGXPhaseFunctionNormalDistribution phase(sigma,mean,std_or);

            // Diffuse phase function
            Point2f u1{r(), r()}, u2{r(), r()};
            // SampledSpectrum albedo_particle;

            for (int depth = 0; depth < maxDepth; ++depth) {
                selected_media = medias[media_index];
            SGGXPhaseFunctionMixedRefl phase(selected_media.platelets_roughness, selected_media.platelets_orientation, selected_media.std_or, selected_media.diff_conc,selected_media.g1,selected_media.g2,selected_media.w_g,selected_media.albedo_flakes,selected_media.albedo_diff);
                thickness = selected_media.thickness;
                // Sample next event for layered BSDF evaluation random walk
                PBRT_DBG("beta: %f %f %f %f, w: %f %f %f, f: %f %f %f %f\n", beta[0],
                         beta[1], beta[2], beta[3], w.x, w.y, w.z, f[0], f[1], f[2],
                         f[3]);
                // Possibly terminate layered BSDF random walk with Russian roulette
                 if (depth > 128 && beta.MaxComponentValue() < 0.01f) {
                 Float q = std::max<Float>(0, 1 - beta.MaxComponentValue());
                 if (r() < q)
                 break;
                beta /= 1 - q;
                PBRT_DBG("After RR with q = %f, beta: %f %f %f %f\n", q, beta[0],
                 beta[1], beta[2], beta[3]);
                }
                // Account for media between layers and possibly scatter
                // if (!true) {
                    //line after was the usual one but in any case this if-else is not used
                if (!true) {
                    // Advance to next layer boundary and update _beta_ for transmittance
                    z = (z == thickness) ? 0 : thickness;
                    beta *= Tr(thickness, w);
                } else {
                    // Sample medium scattering for layered BSDF evaluation
                    // Take into account the projected Area of particles along direction w
                    Float pa = phase.projectedArea(-w, Point2f(r(), r()));
                    //TODO:
                    //change this put the weighting factor
                    // sigma_t *= pa; //Only for SGGX
                    // sigma_t = diff_conc * sigma_t + (1-diff_conc)*(sigma_t * pa) ;
                    Float reduced_st = selected_media.diff_conc * selected_media.sigma_t + (1-selected_media.diff_conc)*(selected_media.sigma_t * pa) ;
                    
                    if (reduced_st <= 0) {
                        std::cout << "Projected area is NOT positive" << std::endl;
                        continue;
                    }
                    // std::cout<<" pa "<< pa << std::endl;
                    // Float distSampled = (sigmaT == 0.f) ? Infinity: -std::log(1.f -
                    // rng.Uniform<Float>() ) / sigmaT;
                    Float dz = SampleExponential(r(), reduced_st / std::abs(w.z));

                    Float zp = w.z > 0 ? (z + dz) : (z - dz);
                    // std::cout<<"zp " << zp << std::endl;
                    DCHECK_RARE(1e-5, z == zp);
                    if (z == zp)
                        continue;
                    if (0 < zp && zp < thickness) {
                        // Handle scattering event in layered BSDF medium
                        // Account for scattering through _exitInterface_ using _wis_
                        
                        Float wt = 1;
                        pstd::optional<PhaseFunctionSample> ps;
                        if (!IsSpecular(exitInterface.Flags()))
                                wt = PowerHeuristic(
                                    1, wis->pdf, 1,
                                    phase.PDF(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()) ) );
                            Float m_pdf = phase.PDF(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));

                            pa = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                            

                            
                            SampledSpectrum old_f = f;
                            Float pf_p = phase.p(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));
                            SampledSpectrum w_albedo = phase.albedoW(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));
                            // std::cout<<"w_albedo"<< w_albedo << std::endl;
                            auto sigmaTAlongWi = selected_media.diff_conc  + (1-selected_media.diff_conc)*( pa) ;
                            if (media_index == exit_media_index){
                                //we only need to traverse one media
                                f += beta * w_albedo *  wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            }else{
                                //we need to traverse the other media as well
                                auto sigmaTAlongWi_exit_media = medias[exit_media_index].diff_conc  + (1-medias[exit_media_index].diff_conc)*( pa) ;
                                f += beta * w_albedo *  wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * Tr(medias[exit_media_index].thickness, wis->wi, sigmaTAlongWi_exit_media) * wis->f /wis->pdf;

                            }
                            // f += beta * w_albedo * pf_p * wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            // f += beta * pf_p * wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;
                            
                            SampledSpectrum albedo_particle;
                            ps = phase.Sample_p(-w, Point2f(r(), r()),Point2f(r(), r()),r(),Point2f(r(), r()));
                            if (ps->p_flag == _MediumFlags::Diffuser){
                                albedo_particle = selected_media.albedo_diff;
                                // std::cout<<"albedo_diff "<< albedo_diff << std::endl;
                            }else{
                                albedo_particle = selected_media.albedo_flakes;
                                // std::cout<<"albedo_flakes "<< albedo_flakes << std::endl;

                            }


                            if (IsNaN(f[0])) {
                                std::cout << " l3170" << std::endl;
                                
                                std::cout << " nan found in f()_s_05" << std::endl;
                                std::cout << "oldF " << old_f << " pf_p " << pf_p
                                          << std::endl;
                                std::cout << "beta " << beta << " albedo " << w_albedo
                                          << " wt " << wt << " wis->f " << wis->f
                                          << " wis->pdf" << wis->pdf << std::endl;
                                std::cout << "Tr "
                                          << Tr(zp - exitZ, wis->wi, sigmaTAlongWi)
                                          << std::endl;
                            }
   

                        if (!ps || ps->pdf == 0 || ps->wi.z == 0 || ps->p == 0.0f)
                            continue;
                        // original line in code was the following one --> replaced since ps->p and ps->pdf cancels out as they are the same
                        // beta *= albedo * ps->p / ps->pdf;
                        beta *= albedo_particle;
                        w = ps->wi;
                        z = zp;

                        // Possibly account for scattering through _exitInterface_
                        if (((z < exitZ && w.z > 0 ) || (z > exitZ && w.z < 0 )) && !IsSpecular(exitInterface.Flags() )){                        
                        // if (((z < exitZ && w.z > 0) || (z > exitZ && w.z < 0)) && !IsSpecular(exitInterface.Flags())) {
                            // Account for scattering through _exitInterface_
                            SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                            if (fExit) {
                                Float exitPDF = exitInterface.PDF(-w, wi, mode, BxDFReflTransFlags::Transmission);
                                Float wt = PowerHeuristic(1, ps->pdf, 1, exitPDF);
                                // auto sigmaTAlongWi = phase.projectedArea(-ps->wi, Point2f(r(), r()));
                                pa = phase.projectedArea(-ps->wi, Point2f(r(), r()));
                                if (w.z > 0){
                                    //in case the ray is pointing upward
                                    sigmaTAlongWi = selected_media.diff_conc + (1-selected_media.diff_conc)*(1 * pa) ;
                                    if (media_index < last_index_media ){
                                        //in case we are in the bottom layer, we have to account for the top layer
                                        auto sigmaTAlongWi_top_media = medias[media_index+1].diff_conc + (1-medias[media_index+1].diff_conc)*(1 * pa) ;

                                        f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi)*Tr(medias[media_index+1].thickness, ps->wi, sigmaTAlongWi_top_media) *fExit * wt;
                                    }else{
                                        f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi) *fExit * wt;
                                    }
                                }
                                if (w.z < 0){
                                    //in case the ray is pointing upward
                                    sigmaTAlongWi = selected_media.diff_conc + (1-selected_media.diff_conc)*(1 * pa) ;
                                    if (media_index > 0 ){
                                        //in case we are in the top layer, we have to account for the bottom layer
                                        auto sigmaTAlongWi_bot_media = medias[media_index-1].diff_conc + (1-medias[media_index-1].diff_conc)*(1 * pa) ;

                                        f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi)*Tr(medias[media_index-1].thickness, ps->wi, sigmaTAlongWi_bot_media) *fExit * wt;
                                    }else{
                                        f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi) *fExit * wt;
                                    }
                                }
                                
                                // f += beta * Tr(zp - exitZ, ps->wi) * fExit * wt;
                            }
                        }

                        continue;
                    }
                    z = Clamp(zp, 0, thickness);
                    if ( z == 0){
                        media_index--;
                        if (media_index>-1){
                            selected_media = medias[media_index];
                            z = selected_media.thickness;
                            if (exit_media_index == 0){
                                exitZ = 0;
                            }
                            continue;
                        }
                    }else{
                        media_index++;
                        if (media_index<2){
                            selected_media = medias[media_index];
                            z = 0;
                            if (exit_media_index == 1 ){
                                exitZ = selected_media.thickness;
                            }
                            continue;
                        }

                    }
                }
                if (media_index == 2 || media_index == -1){
                    media_index = (media_index == 2) ? 1 : 0;
                }else{
                    std::cout<<"media index " << media_index << " thickness"<<thickness<<std::endl;
                }
                //TODO fix the next part
                // Account for scattering at appropriate interface
                if (z == exitZ ) {
                    // if (media_index == exit_media_index){
                    //     std::cout<<"Should not happen"<<std::endl;
                    // }
                    // Account for reflection at _exitInterface_
                    Float uc = r();
                    pstd::optional<BSDFSample> bs = exitInterface.Sample_f(
                        -w, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0){
                        break;
                    }
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;
                    w = bs->wi;

                } else {
                    // Account for scattering at _nonExitInterface_
                    if (!IsSpecular(nonExitInterface.Flags())) {
                        // Add NEE contribution along presampled _wis_ direction
                        Float wt = 1;
                        if (!IsSpecular(exitInterface.Flags()))
                            wt = PowerHeuristic(1.0f, wis->pdf, 1.0f,
                                                nonExitInterface.PDF(-w, -wis->wi, mode));

                        auto pa = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                        
                        auto sigmaTAlongWi = selected_media.diff_conc  + (1.0f-selected_media.diff_conc)*(pa) ;
                        auto sigmaTAlongWi_exit_media = medias[exit_media_index].diff_conc  + (1.0f-medias[exit_media_index].diff_conc)*(pa) ;

                        // f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                        //      AbsCosTheta(wis->wi) * wt *
                        //      Tr(thickness, wis->wi, sigmaTAlongWi) * wis->f / wis->pdf;
                        
                        f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                             AbsCosTheta(wis->wi) * wt *
                             Tr(thickness, wis->wi, sigmaTAlongWi)* Tr(medias[exit_media_index].thickness, wis->wi, sigmaTAlongWi_exit_media) * wis->f / wis->pdf;
                    }
                    // Sample new direction using BSDF at _nonExitInterface_
                    Float uc = r();
                    Point2f u(r(), r());
                    pstd::optional<BSDFSample> bs = nonExitInterface.Sample_f(
                        -w, uc, u, mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        break;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;

                    w = bs->wi;

                    if (!IsSpecular(exitInterface.Flags())) {
                        // Add NEE contribution along direction from BSDF sample
                        SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                        if (fExit) {
                            Float wt = 1;
                            if (!IsSpecular(nonExitInterface.Flags())) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                wt = PowerHeuristic(1, bs->pdf, 1, exitPDF);
                            }
                            auto pa = phase.projectedArea(-bs->wi, Point2f(r(), r()));
                            // auto sigmaTAlongWi = phase.projectedArea(-bs->wi, Point2f(r(), r()));
                            //TODO fix this
                            auto sigmaTAlongWi = selected_media.diff_conc * 1 + (1-selected_media.diff_conc)*(1 * pa) ;

                            f += beta * Tr(thickness, bs->wi, sigmaTAlongWi) * fExit * wt;
                        }
                    }
                }
            }
        }
        if (f.HasNaNs()) {
            std::cout << " nan found in f()" << std::endl;
        }
        return f / nSamples;
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(

        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now

//        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
//           return {};
        bool enteredTop = twoSided || wo.z > 0;

        // Set _wo_ for layered BSDF sampling
        bool flipWi = false;
        if (twoSided && wo.z < 0) {
            wo = -wo;
            flipWi = true;
        }
        //TODO remove, used only for debug
        //return bottom.Sample_f(wo, uc, u, mode);
        // Sample BSDF at entrance interface to get initial direction _w_
        pstd::optional<BSDFSample> bs = enteredTop ? top.Sample_f(wo, uc, u, mode) : bottom.Sample_f(wo, uc, u, mode);
        
       // return bs;

        if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
            return {};
        if (bs->IsReflection()) {
            // if(!enteredTop){
                // std::cout<<"SHOULD NOT HAPPEN"<<std::endl;
            // }
            if (flipWi)
                bs->wi = -bs->wi;
            bs->pdfIsProportional = true;
            return bs;
        }
        
        Vector3f w = bs->wi;
        // std::cout<<"bs->wi "<< w << " wo "<< wo << std::endl ;
        bool specularPath = bs->IsSpecular();

        // Declare _RNG_ for layered BSDF sampling
        RNG rng(Hash(GetOptions().seed, wo), Hash(uc, u));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };
        
        // size_t media_index = 0;
        size_t last_index = 1;
               
        // Declare common variables for layered BSDF sampling
        SampledSpectrum f = bs->f * AbsCosTheta(bs->wi);
                if (IsInf(f[0])){
                    std::cout<< "it's INF line 1900"<< std::endl;
                }
        Float pdf = bs->pdf;
        size_t media_index = enteredTop ? 1 : 0;

        Float thickness = medias[media_index].thickness;
        Float z = enteredTop ? thickness : 0;
        // std::cout<<"media1.diffprec" << medias[0].diff_conc << "media2.diffprce " << medias[1].diff_conc<< std::endl; 
        // phase.toString();
        // SGGXPhaseFunctionMixedRefl phase(media1.sigma, media1.mean, media1.std_or, media1.diff_conc,media1.g1,media1.g2,media1.w_g,media1.albedo_flakes,media1.albedo_diff);
        // std::cout<<"maxDepth "<< maxDepth << " thick "<< thickness << " sigmat " << sigma_t << std::endl;
        if ((medias[0].thickness < MachineEpsilon && medias[0].thickness > -MachineEpsilon) && (medias[1].thickness < MachineEpsilon && medias[1].thickness > -MachineEpsilon)){
            //just get the value of the bottomLayer material
            return bottom.Sample_f(wo, uc, u, mode);
        }
        for (int depth = 0; depth < maxDepth; ++depth) {
            // std::cout<<" d "<< depth << std::endl;
            // Follow random walk through layers to sample layered BSDF
            // Possibly terminate layered BSDF sampling with Russian Roulette
            Float rrBeta = f.MaxComponentValue() / pdf;
            if (depth > 128 && rrBeta < 0.01f) {
            Float q = std::max<Float>(0, 1 - rrBeta);
            if (r() < q)
            return {};
            pdf *= 1 - q;
            }
            if (w.z == 0)
                return {};
            //if one of the two is selected than it makes sense to use this
            
            CosmeticStructBxdf selected_media = medias[media_index];
            thickness = selected_media.thickness;
            //TODO fix here
            if (thickness == 0.0){
                
                    if (w.z > 0){
                        media_index++;
                    }else{
                        media_index--;
                    }
                    if (media_index > -1 && media_index<2){
                        if (w.z>0){
                            z = 0;
                        }else{
                            z = medias[media_index].thickness;
                        }
                        
                                            //this media has no thickness go to the next one
                    continue;
                    }
                    /*else{

                    if (media_index == -1){
                            //TODO
                            //if this reflect, than we should calculate the path going upward 
                            //last media, use the bottom interface
                    bs = bottom.Sample_f(-w, uc, u, mode);
                     if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        return {};
                    w = bs->wi;
                    f *= bs->f;
                    pdf *= bs->pdf;
                    
                   
                    if (bs->IsReflection()){
                        //reflected back?--> then continue
                        continue;
                    }else{
                        //transmission event when bottom layer? --> BSSRDF
                        BxDFFlags flags = BxDFFlags::Transmission;
                        return BSDFSample(f, w, pdf, flags,1.3f);
                    }
                        }
                    else{
                    //Top layer? Transmission is the only way
                    BxDFFlags flags =  BxDFFlags::Transmission;
                    //bool specularPath =false;
                    //specularPath&= bs->IsSpecular();
                    //flags |= specularPath ? BxDFFlags::Specular : BxDFFlags::Glossy;
                   // return bs;
                    bs = top.Sample_f(-w, uc, u, mode);
                         if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        return {};
                    w = bs->wi;
                    f *= bs->f;
                    pdf *= bs->pdf;
                        return BSDFSample(f, w, pdf, flags,1.0f);
                    }
                }*/


            }
            Float sigma_t = selected_media.sigma_t;

            SGGXPhaseFunctionMixedRefl phase(selected_media.platelets_roughness, selected_media.platelets_orientation, selected_media.std_or, selected_media.diff_conc,selected_media.g1,selected_media.g2,selected_media.w_g,selected_media.albedo_flakes,selected_media.albedo_diff);
                   
            if (selected_media.albedo_diff || selected_media.albedo_flakes) {
                // Sample potential scattering event in layered medium
                Float pa = phase.projectedArea(w, Point2f(r(), r()));
                
                Float reduced_st = selected_media.diff_conc * sigma_t + (1-selected_media.diff_conc)*(sigma_t * pa) ;

                Float dz = SampleExponential(r(), reduced_st / AbsCosTheta(w));
                Float zp = w.z > 0 ? (z + dz) : (z - dz);
                CHECK_RARE(1e-5, zp == z);
                if (zp == z)
                    return {};
                //media interaction
                if (0 < zp && zp < thickness) {
                    pstd::optional<PhaseFunctionSample> ps;
                    //get the albedo of the sampled particle
                    SampledSpectrum albedo_particle;
                    ps = phase.Sample_p(-w, Point2f(r(), r()),Point2f(r(), r()),r(),Point2f(r(), r())); 
                    if (ps->p_flag == _MediumFlags::Diffuser){
                                
                                albedo_particle = selected_media.albedo_diff;
                            }else{
                                // std::cout<<"happens"<<std::endl;
                                albedo_particle = selected_media.albedo_flakes;
                    }
                    // if (!albedo_particle){
                        // LOG_ERROR("Albedo particle is not instantiated");
                    // } 
                    if (!ps || ps->pdf == 0 || ps->wi.z == 0){
                        // LOG_ERROR("early termination is not great",pPixel.x, pPixel.y, sampleIndex);
                        return {};

                    }
                    auto old_f = f;
                    
                    if (ps->p_flag == _MediumFlags::Diffuser){
                        f *= albedo_particle;
                    }else{
                        // f *= albedo_particle;
                        
                        f *= albedo_particle * pa;
                    
                    }
                    if (IsInf(f[0])){
                    std::cout<< "it's INF"<< std::endl;
                    std::cout<<"albedo particle" << albedo_particle << " pa "<< pa << media_index << std::endl;
            }
                    
                    specularPath = false;
                    w = ps->wi;
                    z = zp;
                    continue;
                }
                z = Clamp(zp, 0, thickness);
                
                if (z == 0){
                    DCHECK_LT(w.z, 0);
                    media_index --;
                    //going down to 1 media
                    if (media_index > -1){
                    thickness = medias[media_index].thickness;
                    z = thickness;
                        
                    }

                }
                else{
                    media_index++;
                    DCHECK_GT(w.z, 0);
                    if (media_index < last_index+1){
                    //moving upward
                    thickness = medias[media_index].thickness;
                    z = 0;
                        
                    }

                }
                    if (IsInf(f[0])){
                    std::cout<< "it's INF line 3478"<< std::endl;
            }
                } else {
                    // Advance to the other layer interface
                    z = (z == thickness) ? 0 : thickness;
                    f *= Tr(thickness, w);
                }
            if (media_index > -1 || media_index < last_index+1){
                //sample is still inside one of the medias
                continue;
            }else{
                //clip to legal values
                // std::cout<<"media index "<< media_index<< std::endl;
                media_index = (media_index <= -1) ? 0 : 1;
            }
            // Initialize _interface_ for current interface surface
#ifdef interface  // That's enough out of you, Windows.
#undef interface
#endif
            BottomBxDF interface;
            bool isBot = false;
            if ( z == 0){
                isBot = true;
                interface = bottom;
            }
            else{
                interface = top;
            }

            // Sample interface BSDF to determine new path direction
            Float uc = r();
            Point2f u(r(), r());
            
            // std::cout<<"w "<< w << std::endl ;
            pstd::optional<BSDFSample> bs = interface.Sample_f(-w, uc, u, mode);
            if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0){
                // std::cout<<"-w " << -w << " uc "<< uc << " u " << u << " mode " << mode << std::endl;
                return {};
            }

            f *= bs->f;

            pdf *= bs->pdf;
            specularPath &= bs->IsSpecular();
            w = bs->wi;
            
            if (bs->IsTransmission() ) {
                BxDFFlags flags = SameHemisphere(wo, w) ? BxDFFlags::Reflection
                                                        : BxDFFlags::Transmission;
                flags |= specularPath ? BxDFFlags::Specular : BxDFFlags::Glossy;
                if (flipWi)
                    w = -w;
                if (f.HasNaNs()) {
                    std::cout << " nan found in sample()" << std::endl;
                }
                if (IsInf(f[0])){
                    std::cout<< "it's INF"<< std::endl;
                }
                // std::cout<<"exiting f "<< f << " w " << w <<  " depth "<< depth << std::endl; 
                if ( z == thickness){
                    //scatter towards the air
                return BSDFSample(f, w, pdf, flags,1.0f);
                }else{
                    //scatter inside the skin
                return BSDFSample(f, w, pdf, flags,1.3);

                }
                // Old Code
                // return BSDFSample(f, w, pdf, flags, 1.3f, true);
            }
            // Scale _f_ by cosine term after scattering at the interface
            f *= AbsCosTheta(bs->wi);
        }
        //TODO check why this is happening
        //LOG_ERROR("End Reached");
        return {};
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // Set _wo_ and _wi_ for layered BSDF evaluation
        bool enteredTop = twoSided || wo.z > 0;

        //****************************************************
        //                   Debug code
        //****************************************************
        
        //  if(!enteredTop){
        //     // return top.PDF(wo, wi, mode, sampleFlags);
        //     wo = -wo;
        //     if (!(sampleFlags & BxDFReflTransFlags::Reflection))
        //         return 0;
        //     return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;

        //  }

        //****************************************************
        //                   END Debug code
        //****************************************************

        // return bottom.PDF(wo,wi,mode,sampleFlags);
        
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }
                
        // Declare _RNG_ for layered PDF evaluation
        RNG rng(Hash(GetOptions().seed, wi), Hash(wo));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Update _pdfSum_ for reflection at the entrance layer
        Float pdfSum = 0;
        if (SameHemisphere(wo, wi)) {
            auto reflFlag = BxDFReflTransFlags::Reflection;
            pdfSum += enteredTop ? nSamples * top.PDF(wo, wi, mode, reflFlag)
                                 : nSamples * bottom.PDF(wo, wi, mode, reflFlag);
        }

        for (int s = 0; s < nSamples; ++s) {
            // Evaluate layered BSDF PDF sample
            if (SameHemisphere(wo, wi)) {
                // Evaluate TRT term for PDF estimate
                // Old Code
                // TopOrBottomBxDF<TopBxDF, BottomBxDF> rInterface, tInterface;
                // const TopBxDF *rInterface, *tInterface;
                TopBxDF rInterface, tInterface;
                
                if (enteredTop) {
                    rInterface = bottom;
                    tInterface = top;
                } else {
                    rInterface = top;
                    tInterface = bottom;
                }
                // Sample _tInterface_ to get direction into the layers
                auto trans = BxDFReflTransFlags::Transmission;
                pstd::optional<BSDFSample> wos, wis;
                wos = tInterface.Sample_f(wo, r(), {r(), r()}, mode, trans);
                wis = tInterface.Sample_f(wi, r(), {r(), r()}, !mode, trans);

                // Update _pdfSum_ accounting for TRT scattering events
                if (wos && wos->f && wos->pdf > 0 && wis && wis->f && wis->pdf > 0) {
                    if (!IsNonSpecular(tInterface.Flags()))
                        pdfSum += rInterface.PDF(-wos->wi, -wis->wi, mode);
                    else {
                        // Use multiple importance sampling to estimate PDF product
                        pstd::optional<BSDFSample> rs =
                            rInterface.Sample_f(-wos->wi, r(), {r(), r()}, mode);
                        if (rs && rs->f && rs->pdf > 0) {
                            if (!IsNonSpecular(rInterface.Flags()))
                                pdfSum += tInterface.PDF(-rs->wi, wi, mode);
                            else {
                                // Compute MIS-weighted estimate of Equation
                                // (\ref{eq:pdf-triple-canceled-one})
                                Float rPDF = rInterface.PDF(-wos->wi, -wis->wi, mode);
                                Float wt = PowerHeuristic(1, wis->pdf, 1, rPDF);
                                pdfSum += wt * rPDF;

                                Float tPDF = tInterface.PDF(-rs->wi, wi, mode);
                                wt = PowerHeuristic(1, rs->pdf, 1, tPDF);
                                pdfSum += wt * tPDF;
                            }
                        }
                    }
                }

            } else {
                
                // Evaluate TT term for PDF estimate
                // TopOrBottomBxDF<TopBxDF, BottomBxDF> toInterface, tiInterface;
                TopBxDF toInterface, tiInterface;
                
                if (enteredTop) {
                    toInterface = top;
                    tiInterface = bottom;
                } else {
                    toInterface = bottom;
                    tiInterface = top;
                }

                Float uc = r();
                Point2f u(r(), r());
                
                pstd::optional<BSDFSample> wos = toInterface.Sample_f(wo, uc, u, mode);
                if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0 ||
                    wos->IsReflection())
                    continue;

                uc = r();
                u = Point2f(r(), r());
                pstd::optional<BSDFSample> wis = tiInterface.Sample_f(wi, uc, u, !mode);
                if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0 ||
                    wis->IsReflection())
                    continue;

                if (IsSpecular(toInterface.Flags()))
                    pdfSum += tiInterface.PDF(-wos->wi, wi, mode);
                else if (IsSpecular(tiInterface.Flags()))
                    pdfSum += toInterface.PDF(wo, -wis->wi, mode);
                else
                    pdfSum += (toInterface.PDF(wo, -wis->wi, mode) +
                               tiInterface.PDF(-wos->wi, wi, mode)) /
                              2;
            }
        }
        // Return mixture of PDF estimate and constant PDF
        return Lerp(0.9f, 1 / (4 * Pi), pdfSum / nSamples);
    }

  private:
    // LayeredBxDFMixedParticlesTransmittanceAdvanced Private Methods
    PBRT_CPU_GPU
    static Float Tr(Float dz, Vector3f w) {
        if (std::abs(dz) <= std::numeric_limits<Float>::min())
            return 1;
        return FastExp(-std::abs(dz / w.z));
    }
    static Float Tr(Float dz, Vector3f w, Float sigmaTAlongWi) {
        return FastExp(-std::abs(sigmaTAlongWi * dz / w.z));
    }

    // LayeredBxDFMixedParticlesTransmittanceAdvanced Private Members
    TopBxDF top;
    BottomBxDF bottom;

    // Float thickness, g1,g2,w_g, sigma, mean, std_or, std_size, diff_conc;
    // Float sigma_t,eta;
    // Vector3f dir;
    // SampledSpectrum albedo_diff,albedo_flakes;
    Float eta;
    CosmeticStructBxdf media1,media2;
    CosmeticStructBxdf medias[2];
    int maxDepth, nSamples;
    bool use_specular;
};


template <typename TopBxDF, typename BottomBxDF, bool twoSided>
class LayeredBxDFMixedParticlesNormalizedFresnel {
  public:
    // LayeredBxDFMixedParticlesNormalizedFresnel Public Methods
    LayeredBxDFMixedParticlesNormalizedFresnel() = default;
    PBRT_CPU_GPU
    LayeredBxDFMixedParticlesNormalizedFresnel(TopBxDF top, BottomBxDF bottom, Float thickness,
                              Vector3f dir, Float sigma, Float mean, Float std_or,
                              Float std_size, Float diff_conc, bool use_refl,
                              const SampledSpectrum &albedo_flakes,const SampledSpectrum &albedo_diff, Float g1,Float g2,Float w_g, int maxDepth,
                              int nSamples)
        : top(top),
          bottom(bottom),
          thickness(std::max(thickness, std::numeric_limits<Float>::min())),
          dir(dir),
          sigma(sigma),
          mean(mean),
          std_or(std_or),
          std_size(std_size),
          use_specular(use_refl),
          g1(g1),
          g2(g2),
          w_g(w_g),
          albedo_flakes(albedo_flakes),
          albedo_diff(albedo_diff),
          maxDepth(maxDepth),
          nSamples(nSamples),
          diff_conc(diff_conc),
          sigma_t(1.0f),
          eta(1.3f) {
            
        // std::cout<<"std SIZE is " << std_size << std::endl;
    }

    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() {
        top.Regularize();
        bottom.Regularize();
    }

    PBRT_CPU_GPU
    BxDFFlags Flags() const {
        //TODO - remove this        
        // SampledSpectrum albedo(0.0f);
        BxDFFlags topFlags = top.Flags(), bottomFlags = bottom.Flags();
        CHECK(IsTransmissive(topFlags) ||
              IsTransmissive(bottomFlags));  // otherwise, why bother?

        BxDFFlags flags = BxDFFlags::Reflection;
        if (IsSpecular(topFlags))
            flags = flags | BxDFFlags::Specular;

        if (IsDiffuse(topFlags) || IsDiffuse(bottomFlags) || albedo_diff || albedo_flakes)
            flags = flags | BxDFFlags::Diffuse;
        else if (IsGlossy(topFlags) || IsGlossy(bottomFlags))
            flags = flags | BxDFFlags::Glossy;

        if (IsTransmissive(topFlags) && IsTransmissive(bottomFlags))
            flags = flags | BxDFFlags::Transmission;

        return flags;
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const {
        SampledSpectrum w_albedo2;
        



         w_albedo2 = diff_conc * albedo_diff + (1-diff_conc) * albedo_flakes;
                            // if (ps->p_flag == _MediumFlags::Diffuser){
                            //     albedo_particle = albedo_diff;
                            // }else{
                            //     albedo_particle = albedo_flakes;
                            // }
        // SampledSpectrum f(0.);
        // Next line is for debugging purposes only, D.L.
        //  wo = Vector3f(0,0,1);
        //  wi = Vector3f(0,0,1);
        //  Estimate _LayeredBxDFMixedParticlesNormalizedFresnel_ value _f_ using random sampling
        //  Set _wo_ and _wi_ for layered BSDF evaluation
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }
        
        // Determine entrance interface for layered BSDF
        //Old Code
        TopOrBottomBxDF<TopBxDF, BottomBxDF> enterInterface;
        // TopBxDF enterInterface;

        bool enteredTop = twoSided || wo.z > 0;
        if (enteredTop)
            enterInterface = &top;
        else
            enterInterface = &bottom;
        
        //****************************************************
        //                   Debug code
        //****************************************************
        //  if(!enteredTop){
            // return top.f(wo, wi, mode);
            // wo = -wo;
            // if (!SameHemisphere(wo, wi))
            //             return SampledSpectrum(0.f);
                    // Compute $\Sw$ factor for BSSRDF value


                    // return f;
        //  }
        //****************************************************
        //                   END Debug code
        //****************************************************

        // Determine exit interface and exit $z$ for layered BSDF
        
        //Old Code
        TopOrBottomBxDF<TopBxDF, BottomBxDF> exitInterface, nonExitInterface;
                // if (SameHemisphere(wo, wi) ^ enteredTop) {
        //     exitInterface = &top;
        //     nonExitInterface = &top;
        // } else {
        //     exitInterface = &top;
        //     nonExitInterface = &top;
        // }
        Float c = 1 - 2 * FresnelMoment1(1 / eta);
        SampledSpectrum f((1 - FrDielectric(CosTheta(wi), eta)) / (c * Pi));

                    // Update BSSRDF transmission term to account for adjoint light transport
        if (mode == TransportMode::Radiance)
            f *= Sqr(eta);


        // TopBxDF exitInterface, nonExitInterface;
        if (SameHemisphere(wo, wi) ^ enteredTop) {
            exitInterface = &bottom;
            nonExitInterface = &top;
        } else {
            exitInterface = &top;
            nonExitInterface = &bottom;
        }
        Float exitZ = (SameHemisphere(wo, wi) ^ enteredTop) ? 0 : thickness;

        // Account for reflection at the entrance interface
        if (SameHemisphere(wo, wi))
        {   
            //changed here
            f *= nSamples * enterInterface.f(wo, wi, mode);
        }        

        // Declare _RNG_ for layered BSDF evaluation
        RNG rng(Hash(GetOptions().seed, wo), Hash(wi));

        
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };
        // std::cout << "wi "<< wi << std::endl;
        // std::cout << "wo "<< wo << std::endl;
        // auto t_seed = GetOptions().seed;
        // std::cout<< t_seed << std::endl;
        // std::cout<<"first r " << r() <<std::endl;
        for (int s = 0; s < nSamples; ++s) {
            // Sample random walk through layers to estimate BSDF value
            // Sample transmission direction through entrance interface
            Float uc = r();
            pstd::optional<BSDFSample> wos = enterInterface.Sample_f(
                wo, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Transmission);
            if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0)
                continue;

            // Sample BSDF for virtual light from _wi_
            uc = r();
            pstd::optional<BSDFSample> wis = exitInterface.Sample_f(
                wi, uc, Point2f(r(), r()), !mode, BxDFReflTransFlags::Transmission);
            if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0)
                continue;
            // std::cout<<"wos->wi" << wos->wi <<" wi " << wi << std::endl;
            // std::cout<<"wis->wi" << wis->wi <<" wo " << wo << std::endl;



            // Declare state for random walk through BSDF layers
            SampledSpectrum beta = wos->f * AbsCosTheta(wos->wi) / wos->pdf;
            Float z = enteredTop ? thickness : 0;
            // std::cout<<" z is "<< z << std::endl;
            // Debugging D.L.
            //  wis->wi = -wi;
            //  wos->wi = -wo;
            // std::cout<< " " << AbsCosTheta(wos->wi) << "  " << AbsCosTheta(wo) << std::endl;
            
            Vector3f w = wos->wi;

            // SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g);
            SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g,albedo_flakes,albedo_diff);
            
            // SGGXPhaseFunctionNormalDistribution phase(sigma,mean,std_or);

            // Diffuse phase function
            Point2f u1{r(), r()}, u2{r(), r()};
            bool is_specular_phase = use_specular;
            // SampledSpectrum albedo_particle;

            for (int depth = 0; depth < maxDepth; ++depth) {
                // Sample next event for layered BSDF evaluation random walk
                PBRT_DBG("beta: %f %f %f %f, w: %f %f %f, f: %f %f %f %f\n", beta[0],
                         beta[1], beta[2], beta[3], w.x, w.y, w.z, f[0], f[1], f[2],
                         f[3]);
                // Possibly terminate layered BSDF random walk with Russian roulette
                 if (depth > 128 && beta.MaxComponentValue() < 0.01f) {
                 Float q = std::max<Float>(0, 1 - beta.MaxComponentValue());
                 if (r() < q)
                 break;
                beta /= 1 - q;
                PBRT_DBG("After RR with q = %f, beta: %f %f %f %f\n", q, beta[0],
                 beta[1], beta[2], beta[3]);
                }
                // Account for media between layers and possibly scatter
                // if (!true) {
                    //line after was the usual one but in any case this if-else is not used
                if (!true) {
                    // Advance to next layer boundary and update _beta_ for transmittance
                    z = (z == thickness) ? 0 : thickness;
                    beta *= Tr(thickness, w);
                } else {
                    // Sample medium scattering for layered BSDF evaluation
                    // Take into account the projected Area of particles along direction w
                    Float pa = phase.projectedArea(-w, Point2f(r(), r()));
                    //TODO:
                    //change this put the weighting factor
                    // sigma_t *= pa; //Only for SGGX
                    // sigma_t = diff_conc * sigma_t + (1-diff_conc)*(sigma_t * pa) ;
                    Float reduced_st = diff_conc * sigma_t + (1-diff_conc)*(sigma_t * pa) ;
                    
                    if (reduced_st <= 0) {
                        std::cout << "Projected area is NOT positive" << std::endl;
                        continue;
                    }
    
                    Float dz = SampleExponential(r(), reduced_st / std::abs(w.z));

                    Float zp = w.z > 0 ? (z + dz) : (z - dz);
                    DCHECK_RARE(1e-5, z == zp);
                    if (z == zp)
                        continue;
                    if (0 < zp && zp < thickness) {
                        // Handle scattering event in layered BSDF medium
                        // Account for scattering through _exitInterface_ using _wis_
                        Float wt = 1;
                        pstd::optional<PhaseFunctionSample> ps;
                        if (!IsSpecular(exitInterface.Flags()))
                                wt = PowerHeuristic(
                                    1, wis->pdf, 1,
                                    phase.PDF(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()) ) );
                            Float m_pdf = phase.PDF(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));

                            pa = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                            auto sigmaTAlongWi = diff_conc  + (1-diff_conc)*( pa) ;
                            
                            SampledSpectrum old_f = f;
                            Float pf_p = phase.p(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));
                            SampledSpectrum w_albedo = phase.albedoW(-w, -wis->wi, Point2f(r(), r()),Point2f(r(), r()));
                            f += beta * w_albedo *  wt *Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /wis->pdf;

                            
                            SampledSpectrum albedo_particle;
                            ps = phase.Sample_p(-w, Point2f(r(), r()),Point2f(r(), r()),r(),Point2f(r(), r()));
                            if (ps->p_flag == _MediumFlags::Diffuser){
                                albedo_particle = albedo_diff;
                                // std::cout<<"albedo_diff "<< albedo_diff << std::endl;
                            }else{
                                albedo_particle = albedo_flakes;
                                // std::cout<<"albedo_flakes "<< albedo_flakes << std::endl;

                            }

                            if (IsNaN(f[0])) {
                                std::cout<<" line3916"<<std::endl;
                                std::cout << " nan found in f()_s_05" << std::endl;
                                std::cout << "oldF " << old_f << " pf_p " << pf_p
                                          << std::endl;
                                std::cout << "beta " << beta << " albedo " << w_albedo
                                          << " wt " << wt << " wis->f " << wis->f
                                          << " wis->pdf" << wis->pdf << std::endl;
                                std::cout << "Tr "
                                          << Tr(zp - exitZ, wis->wi, sigmaTAlongWi)
                                          << std::endl;
                                // std::cout<<" m_pdf is zero " << m_pdf << " -w " << -w
                                // << " -wis->wi " << -wis->wi<< " Good luck! "
                                // <<std::endl;
                            }
                            // std::cout<<"tdebug " << tdebug << std::endl;
                            // std::cout<<"ps->wi " << ps->wi << std::endl;
                         /*else {
                            // Diffuse phase function
                            if (!IsSpecular(exitInterface.Flags()))
                                wt = PowerHeuristic(
                                    1, wis->pdf, 1,
                                    phase.PDF(-w, -wis->wi, u1, Point2f(r(), r())));
                            auto sigmaTAlongWi = phase.projectedArea(
                                -wis->wi, Point2f(r(), r()), Point2f(r(), r()));
                            f += beta * albedo *
                                 phase.p(-w, -wis->wi, Point2f(r(), r()),
                                         Point2f(r(), r())) *
                                 wt * Tr(zp - exitZ, wis->wi, sigmaTAlongWi) * wis->f /
                                 wis->pdf;
                            // f += beta * albedo * phase.p(-w, -wis->wi,
                            // u1,Point2f(r(),r())) * wt *
                            //      Tr(zp - exitZ, wis->wi) * wis->f / wis->pdf;

                            // Sample phase function and update layered path state
                            ps = phase.Sample_p(-w, Point2f(r(), r()), Point2f(r(), r()),
                                                Point2f(r(), r()));
                        }*/

                        if (!ps || ps->pdf == 0 || ps->wi.z == 0 || ps->p == 0.0f)
                            continue;
                        // original line
                        // beta *= albedo * ps->p / ps->pdf;
                        beta *= albedo_particle;
                        //beta *=  w_albedo/pf_p;
                        // beta *=  w_albedo;
                        w = ps->wi;
                        z = zp;

                        // Possibly account for scattering through _exitInterface_
                        if (((z < exitZ && w.z > 0) || (z > exitZ && w.z < 0)) &&
                            !IsSpecular(exitInterface.Flags())) {
                            // Account for scattering through _exitInterface_
                            SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                            if (fExit) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                Float wt = PowerHeuristic(1, ps->pdf, 1, exitPDF);
                                // auto sigmaTAlongWi = phase.projectedArea(-ps->wi, Point2f(r(), r()));
                                pa = phase.projectedArea(-ps->wi, Point2f(r(), r()));
                                sigmaTAlongWi = diff_conc * 1 + (1-diff_conc)*(1 * pa) ;

                                f += beta * Tr(zp - exitZ, ps->wi, sigmaTAlongWi) *
                                     fExit * wt;

                                // f += beta * Tr(zp - exitZ, ps->wi) * fExit * wt;
                            }
                        }

                        continue;
                    }
                    z = Clamp(zp, 0, thickness);
                }

                // Account for scattering at appropriate interface
                if (z == exitZ) {
                    // Account for reflection at _exitInterface_
                    Float uc = r();
                    pstd::optional<BSDFSample> bs = exitInterface.Sample_f(
                        -w, uc, Point2f(r(), r()), mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0){
                        //debug
                        // std::cout<<"reached the end, break because we can't reflect"<<std::endl;
                        break;
                    }
                    // std::cout << "Reflection back: "<< bs->wi << " w "<< w <<  std::endl;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;
                    w = bs->wi;

                } else {
                    // Account for scattering at _nonExitInterface_
                    if (!IsSpecular(nonExitInterface.Flags())) {
                        // Add NEE contribution along presampled _wis_ direction
                        Float wt = 1;
                        if (!IsSpecular(exitInterface.Flags()))
                            wt = PowerHeuristic(1.0f, wis->pdf, 1.0f,
                                                nonExitInterface.PDF(-w, -wis->wi, mode));

                        // f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                        //  AbsCosTheta(wis->wi) * wt * Tr(thickness, wis->wi) * wis->f /
                        //  wis->pdf;
                        // auto sigmaTAlongWi = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                        auto pa = phase.projectedArea(-wis->wi, Point2f(r(), r()));
                        auto sigmaTAlongWi = diff_conc  + (1.0f-diff_conc)*(pa) ;
                        
                        f += beta * nonExitInterface.f(-w, -wis->wi, mode) *
                             AbsCosTheta(wis->wi) * wt *
                             Tr(thickness, wis->wi, sigmaTAlongWi) * wis->f / wis->pdf;
                    }
                    // Sample new direction using BSDF at _nonExitInterface_
                    Float uc = r();
                    Point2f u(r(), r());
                    pstd::optional<BSDFSample> bs = nonExitInterface.Sample_f(
                        -w, uc, u, mode, BxDFReflTransFlags::Reflection);
                    if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
                        break;
                    beta *= bs->f * AbsCosTheta(bs->wi) / bs->pdf;

                    w = bs->wi;

                    if (!IsSpecular(exitInterface.Flags())) {
                        // Add NEE contribution along direction from BSDF sample
                        SampledSpectrum fExit = exitInterface.f(-w, wi, mode);
                        if (fExit) {
                            Float wt = 1;
                            if (!IsSpecular(nonExitInterface.Flags())) {
                                Float exitPDF = exitInterface.PDF(
                                    -w, wi, mode, BxDFReflTransFlags::Transmission);
                                wt = PowerHeuristic(1, bs->pdf, 1, exitPDF);
                            }
                            auto pa = phase.projectedArea(-bs->wi, Point2f(r(), r()));
                            // auto sigmaTAlongWi = phase.projectedArea(-bs->wi, Point2f(r(), r()));
                            auto sigmaTAlongWi = diff_conc * 1 + (1-diff_conc)*(1 * pa) ;

                            f += beta * Tr(thickness, bs->wi, sigmaTAlongWi) * fExit * wt;
                        }
                    }
                }
            }
        }
        if (f.HasNaNs()) {
            std::cout << " nan found in f()" << std::endl;
        }
        return f / nSamples;
    }

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(
        Vector3f wo, Float uc, Point2f u, TransportMode mode,
        BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // return bottom.Sample_f(wo, uc, u, mode);

        if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            return {};
        bool enteredTop = twoSided || wo.z > 0;

        // Cosine-sample the hemisphere, flipping the direction if necessary
        //****************************************************
        //                   Debug code
        //****************************************************
        // if (!enteredTop){
            // return top.Sample_f(wo, uc, u, mode);
        // auto wo_temp = -wo;
        // Vector3f wi = SampleCosineHemisphere(u);
        // if (wo.z < 0)
        //     wi.z *= -1;
        
        // return BSDFSample(f(wo, wi, mode), wi, PDF(wo, wi, mode, sampleFlags),
        //                   BxDFFlags::DiffuseReflection);
        //  }
        //****************************************************
        //                   Debug code
        //****************************************************

        // Set _wo_ for layered BSDF sampling
        bool flipWi = false;
        if (twoSided && wo.z < 0) {
            wo = -wo;
            flipWi = true;
        }

        // Sample BSDF at entrance interface to get initial direction _w_
        
        // pstd::optional<BSDFSample> bs =enteredTop ? top.Sample_f(wo, uc, u, mode) : bottom.Sample_f(wo, uc, u, mode,BxDFReflTransFlags::Reflection);
        
        //For this material (that is used only with a BSSRDF), we always sample the bottom interface as starting point because we
        pstd::optional<BSDFSample> bs = bottom.Sample_f(wo, uc, u, mode,BxDFReflTransFlags::Reflection);

        // if (!enteredTop){
        //     // std::cout<<"this case"<<std::endl;
        //     bs = bottom.Sample_f(wo, uc, u, mode,BxDFReflTransFlags::Transmission);
        // }
        if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0)
            return {};
        /*if (bs->IsReflection()) {
            if (flipWi)
                bs->wi = -bs->wi;
            bs->pdfIsProportional = true;
            return bs;
        }*/
        Vector3f w = bs->wi;
        // std::cout<<"bs->wi "<< w << " wo "<< wo << std::endl ;
        bool specularPath = bs->IsSpecular();

        // Declare _RNG_ for layered BSDF sampling
        RNG rng(Hash(GetOptions().seed, wo), Hash(uc, u));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Declare common variables for layered BSDF sampling
        SampledSpectrum f = bs->f * AbsCosTheta(bs->wi);
                        if (IsInf(f[0])){
                    std::cout<< "it's INF line 1900"<< std::endl;
                }
        Float pdf = bs->pdf;
        Float z = enteredTop ? thickness : 0;

        // std::cout<< "z "<< z << " thickness "<< thickness << std::endl;
        // HGPhaseFunction phase(g);//Original code: Henyey-Greenstein phase function
        SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, diff_conc,g1,g2,w_g,albedo_flakes,albedo_diff);
        // SGGXPhaseFunctionMixedRefl phase(sigma, mean, std_or, std_size,diff_conc,g);
        // SGGXPhaseFunctionNormalDistribution phase(sigma,mean,std_or);

        // phase.SetMicroflakes(dir, sigma);//Fiber-like distributions

        // Diffuse phase function
        //  Point2f u1{r(), r()}, u2{r(), r()};
        bool is_specular_phase = use_specular;
        
        for (int depth = 0; depth < maxDepth; ++depth) {
            // Follow random walk through layers to sample layered BSDF
            // Possibly terminate layered BSDF sampling with Russian Roulette
            Float rrBeta = f.MaxComponentValue() / pdf;
            if (depth > 128 && rrBeta < 0.01f) {
            Float q = std::max<Float>(0, 1 - rrBeta);
            if (r() < q)
            return {};
            pdf *= 1 - q;
            }
            if (w.z == 0)
                return {};

            if (albedo_diff || albedo_flakes) {
                // Sample potential scattering event in layered medium
                Float pa = phase.projectedArea(w, Point2f(r(), r()));
                
                Float reduced_st = diff_conc * sigma_t + (1-diff_conc)*(sigma_t * pa) ;

                Float dz = SampleExponential(r(), reduced_st / AbsCosTheta(w));
                Float zp = w.z > 0 ? (z + dz) : (z - dz);
                CHECK_RARE(1e-5, zp == z);
                if (zp == z)
                    return {};
                if (0 < zp && zp < thickness) {
                    // std::cout<<"should not happen "<<std::endl;
                    // Update path state for valid scattering event between interfaces

                    pstd::optional<PhaseFunctionSample> ps;
                    //get the albedo of the sampled particle
                    SampledSpectrum albedo_particle;
                    ps = phase.Sample_p(-w, Point2f(r(), r()),Point2f(r(), r()),r(),Point2f(r(), r())); 
                    if (ps->p_flag == _MediumFlags::Diffuser){
                                
                                albedo_particle = albedo_diff;
                            }else{
                                albedo_particle = albedo_flakes * pa;
                    }
                    // if (!albedo_particle){
                        // LOG_ERROR("Albedo particle is not instantiated");
                    // } 
                    if (!ps || ps->pdf == 0 || ps->wi.z == 0){
                        // LOG_ERROR("early termination is not great",pPixel.x, pPixel.y, sampleIndex);
                        return {};

                    }
                    auto old_f = f;

                    f *= albedo_particle;

                    
                    specularPath = false;
                    w = ps->wi;
                    z = zp;

                    continue;
                }
                z = Clamp(zp, 0, thickness);
                if (z == 0)
                    DCHECK_LT(w.z, 0);
                else
                    DCHECK_GT(w.z, 0);

            } else {
                // Advance to the other layer interface
                z = (z == thickness) ? 0 : thickness;
                f *= Tr(thickness, w);
            }
            // Initialize _interface_ for current interface surface
#ifdef interface  // That's enough out of you, Windows.
#undef interface
#endif
            // Old Code
            TopOrBottomBxDF<TopBxDF, BottomBxDF> interface;
            // BottomBxDF interface;
            // std::cout << "z is "<< z << std::endl;
            
            if ( z == 0){
                // std::cout<<"bot depth "<< depth << " z "<< z << " thickness "<<thickness <<std::endl;
                interface = &bottom;
            }
            else{
                interface = &top;

            }

            


            // Sample interface BSDF to determine new path direction
            Float uc = r();
            Point2f u(r(), r());
            // return bottom.Sample_f(-wo, uc, u, mode);
            //old Code
            //pstd::optional<BSDFSample> bs = interface.Sample_f(-w, uc, u, mode);
            
            // std::cout<<"w "<< w << std::endl ;
            pstd::optional<BSDFSample> bs = interface.Sample_f(-w, uc, u, mode);
            if (!bs || !bs->f || bs->pdf == 0 || bs->wi.z == 0){
                        // LOG_ERROR("early termination is not great");
                return {};

            }
            


            // Old Code 
            f *= bs->f;
            // f = bs->f;
            if (IsInf(f[0])){
                    std::cout<< "it's INF line 1990"<< std::endl;
            }
            //Old code
            pdf *= bs->pdf;
            // pdf = bs->pdf;
            specularPath &= bs->IsSpecular();
            w = bs->wi;
            
            // bool specularPath2 = bs->IsSpecular();
            
            // BxDFFlags flags2 = SameHemisphere(wo, bs->wi) ? BxDFFlags::Reflection : BxDFFlags::Transmission;
            
            // flags2 |= specularPath2 ? BxDFFlags::Specular : BxDFFlags::Glossy;
            // if (flipWi){
            // return BSDFSample(f, -w,pdf, flags2);
            // }else{
            // return BSDFSample(f, w, pdf, flags2);
            // }
            //Here lies the problem! If the thickness is zero it will always sample the bot interface meanwhile it should also sample the top one, allowing the reflected ray to be transmitted
            // Since this accept only transmission events it will be always closed
            // Return _BSDFSample_ if path has left the layers
            if (bs->IsTransmission() ) {
                BxDFFlags flags = SameHemisphere(wo, w) ? BxDFFlags::Reflection
                                                        : BxDFFlags::Transmission;
                flags |= specularPath ? BxDFFlags::Specular : BxDFFlags::Glossy;
                if (flipWi)
                    w = -w;
                if (f.HasNaNs()) {
                    std::cout << " nan found in sample()" << std::endl;
                }
                if (IsInf(f[0])){
                    std::cout<< "it's INF"<< std::endl;
                }
                // std::cout<<"exiting f "<< f << " w " << w <<  " depth "<< depth << std::endl; 
                if ( z == thickness){

                    //scatter towards the air
                    return BSDFSample(f, w, pdf, flags,1.0f);
                }else{
                    //scatter inside the skin
                return BSDFSample(f, w, pdf, flags,1.3);

                }
                // Old Code
                // return BSDFSample(f, w, pdf, flags, 1.0f, true);
            }
            // Scale _f_ by cosine term after scattering at the interface
            f *= AbsCosTheta(bs->wi);
        }
        // LOG_ERROR("Reaching End");
        return {};
    }

    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
        CHECK(sampleFlags == BxDFReflTransFlags::All);  // for now
        // Set _wo_ and _wi_ for layered BSDF evaluation
        bool enteredTop = twoSided || wo.z > 0;

        //****************************************************
        //                   Debug code
        //****************************************************
        
        //  if(!enteredTop){
            // return top.PDF(wo, wi, mode, sampleFlags);
            // wo = -wo;
            // if (!(sampleFlags & BxDFReflTransFlags::Reflection))
            //     return 0;
            // return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;

        //  }

        //****************************************************
        //                   END Debug code
        //****************************************************

        // return bottom.PDF(wo,wi,mode,sampleFlags);
        
        if (twoSided && wo.z < 0) {
            wo = -wo;
            wi = -wi;
        }
                
        // Declare _RNG_ for layered PDF evaluation
        RNG rng(Hash(GetOptions().seed, wi), Hash(wo));
        auto r = [&rng]() {
            return std::min<Float>(rng.Uniform<Float>(), OneMinusEpsilon);
        };

        // Update _pdfSum_ for reflection at the entrance layer
        Float pdfSum = 0;
        if (SameHemisphere(wo, wi)) {
            auto reflFlag = BxDFReflTransFlags::Reflection;
            pdfSum += enteredTop ? nSamples * top.PDF(wo, wi, mode, reflFlag)
                                 : nSamples * bottom.PDF(wo, wi, mode, reflFlag);
        }

        for (int s = 0; s < nSamples; ++s) {
            // Evaluate layered BSDF PDF sample
            if (SameHemisphere(wo, wi)) {
                // Evaluate TRT term for PDF estimate
                // Old Code
                TopOrBottomBxDF<TopBxDF, BottomBxDF> rInterface, tInterface;
                // const TopBxDF *rInterface, *tInterface;
                // TopBxDF rInterface, tInterface;
                
                if (enteredTop) {
                    rInterface = &bottom;
                    tInterface = &top;
                } else {
                    rInterface = &top;
                    tInterface = &bottom;
                }
                // Sample _tInterface_ to get direction into the layers
                auto trans = BxDFReflTransFlags::Transmission;
                pstd::optional<BSDFSample> wos, wis;
                wos = tInterface.Sample_f(wo, r(), {r(), r()}, mode, trans);
                wis = tInterface.Sample_f(wi, r(), {r(), r()}, !mode, trans);

                // Update _pdfSum_ accounting for TRT scattering events
                if (wos && wos->f && wos->pdf > 0 && wis && wis->f && wis->pdf > 0) {
                    if (!IsNonSpecular(tInterface.Flags()))
                        pdfSum += rInterface.PDF(-wos->wi, -wis->wi, mode);
                    else {
                        // Use multiple importance sampling to estimate PDF product
                        pstd::optional<BSDFSample> rs =
                            rInterface.Sample_f(-wos->wi, r(), {r(), r()}, mode);
                        if (rs && rs->f && rs->pdf > 0) {
                            if (!IsNonSpecular(rInterface.Flags()))
                                pdfSum += tInterface.PDF(-rs->wi, wi, mode);
                            else {
                                // Compute MIS-weighted estimate of Equation
                                // (\ref{eq:pdf-triple-canceled-one})
                                Float rPDF = rInterface.PDF(-wos->wi, -wis->wi, mode);
                                Float wt = PowerHeuristic(1, wis->pdf, 1, rPDF);
                                pdfSum += wt * rPDF;

                                Float tPDF = tInterface.PDF(-rs->wi, wi, mode);
                                wt = PowerHeuristic(1, rs->pdf, 1, tPDF);
                                pdfSum += wt * tPDF;
                            }
                        }
                    }
                }

            } else {
                
                // Evaluate TT term for PDF estimate
                TopOrBottomBxDF<TopBxDF, BottomBxDF> toInterface, tiInterface;
                // TopBxDF toInterface, tiInterface;
                
                if (enteredTop) {
                    toInterface = &top;
                    tiInterface = &bottom;
                } else {
                    toInterface = &bottom;
                    tiInterface = &top;
                }

                Float uc = r();
                Point2f u(r(), r());
                
                pstd::optional<BSDFSample> wos = toInterface.Sample_f(wo, uc, u, mode);
                if (!wos || !wos->f || wos->pdf == 0 || wos->wi.z == 0 ||
                    wos->IsReflection())
                    continue;

                uc = r();
                u = Point2f(r(), r());
                pstd::optional<BSDFSample> wis = tiInterface.Sample_f(wi, uc, u, !mode);
                if (!wis || !wis->f || wis->pdf == 0 || wis->wi.z == 0 ||
                    wis->IsReflection())
                    continue;

                if (IsSpecular(toInterface.Flags()))
                    pdfSum += tiInterface.PDF(-wos->wi, wi, mode);
                else if (IsSpecular(tiInterface.Flags()))
                    pdfSum += toInterface.PDF(wo, -wis->wi, mode);
                else
                    pdfSum += (toInterface.PDF(wo, -wis->wi, mode) +
                               tiInterface.PDF(-wos->wi, wi, mode)) /
                              2;
            }
        }
        // Return mixture of PDF estimate and constant PDF
        return Lerp(0.9f, 1 / (4 * Pi), pdfSum / nSamples);
    }

  private:
    // LayeredBxDFMixedParticlesNormalizedFresnel Private Methods
    PBRT_CPU_GPU
    static Float Tr(Float dz, Vector3f w) {
        if (std::abs(dz) <= std::numeric_limits<Float>::min())
            return 1;
        return FastExp(-std::abs(dz / w.z));
    }
    static Float Tr(Float dz, Vector3f w, Float sigmaTAlongWi) {
        return FastExp(-std::abs(sigmaTAlongWi * dz / w.z));
    }

    // LayeredBxDFMixedParticlesNormalizedFresnel Private Members
    TopBxDF top;
    BottomBxDF bottom;
    Float thickness, g1,g2,w_g, sigma, mean, std_or, std_size, diff_conc;
    Float sigma_t,eta;

    Vector3f dir;
    
    SampledSpectrum albedo_diff,albedo_flakes;
    int maxDepth, nSamples;
    bool use_specular;
};



// CoatedDiffuseBxDF Definition
class CoatedDiffuseBxDF : public LayeredBxDF<DielectricBxDF, DiffuseBxDF, true> {
  public:
    // CoatedDiffuseBxDF Public Methods
    using LayeredBxDF::LayeredBxDF;
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "CoatedDiffuseBxDF"; }
};

class MyCoatedDiffuseBxDF : public MyLayeredBxDF<DielectricBxDF, DiffuseBxDF, true> {
  public:
    // CoatedDiffuseBxDF Public Methods
    using MyLayeredBxDF::MyLayeredBxDF;
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "MyCoatedDiffuseBxDF"; }
};
class CoatedDiffuseBxDFMixedParticles
    : public LayeredBxDFMixedParticles<DielectricBxDF, DiffuseBxDF, true> {
  public:
    // CoatedDiffuseBxDF Public Methods
    using LayeredBxDFMixedParticles::LayeredBxDFMixedParticles;
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "CoatedDiffuseBxDFMixedParticles"; }
};
class TransmittanceBxDFMixedParticles: public LayeredBxDFMixedParticlesTransmittance<DielectricBxDF, DielectricBxDF, false> {
  public:
    // Transmittance Public Methods
    using LayeredBxDFMixedParticlesTransmittance::LayeredBxDFMixedParticlesTransmittance;
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "Transmittance"; }
};

class TransmittanceBxDFMixedParticlesAdvanced: public LayeredBxDFMixedParticlesTransmittanceAdvanced<DielectricBxDF, DielectricBxDF, false> {
  public:
    // Transmittance Public Methods
    using LayeredBxDFMixedParticlesTransmittanceAdvanced::LayeredBxDFMixedParticlesTransmittanceAdvanced;
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "Transmittance"; }
};


class TransmittanceBxDFNormFresnel: public LayeredBxDFMixedParticlesNormalizedFresnel<DielectricBxDF, NormalizedFresnelBxDF, false> {
  public:
    // Transmittance Public Methods
    using LayeredBxDFMixedParticlesNormalizedFresnel::LayeredBxDFMixedParticlesNormalizedFresnel;
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "LayeredBxDFMixedParticlesNormalizedFresnel"; }
};

class TransmittanceBxDFNormFresnelAdvanced: public LayeredBxDFMixedParticlesNormalizedFresnelAdvanced<DielectricBxDF, NormalizedFresnelBxDF, false> {
  public:
    // Transmittance Public Methods
    using LayeredBxDFMixedParticlesNormalizedFresnelAdvanced::LayeredBxDFMixedParticlesNormalizedFresnelAdvanced;
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "LayeredBxDFMixedParticlesNormalizedFresnelAdvanced"; }
};


// CoatedConductorBxDF Definition
class CoatedConductorBxDF : public LayeredBxDF<DielectricBxDF, ConductorBxDF, true> {
  public:
    // CoatedConductorBxDF Public Methods
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "CoatedConductorBxDF"; }
    using LayeredBxDF::LayeredBxDF;
};

// HairBxDF Definition
class HairBxDF {
  public:
    // HairBxDF Public Methods
    HairBxDF() = default;
    PBRT_CPU_GPU
    HairBxDF(Float h, Float eta, const SampledSpectrum &sigma_a, Float beta_m,
             Float beta_n, Float alpha);
    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const;
    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(Vector3f wo, Float uc, Point2f u,
                                        TransportMode mode,
                                        BxDFReflTransFlags sampleFlags) const;
    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags) const;

    PBRT_CPU_GPU
    void Regularize() {}

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "HairBxDF"; }
    std::string ToString() const;

    PBRT_CPU_GPU
    BxDFFlags Flags() const { return BxDFFlags::GlossyReflection; }

    PBRT_CPU_GPU
    static RGBUnboundedSpectrum SigmaAFromConcentration(Float ce, Float cp);
    PBRT_CPU_GPU
    static SampledSpectrum SigmaAFromReflectance(const SampledSpectrum &c, Float beta_n,
                                                 const SampledWavelengths &lambda);

  private:
    // HairBxDF Constants
    static constexpr int pMax = 3;

    // HairBxDF Private Methods
    PBRT_CPU_GPU static Float Mp(Float cosTheta_i, Float cosTheta_o, Float sinTheta_i,
                                 Float sinTheta_o, Float v) {
        Float a = cosTheta_i * cosTheta_o / v, b = sinTheta_i * sinTheta_o / v;
        Float mp = (v <= .1)
                       ? (FastExp(LogI0(a) - b - 1 / v + 0.6931f + std::log(1 / (2 * v))))
                       : (FastExp(-b) * I0(a)) / (std::sinh(1 / v) * 2 * v);
        DCHECK(!IsInf(mp) && !IsNaN(mp));
        return mp;
    }

    PBRT_CPU_GPU static pstd::array<SampledSpectrum, pMax + 1> Ap(Float cosTheta_o,
                                                                  Float eta, Float h,
                                                                  SampledSpectrum T) {
        pstd::array<SampledSpectrum, pMax + 1> ap;
        // Compute $p=0$ attenuation at initial cylinder intersection
        Float cosGamma_o = SafeSqrt(1 - Sqr(h));
        Float cosTheta = cosTheta_o * cosGamma_o;
        Float f = FrDielectric(cosTheta, eta);
        ap[0] = SampledSpectrum(f);

        // Compute $p=1$ attenuation term
        ap[1] = Sqr(1 - f) * T;

        // Compute attenuation terms up to $p=_pMax_$
        for (int p = 2; p < pMax; ++p)
            ap[p] = ap[p - 1] * T * f;

        // Compute attenuation term accounting for remaining orders of scattering
        if (1 - T * f)
            ap[pMax] = ap[pMax - 1] * f * T / (1 - T * f);

        return ap;
    }

    PBRT_CPU_GPU static inline Float Phi(int p, Float gamma_o, Float gamma_t) {
        return 2 * p * gamma_t - 2 * gamma_o + p * Pi;
    }

    PBRT_CPU_GPU static inline Float Np(Float phi, int p, Float s, Float gamma_o,
                                        Float gamma_t) {
        Float dphi = phi - Phi(p, gamma_o, gamma_t);
        // Remap _dphi_ to $[-\pi,\pi]$
        while (dphi > Pi)
            dphi -= 2 * Pi;
        while (dphi < -Pi)
            dphi += 2 * Pi;

        return TrimmedLogistic(dphi, s, -Pi, Pi);
    }

    PBRT_CPU_GPU
    pstd::array<Float, pMax + 1> ApPDF(Float cosTheta_o) const;

    // HairBxDF Private Members
    Float h, eta;
    SampledSpectrum sigma_a;
    Float beta_m, beta_n;
    Float v[pMax + 1];
    Float s;
    Float sin2kAlpha[pMax], cos2kAlpha[pMax];
};

// MeasuredBxDF Definition
class MeasuredBxDF {
  public:
    // MeasuredBxDF Public Methods
    MeasuredBxDF() = default;
    PBRT_CPU_GPU
    MeasuredBxDF(const MeasuredBxDFData *brdf, const SampledWavelengths &lambda)
        : brdf(brdf), lambda(lambda) {}

    static MeasuredBxDFData *BRDFDataFromFile(const std::string &filename,
                                              Allocator alloc);

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const;

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(Vector3f wo, Float uc, Point2f u,
                                        TransportMode mode,
                                        BxDFReflTransFlags sampleFlags) const;
    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags) const;

    PBRT_CPU_GPU
    void Regularize() {}

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "MeasuredBxDF"; }

    std::string ToString() const;

    PBRT_CPU_GPU
    BxDFFlags Flags() const { return (BxDFFlags::Reflection | BxDFFlags::Glossy); }

  private:
    // MeasuredBxDF Private Methods
    PBRT_CPU_GPU
    static Float theta2u(Float theta) { return std::sqrt(theta * (2 / Pi)); }
    PBRT_CPU_GPU
    static Float phi2u(Float phi) { return phi * (1 / (2 * Pi)) + .5f; }

    PBRT_CPU_GPU
    static Float u2theta(Float u) { return Sqr(u) * (Pi / 2.f); }
    PBRT_CPU_GPU
    static Float u2phi(Float u) { return (2.f * u - 1.f) * Pi; }

    // MeasuredBxDF Private Members
    const MeasuredBxDFData *brdf;
    SampledWavelengths lambda;
};

class MyMeasuredBxDF {
  public:
    // MyMeasuredBxDF Public Methods
    MyMeasuredBxDF() = default;
    PBRT_CPU_GPU
    MyMeasuredBxDF(const MyMeasuredBxDFData *brdf, const SampledWavelengths &lambda)
        : brdf(brdf), lambda(lambda) {theta_o_res = 90;
    theta_i_res = 90;
    wl_res = 80;
    }

    static MyMeasuredBxDFData *BRDFDataFromFile(const std::string &filename,
                                              Allocator alloc);

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const;

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(Vector3f wo, Float uc, Point2f u,
                                        TransportMode mode,
                                        BxDFReflTransFlags sampleFlags) const;
    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags) const;

    PBRT_CPU_GPU
    void Regularize() {}

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "MyMeasuredBxDF"; }

    std::string ToString() const;

    PBRT_CPU_GPU
    BxDFFlags Flags() const { return (BxDFFlags::Reflection | BxDFFlags::Glossy); }

  private:
    // MyMeasuredBxDF Private Methods
    PBRT_CPU_GPU
    static Float theta2u(Float theta) { return std::sqrt(theta * (2 / Pi)); }
    PBRT_CPU_GPU
    static Float phi2u(Float phi) { return phi * (1 / (2 * Pi)) + .5f; }

    PBRT_CPU_GPU
    static Float u2theta(Float u) { return Sqr(u) * (Pi / 2.f); }
    PBRT_CPU_GPU
    static Float u2phi(Float u) { return (2.f * u - 1.f) * Pi; }
    size_t  theta_o_res;
    size_t  theta_i_res;
    size_t wl_res;
    // MyMeasuredBxDF Private Members
    const MyMeasuredBxDFData *brdf;
    SampledWavelengths lambda;
};


// SheenVolumeBxDF Definition (Juan Raul testing stuff!!!!)
// Adaption of Tizian implementation: Practical Multiple-Scattering Sheen Using Linearly
// Transformed Cosines (https://github.com/tizian/ltc-sheen/tree/master)
class SheenVolumeBxDF {
  public:
    SheenVolumeBxDF() = default;
    PBRT_CPU_GPU
    SheenVolumeBxDF(int maxBounces, SampledSpectrum albedo, Float density, Float sigma,
                    uint64_t seed)
        : maxBounces(maxBounces),
          albedo(albedo),
          density(density),
          sigma(sigma),
          seed(seed) {
        // logSGGX();
    }

    PBRT_CPU_GPU
    SampledSpectrum f(Vector3f wo, Vector3f wi, TransportMode mode) const;

    PBRT_CPU_GPU
    pstd::optional<BSDFSample> Sample_f(Vector3f wo, Float uc, Point2f u,
                                        TransportMode mode,
                                        BxDFReflTransFlags sampleFlags) const;
    PBRT_CPU_GPU
    Float PDF(Vector3f wo, Vector3f wi, TransportMode mode,
              BxDFReflTransFlags sampleFlags) const;

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "SheenVolumeBxDF"; }

    std::string ToString() const;

    PBRT_CPU_GPU
    void Regularize() {}

    PBRT_CPU_GPU
    BxDFFlags Flags() const { return (BxDFFlags::Reflection | BxDFFlags::Glossy); }

  private:
    int maxBounces;
    SampledSpectrum albedo;
    Float density, sigma;
    uint64_t seed;
};

inline SampledSpectrum BxDF::f(Vector3f wo, Vector3f wi, TransportMode mode) const {
    auto f = [&](auto ptr) -> SampledSpectrum { return ptr->f(wo, wi, mode); };
    return Dispatch(f);
}

inline pstd::optional<BSDFSample> BxDF::Sample_f(Vector3f wo, Float uc, Point2f u,
                                                 TransportMode mode,
                                                 BxDFReflTransFlags sampleFlags) const {
    auto sample_f = [&](auto ptr) -> pstd::optional<BSDFSample> {
        return ptr->Sample_f(wo, uc, u, mode, sampleFlags);
    };
    return Dispatch(sample_f);
}

inline Float BxDF::PDF(Vector3f wo, Vector3f wi, TransportMode mode,
                       BxDFReflTransFlags sampleFlags) const {
    auto pdf = [&](auto ptr) { return ptr->PDF(wo, wi, mode, sampleFlags); };
    return Dispatch(pdf);
}

inline BxDFFlags BxDF::Flags() const {
    auto flags = [&](auto ptr) { return ptr->Flags(); };
    return Dispatch(flags);
}

inline void BxDF::Regularize() {
    auto regularize = [&](auto ptr) { ptr->Regularize(); };
    return Dispatch(regularize);
}

extern template class LayeredBxDF<DielectricBxDF, DiffuseBxDF, true>;
extern template class LayeredBxDF<DielectricBxDF, ConductorBxDF, true>;

extern template class MyLayeredBxDF<DielectricBxDF, DiffuseBxDF, true>;
extern template class LayeredBxDFMixedParticles<DielectricBxDF, DiffuseBxDF, true>;
extern template class LayeredBxDFMixedParticlesTransmittance<DielectricBxDF, DielectricBxDF, false>;
extern template class LayeredBxDFMixedParticlesTransmittanceAdvanced<DielectricBxDF, DielectricBxDF, false>;
extern template class LayeredBxDFMixedParticlesNormalizedFresnel<DielectricBxDF, NormalizedFresnelBxDF, false>;
extern template class LayeredBxDFMixedParticlesNormalizedFresnelAdvanced<DielectricBxDF, NormalizedFresnelBxDF, false>;


}  // namespace pbrt

#endif  // PBRT_BXDFS_H
