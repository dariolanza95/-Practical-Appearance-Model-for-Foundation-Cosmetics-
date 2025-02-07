// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

// Unit tests to virtually measure the reflectance (eval), pdf and samples of materials
// Remember to include these file inside CMakeLists.txt (add
// src/pbrt/measurements_test.cpp)

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
#include <pbrt/util/math.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>

#include <json/json.h>
#include <cmath>
#include <numeric>

using namespace pbrt;

/*
inline SampledSpectrum Pow(const SampledSpectrum &s, Float e) {
    SampledSpectrum ret;
    for (int i = 0; i < NSpectrumSamples; ++i)
        ret[i] = std::pow(s[i], e);
    return ret;
}
*/

// Conversion from spherical coordinates to unit directions.
inline Vector3f sphericalDirection(Float theta, Float phi) {
    double sin_theta = std::sin(theta), cos_theta = std::cos(theta),
           sin_phi = std::sin(phi), cos_phi = std::cos(phi);
    return Vector3f(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
}

// Conversion from unit directions to spherical coordinates.
inline std::pair<Float, Float> sphericalCoordinates(const Vector3f& v) {
    Float theta = std::acos(v.z), phi = std::atan2(v.y, v.x);
    if (phi < 0.0) {
        phi += 2.0 * Pi;
    }
    return {theta, phi};
}

// Analogous to np.linspace
template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N - 1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

void EvalBsdf(
    std::function<BSDF*(const SurfaceInteraction&, Allocator, Float*)> createBSDF,
    const char* description, std::string out_path, Float* extra_params) {
    // const int thetaRes = 20;
    // const int phiRes = 20;
    // const int sampleCount = 20;
    // Float* frequencies = new Float[thetaRes * phiRes];
    // Float* expFrequencies = new Float[thetaRes * phiRes];
    RNG rng;
    Float lambda_min = 380, lambda_max = 830;

    int index = 0;
    std::cout.precision(3);

    // Create BSDF, which requires creating a Shape, casting a Ray that
    // hits the shape to get a SurfaceInteraction object.
    BSDF* bsdf = nullptr;
    auto t = std::make_shared<const Transform>(RotateX(-180));
    auto tInv = std::make_shared<const Transform>(Inverse(*t));

    bool reverseOrientation = true;

    std::shared_ptr<Disk> disk =
        std::make_shared<Disk>(t.get(), tInv.get(), reverseOrientation, 0., 1., 0, 360.);
    Point3f origin(0.1, 0, 1);  // offset slightly so we don't hit center of disk
    Vector3f direction(0, 0, -1);
    Ray r(origin, direction);
    auto si = disk->Intersect(r);
    ASSERT_TRUE(si.has_value());
    bsdf = createBSDF(si->intr, Allocator(), extra_params);

    // --------------------Uniform sampling code ----------------------------
    
    // Original code
    /*
    int theta_i_res = 101;
    int phi_i_res = 1;
    int theta_o_res = 101;
    int phi_o_res = 1;
    */

    int theta_i_res = 181;
    int phi_i_res = 1;
    int theta_o_res = 181;
    int phi_o_res = 1;

    // Angular sampling
    //std::vector<Float> theta_i_values = linspace(0.0f, 100.0f, theta_i_res);
    //std::vector<Float> theta_o_values = linspace(0.0f, 100.0f, theta_o_res);//Original code
    
    std::vector<Float> theta_i_values = linspace(0.0f, 180.0f, theta_i_res);
    std::vector<Float> theta_o_values = linspace(0.0f, 180.0f, theta_o_res);

    std::vector<Float> phi_i_values = linspace(-180.0f, 180.0f, phi_i_res);
    // original code
    //  std::vector<Float> phi_o_values = linspace(-180.0f, 180.0f, phi_o_res);
    std::vector<Float> phi_o_values{-180, 0.0f};
    int n_samples = 1;
    FILE* fp = fopen(out_path.c_str(), "w");
    //    FILE *fp = fopen("../../scripts/measurements/conductor_pbrt.txt", "w");
    Float phi_i = 0.0f;
    for (auto theta_i : theta_i_values) {
        // std::cout<<"theta_i "<< theta_i<<std::endl;
        Vector3f wi = sphericalDirection(Radians(theta_i), Radians(phi_i));
        for (auto phi_o : phi_o_values) {
            for (auto theta_o : theta_o_values) {
                SampledWavelengths lambda = SampledWavelengths::SampleUniform(
                    rng.Uniform<Float>(), lambda_min, lambda_max);
                // Float fake_res(0.8);
                // if (phi_o == Float(0.0f)){
                //     fake_res = 0.2;
                // }
                // fake_res *= theta_o;

                // if (theta_i == Float(30.0f)){
                //     fake_res = 520;
                //     if (phi_o == Float(0.0f)){
                //     fake_res = 832;
                // }
                // }
                SampledSpectrum res(0.0f);
                for (int i = 0; i < n_samples; i++) {
                    // SampledWavelengths lambda =
                    SampledWavelengths::SampleVisible(
                        rng.Uniform<Float>());  // Random spectral samples
                    // wo is the Observation Vector.
                    Vector3f wo = sphericalDirection(Radians(theta_o), Radians(phi_o));

                    // Evaluate BSDF
                    // Fixed Camera, moving light
                    SampledSpectrum f = bsdf->f(wo, wi, TransportMode::Radiance) * Pi;  // *AbsCosTheta(wi);
                    // debugging purposes
                    // SampledSpectrum f(fake_res/100);
                    res += f;

                    // Fixed Light, moving camera
                    // SampledSpectrum f = bsdf->f(wi, wo, TransportMode::Radiance)
                    // *AbsCosTheta(wo);

                    // Old code
                    // SampledSpectrum f = bsdf.f(wo, wi, TransportMode::Radiance)
                    // *AbsCosTheta(wo);
                    //  SampledSpectrum f = bsdf.f(wi, wo, TransportMode::Radiance)
                    //  *AbsCosTheta(wo);

                    // std::cout << lambda[0] << std::endl;
                    // std::cout << f[0] << std::endl;

                    // Save reflectance (1 entry per wavelength)
                }

                res /= (float)n_samples;

                //Diagonal Reflectance sanity check
                
                /*
                if (int(theta_i) == int(theta_o)){
                    
                    //std::cout << "reflectance: " << theta_i << ", " << theta_o << std::endl;

                    res = SampledSpectrum(100.0f);
                }else{
                    res = SampledSpectrum(0.0f);
                }
                */

                //Gamma correction
                Float gamma = 1.0/2.2;
                //Float gamma = 2.2;
                //res = Pow(res, gamma);

                // char rgb[3];
                // rgb[0] = 'r';
                // rgb[1] = 'g';
                // rgb[2] = 'b';
                /*
                for (int k = 0; k < 3; k++) {
                    fprintf(fp, "%f %f %f %f %f %f\n", phi_i, theta_i, phi_o, theta_o,
                            lambda[k], res[k] / n_samples);
                }
                */

                for (int k = 0; k < 3; k++) {
                    fprintf(fp, "%f %f %f %f %f %f\n", phi_i, theta_i, phi_o, theta_o,
                            lambda[k], res[k]);
                }

                //                 for (int k = 0; k < NSpectrumSamples; k++) {
                //     fprintf(fp, "%f %f %f %f %f %f\n", phi_i, theta_i, phi_o, theta_o,
                //             lambda[k], res[k] / n_samples);
                // }
            }
        }
    }

    fclose(fp);
    /*

    // -------------------Matching measurements setup
    // Adaptative sampling for theta_i and theta_o --> we can activate it with a boolean
    // parameter
    std::cout << "out_p" << out_path << std::endl;
    FILE* fp = fopen(out_path.c_str(), "w");
    Float phi_i = 0.0f;
    Float phi_o = -180.0f;
    Float new_theta_o = 0.0f;
    int n_samples = 1;

    // Angles is a dictionary where we store for each theta_i a list of theta_o angles in
    // degrees
    std::ifstream file("../../scripts/measurements/angles.json");
    // std::ifstream file("../../scripts/measurements/angles.json");
    Json::Value root;  // Root of JSON object
    Json::Reader reader;

    bool parsingSuccessful = reader.parse(file, root);
    if (!parsingSuccessful) {
        // Report the failure and their locations in the document.
        std::cout << "Failed to parse JSON: " << reader.getFormattedErrorMessages()
                  << std::endl;
        return;
    }

    for (auto it = root.begin(); it != root.end(); ++it) {
        std::string name = it.key().asString();
        Float theta_i = std::strtof(name.c_str(), NULL);
        Vector3f wi = sphericalDirection(Radians(theta_i), Radians(phi_i));

        auto value = root[name];

        for (int i = 0; i < value.size(); i++) {
            // std::cout << "theta_o: " << value[i] << '\n';

            std::string aux = value[i].asString();
            Float theta_o = std::strtof(aux.c_str(), NULL);

            // Float phi_o = theta_o > 0 ? -180.0f : 0.0f;

            if (theta_o > 0) {
                // Normal reflection scenario [0, 90]
                phi_o = -180.0f;
                new_theta_o = theta_o;
            } else {
                // Back reflection scenario [-90, 0]
                phi_o = 0.0f;
                // new_theta_o = 90 - theta_o;
                new_theta_o =  - theta_o;
            }

            SampledWavelengths lambda = SampledWavelengths::SampleUniform(
                rng.Uniform<Float>(), lambda_min, lambda_max);

            SampledSpectrum res(0.0f);

            Vector3f wo = sphericalDirection(Radians(new_theta_o), Radians(phi_o));

            for (int i = 0; i < n_samples; i++) {
                // Evaluate BSDF
                // Fixed Camera, moving light
                // SampledSpectrum f = bsdf->f(wo, wi,TransportMode::Radiance) *
                // AbsCosTheta(wi);
                // SampledSpectrum f = bsdf->f(wo, wi,TransportMode::Radiance);  //
                // Measurements are BRDF --> Ignore cosine
                SampledSpectrum f = bsdf->f(wo, wi,TransportMode::Radiance) * Pi ;  //
                                         // Measurements are BRDF --> Ignore cosine
                res += f;
            }

            for (int k = 0; k < NSpectrumSamples; k++) {
                fprintf(fp, "%f %f %f %f %f %f\n", phi_i, theta_i, -180.0f, theta_o,
                        lambda[k], res[k] / n_samples);
            }
        }
    }

    fclose(fp);
*/
    EXPECT_TRUE(true);
}

BSDF* createConductor(const SurfaceInteraction& si, Allocator alloc, float eta, float k,
                      Float roughx) {
    TrowbridgeReitzDistribution distrib(roughx, roughx);
    return alloc.new_object<BSDF>(Normal3f(0, 0, 1), si.shading.dpdu,
                                  alloc.new_object<ConductorBxDF>(
                                      distrib, SampledSpectrum(eta), SampledSpectrum(k)));
}

BSDF* createDiffuse(const SurfaceInteraction& si, Allocator alloc) {
    SampledSpectrum Kd(1.);
    return alloc.new_object<BSDF>(si.shading.n, si.shading.dpdu,
                                  alloc.new_object<DiffuseBxDF>(Kd));
}

BSDF* createDielectric(const SurfaceInteraction& si, Allocator alloc, Float ior,
                       Float roughx) {
    TrowbridgeReitzDistribution distrib(roughx, roughx);
    return alloc.new_object<BSDF>(Normal3f(0, 0, 1), si.shading.dpdu,
                                  alloc.new_object<DielectricBxDF>(ior, distrib));
}

BSDF* createLayeredMaterial(const SurfaceInteraction& si, Allocator alloc, Float mean,
                            Float std_or, Float std_size, Float surf_rough, Float thick,
                            Float sigma, std::vector<Float> albedo_diff,
                            std::vector<Float> albedo_flake, Float ior,
                            Float diff_particles_perc, Float g1, Float g2, Float w_g,
                            bool use_refl, int nSamples = 128) {
    int maxDepth = 32;
    // Float sigma(0.1);
    // Float mean(45.0);
    // Float std(0.0);
    Vector3f dir(0, 0, 1);
    SampledSpectrum r(0.0f);
    // Float roughx(0.5f);
    //  bool use_refl = false;
    TrowbridgeReitzDistribution distrib(surf_rough, surf_rough);
    // BeckmannDistribution distrib(surf_rough, surf_rough);

    // Float thick(0.5f);
    Float sampledEta = ior;
    pstd::span<Float> naffo(albedo_flake);
    SampledSpectrum albedo_flakes(naffo);
    SampledSpectrum albedo_diffs(albedo_diff);
    // MyLayeredBxDF<DielectricBxDF,DiffuseBxDF,true> m_layered(DielectricBxDF(sampledEta,
    // distrib), DiffuseBxDF(r), thick, dir,sigma,mean,std,use_refl,albedo, gg, maxDepth,
    // nSamples);
    std::cout << "nSamples" << nSamples << std::endl;
    return alloc.new_object<BSDF>(
        Normal3f(0, 0, 1), si.shading.dpdu,
        alloc.new_object<CoatedDiffuseBxDFMixedParticles>(
            DielectricBxDF(sampledEta, distrib), DiffuseBxDF(r), thick, dir, sigma, mean,
            std_or, std_size, diff_particles_perc, use_refl, albedo_flakes, albedo_diffs,
            g1, g2, w_g, maxDepth, nSamples));
    // return alloc.new_object<BSDF>(Normal3f(0,0,1),
    // si.shading.dpdu,alloc.new_object<DielectricBxDF>(sampledEta, distrib));
    //  return alloc.new_object<BSDF>(Normal3f(0,0,1),
    //  si.shading.dpdu,alloc.new_object<MyCoatedDiffuseBxDF>(DielectricBxDF(sampledEta,
    //  distrib), DiffuseBxDF(r), thick, dir,sigma,mean,std_or,std_size,use_refl,m_albedo,
    //  gg, maxDepth, nSamples));
}

TEST(Measurements, EvalConductor) {
    std::vector<Float> rough_vals{0.05, 0.1, 0.2, 0.3};
    for (auto rough : rough_vals) {
        std::ostringstream out_path;
        out_path << "../../scripts/measurements/conductor_r_" << std::setprecision(2)
                 << rough << ".txt";
        // std::string out_path  = "../../scripts/measurements/conductor_r_" +
        // std::to_string(rough) + ".txt";
        std::cout << out_path.str() << std::endl;
        EvalBsdf(
            [](const SurfaceInteraction& si, Allocator alloc,
               Float* extra_params) -> BSDF* {
                Float eta(1.0f);
                Float k(1.5f);
                Float rough = extra_params[0];
                return createConductor(si, alloc, eta, k, rough);
            },
            "EvalConductor", out_path.str(), &rough);
    }
}

TEST(Measurements, EvalDiffuse) {
    std::ostringstream out_path;

    Float rough(0.5f);
    out_path << "../../scripts/measurements/layered/diffuse.txt";
    // std::string out_path  = "../../scripts/measurements/conductor_r_" +
    // std::to_string(rough) + ".txt";
    std::cout << out_path.str() << std::endl;
    EvalBsdf([](const SurfaceInteraction& si, Allocator alloc,
                Float* extra_params) -> BSDF* { return createDiffuse(si, alloc); },
             "EvalDiffuse", out_path.str(), nullptr);
}

TEST(Measurements, EvalDielectric) {
    std::vector<Float> rough_vals{0.05, 0.1, 0.2, 0.3};
    std::vector<Float> IoRs{1.0, 1.3, 1.5};
    for (auto eta : IoRs) {
        for (auto rough : rough_vals) {
            std::ostringstream out_path;
            out_path << "../../scripts/measurements/dielectric_r_" << std::setprecision(2)
                     << rough << "_eta_" << eta << ".txt";
            // std::string out_path  = "../../scripts/measurements/conductor_r_" +
            // std::to_string(rough) + ".txt";
            Float extraParams[2];
            extraParams[0] = eta;
            extraParams[1] = rough;
            std::cout << out_path.str() << std::endl;
            EvalBsdf(
                [](const SurfaceInteraction& si, Allocator alloc,
                   Float* extra_params) -> BSDF* {
                    Float eta(0.0f);
                    Float ior = extra_params[0];
                    Float rough = extra_params[1];
                    return createDielectric(si, alloc, ior, rough);
                },

                "EvalConductor", out_path.str(), extraParams);
        }
    }
}
std::vector<Float> json_to_vector(const Json::Value input_js) {
    std::vector<Float> res;
    for (int i = 0; i < input_js.size(); ++i) {
        auto val = input_js[i].asFloat();
        res.push_back(val);
    }
    return res;
}

TEST(Measurements, EvalLayered) {
    /*
    std::vector<Float> means{0,5,10};
    std::vector<Float> stds_or{0.0};
    std::vector<Float> stds_size{0.0};
    // std::vector<Float> stds_size{0.1,0.5,1.5};
    std::vector<Float> thickness{1.0,2.0,5.0};
    // std::vector<Float> thickness{0.5,1.5,2,4};
    std::vector<Float> sigmas{0.01,0.02,0.04,0.06,0.08};
    // std::vector<Float> sigmas{0.05,0.1,0.2};
    std::vector<Float> refl_levels{0};
    std::vector<Float> surf_roughs{0.04,0.05,0.06,0.07,0.08,0.09,0.1};
    // std::vector<Float> albedos{0.5,0.6,0.7,0.8,0.9};
    */
    std::ifstream file("../../scripts/measurements/params_range.json");
    // std::ifstream file("../../scripts/measurements/bluffo.json");
    Json::Value root;  // Root of JSON object
    Json::Reader reader;

    bool parsingSuccessful = reader.parse(file, root);
    if (!parsingSuccessful) {
        // Report the failure and their locations in the document.
        std::cout << "Failed to parse JSON: " << reader.getFormattedErrorMessages()
                  << std::endl;
        return;
    }

    /*const Json::Value albedo_js = root["surf_rough"];
    const Json::Value albedo_js = root["thickness"];
    const Json::Value albedo_js = root["refl_levels"];

    const Json::Value albedo_js = root["mean_or"];
    const Json::Value albedo_js = root["std_or"];
    const Json::Value albedo_js = root["mean_size"];
    const Json::Value albedo_js = root["std_size"];*/

    // Iterate through the array and print its elements
    // std::vector<Float> albedos,surf_roughs,thickness,refl_levels;
    // std::vector<Float> means_or,stds_or,mean_size,stds_size;

    // std::vector<Float> albedo_flakes = json_to_vector(root["albedo_flake"]);
    //    std::vector<Float> albedo_flakes_g = json_to_vector(root["albedo_flake_g"]);
    //    std::vector<Float> albedo_flakes_b = json_to_vector(root["albedo_flake_b"]);
    // std::vector<Float> albedo_diffs = json_to_vector(root["albedo_diff"]);
    std::cout<<"Reading params_range.json"<<std::endl;
    std::vector<Float> surf_roughs = json_to_vector(root["surf_roughs"]);
    std::vector<Float> thickness = json_to_vector(root["thickness"]);
    std::vector<Float> refl_levels = json_to_vector(root["refl_levels"]);
    std::vector<Float> means_or = json_to_vector(root["means_or"]);
    std::vector<Float> stds_or = json_to_vector(root["stds_or"]);
    std::vector<Float> mean_size = json_to_vector(root["mean_size"]);
    std::vector<Float> stds_size = json_to_vector(root["stds_size"]);
    std::vector<Float> iors = json_to_vector(root["ior"]);
    std::vector<Float> diff_particles_percs =
        json_to_vector(root["diff_particles_percs"]);
    std::vector<Float> g1s = json_to_vector(root["g1"]);
    std::vector<Float> g2s = json_to_vector(root["g2"]);
    std::vector<Float> wgs = json_to_vector(root["w_g"]);
    std::vector<Float> n_samples = json_to_vector(root["n_samples"]);
    if (n_samples.empty()) {
        n_samples.push_back(8.f);
    }
    std::vector<Float> albedo_diff = json_to_vector(root["albedo_diff"]);
    // achromatic albedo
    if (albedo_diff.size() == 1) {
        std::cout << "Achromatic diffusers" << std::endl;

        albedo_diff.push_back(albedo_diff[0]);
        albedo_diff.push_back(albedo_diff[0]);
    }
    std::cout<<"albedo_diff r "<< albedo_diff[0] << " g " << albedo_diff[1] << " b " << albedo_diff[2]<<  std::endl;
    std::vector<Float> albedo_flake = json_to_vector(root["albedo_flake"]);
    // achromatic albedo
    if (albedo_flake.size() == 1) {
        std::cout << "Achromatic flakes" << std::endl;
        albedo_flake.push_back(albedo_flake[0]);
        albedo_flake.push_back(albedo_flake[0]);
    }

    std::cout << "Measure test" << std::endl;

    // Access values in JSON
    // std::vector<std::vector<Float>>
    for (auto mean : means_or) {
        for (auto std_or : stds_or) {
            for (auto thick : thickness) {
                for (auto sigma : mean_size) {
                    for (auto refl : refl_levels) {
                        for (auto std_size : stds_size) {
                            for (auto surf_rough : surf_roughs) {
                                for (auto ior : iors) {
                                    for (auto diff_particles_perc :
                                         diff_particles_percs) {
                                        for (auto g1 : g1s) {
                                            for (auto g2 : g2s) {
                                                for (auto w_g : wgs) {
                                                    std::ostringstream out_path;
                                                    auto refl_type =
                                                        refl == 0 ? "_spec" : "_diff";
                                                    out_path
                                                        << "../../scripts/"
                                                           "measurements/layered/"
                                                           "layered_mean_"
                                                        << std::setprecision(4) << mean
                                                        << "_std_or_" << std_or
                                                        << "_std_sz_" << std_size << "_t_"
                                                        << thick << refl_type << "_s_"
                                                        << sigma << "_r_" << surf_rough
                                                        << "_af_r" << albedo_flake[0]<< "_g_"<< albedo_flake[1] << "_b_"<<albedo_flake[2]
                                                        << "_ad_r" << albedo_diff[0]<< "_g_"<< albedo_diff[1] << "_b_"<<albedo_diff[2]
                                                        << "_ior_" << ior << "_conc_"
                                                        << diff_particles_perc << "_g1_"
                                                        << g1 << "_g2_" << g2 << "_wg_"
                                                        << w_g << ".txt";
                                                    // layered_mean_60_std_0_diffuse_0.15
                                                    // std::string out_path  =
                                                    // "../../scripts/measurements/conductor_r_"
                                                    // + std::to_string(rough) +
                                                    // ".txt";
                                                    Float extraParams[21];
                                                    extraParams[0] = mean;
                                                    extraParams[1] = std_or;
                                                    extraParams[2] = surf_rough;  // Interface's
                                                                     // roughness
                                                    extraParams[3] = thick;
                                                    extraParams[4] = sigma;
                                                    extraParams[5] = albedo_diff[0];  // albedo
                                                    extraParams[6] = albedo_flake[0];  // albedo
                                                    extraParams[7] = refl;  // Specular (0) or
                                                               // Diffuse (1)
                                                               // Reflectance?
                                                    extraParams[8] = std_size;
                                                    extraParams[9] = ior;
                                                    extraParams[10] = diff_particles_perc;
                                                    extraParams[11] = g1;
                                                    extraParams[12] = g2;
                                                    extraParams[13] = w_g;
                                                    extraParams[14] = n_samples[0];

                                                    extraParams[15] = albedo_diff[0];
                                                    extraParams[16] = albedo_diff[1];
                                                    extraParams[17] = albedo_diff[2];

                                                    extraParams[18] = albedo_flake[0];
                                                    extraParams[19] = albedo_flake[1];
                                                    extraParams[20] = albedo_flake[2];

                                                    std::cout << out_path.str()
                                                              << std::endl;
                                                    std::cout << "refl is " << refl
                                                              << std::endl;
                                                    EvalBsdf(
                                                        [](const SurfaceInteraction& si,
                                                           Allocator alloc,
                                                           Float* extra_params) -> BSDF* {
                                                            Float eta(0.0f);
                                                            Float mean = extra_params[0];
                                                            Float std_or =
                                                                extra_params[1];
                                                            Float rough = extra_params[2];
                                                            Float thick = extra_params[3];
                                                            Float sigma = extra_params[4];
                                                            Float albedo_diff_a =
                                                                extra_params[5];
                                                            Float albedo_flake_a =
                                                                extra_params[6];
                                                            bool use_specular =
                                                                extra_params[7] == 0
                                                                    ? true
                                                                    : false;
                                                            Float std_size =
                                                                extra_params[8];
                                                            Float ior = extra_params[9];
                                                            Float diff_particles_perc =
                                                                extra_params[10];
                                                            Float g1 = extra_params[11];
                                                            Float g2 = extra_params[12];
                                                            Float w_g = extra_params[13];
                                                            Float n_samples =
                                                                extra_params[14];
                                                            for (int i; i < 14; i++) {
                                                                std::cout
                                                                    << "val is "
                                                                    << extra_params[i]
                                                                    << std::endl;
                                                            }
                                                            Float albedo_flake_R =
                                                                extra_params[15];
                                                            Float albedo_flake_G =
                                                                extra_params[16];
                                                            Float albedo_flake_B =
                                                                extra_params[17];

                                                            Float albedo_diff_R =
                                                                extra_params[18];
                                                            Float albedo_diff_G =
                                                                extra_params[19];
                                                            Float albedo_diff_B =
                                                                extra_params[20];
                                                            std::vector<Float>
                                                                albedo_diff = {
                                                                    albedo_diff_R,
                                                                    albedo_diff_G,
                                                                    albedo_diff_B};
                                                            std::vector<Float>
                                                                albedo_flake = {
                                                                    albedo_flake_R,
                                                                    albedo_flake_G,
                                                                    albedo_flake_B};
                                                            return createLayeredMaterial(
                                                                si, alloc, mean, std_or,
                                                                std_size, rough, thick,
                                                                sigma, albedo_diff,
                                                                albedo_flake, ior,
                                                                diff_particles_perc, g1,
                                                                g2, w_g, use_specular,
                                                                int(n_samples));

                                                            // return
                                                            // createDielectric(si,
                                                            // alloc, ior,rough);
                                                        },

                                                        "EvalLayeredMaterial",
                                                        out_path.str(), extraParams);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

TEST(Measurements, Eval) {
    printf("Measuring the reflectance of a particular material\n");
    RNG rng;
    Float lambda_min = 380, lambda_max = 830;

    // Diffuse BSDF parameters
    SampledSpectrum reflectance(1.0f);
    // DiffuseBxDF bsdf = DiffuseBxDF(reflectance);

    // Conductor parameters

    SampledSpectrum eta(0.0f);
    SampledSpectrum k(1.0f);

    Float roughx = 0.4;

    // Float alphax = TrowbridgeReitzDistribution::RoughnessToAlpha(roughx);
    // Float alphay = TrowbridgeReitzDistribution::RoughnessToAlpha(roughy);
    TrowbridgeReitzDistribution distrib(roughx, roughx);
    // BeckmannDistribution distrib(roughx, roughx);

    bool reverseOrientation = false;
    auto t = std::make_shared<const Transform>(RotateX(360));
    auto tInv = std::make_shared<const Transform>(Inverse(*t));

    std::shared_ptr<Disk> disk =
        std::make_shared<Disk>(t.get(), tInv.get(), reverseOrientation, 0., 1., 0, 360.);
    Point3f origin(0.1, 1, 0);  // offset slightly so we don't hit center of disk
    Vector3f direction(0, -1, 0);
    Ray r(origin, direction);
    auto si = disk->Intersect(r);
    ASSERT_TRUE(si.has_value());
    // bsdf = createBSDF(si->intr, Allocator());
    auto m_si = si->intr;
    std::cout << "m_si " << m_si.shading.n << std::endl;
    auto m_alloc = Allocator();
    Float ior = 1.5;
    auto bsdf =
        m_alloc.new_object<BSDF>(Normal3f(0, 0, 1), m_si.shading.dpdu,
                                 m_alloc.new_object<ConductorBxDF>(distrib, eta, k));

    // ConductorBxDF bsdf = ConductorBxDF(distrib, eta, k);

    // Dielectric parameters
    /*Float eta = 1.5;

    Float roughx = 0.05, roughy = 0.05;

    Float alphax = TrowbridgeReitzDistribution::RoughnessToAlpha(roughx);
    Float alphay = TrowbridgeReitzDistribution::RoughnessToAlpha(roughy);

    TrowbridgeReitzDistribution distrib(alphax, alphay);

    DielectricBxDF bsdf = DielectricBxDF(eta, distrib);
*/
    // Low resolution test

    /*
    int theta_i_res = 32;
    int phi_i_res   = 32;
    int theta_o_res = 32;
    int phi_o_res   = 32;
    */

    // High resolution test
    int theta_i_res = 101;
    int phi_i_res = 1;
    int theta_o_res = 101;
    int phi_o_res = 1;

    // Angular sampling
    std::vector<Float> theta_i_values = linspace(0.0f, 100.0f, theta_i_res);
    std::vector<Float> theta_o_values = linspace(0.0f, 100.0f, theta_o_res);

    /*
    //Some issues for latlong and polar plots
    std::vector<Float> phi_i_values = linspace(0.0f, 360.0f, phi_i_res);
    std::vector<Float> phi_o_values = linspace(0.0f, 360.0f, phi_o_res);
    */

    std::vector<Float> phi_i_values = linspace(-180.0f, 180.0f, phi_i_res);
    std::vector<Float> phi_o_values = linspace(-180.0f, 180.0f, phi_o_res);

    FILE* fp = fopen("../../scripts/measurements/conductor_pbrt.txt", "w");

    ////#pragma omp parallel for

    /*
    // Full BSDF measurements
    for (auto phi_i: phi_i_values)
    {
        for (auto theta_i: theta_i_values)
        {
            //Vector3f wi = sphericalDirection(Radians(theta_i), Radians(phi_i));
            Vector3f wi = sphericalDirection(Radians(-theta_i), Radians(phi_i));

            for (auto phi_o: phi_o_values)
            {
                for (auto theta_o: theta_o_values)
                {

                    //SampledWavelengths lambda =
    SampledWavelengths::SampleVisible(rng.Uniform<Float>());//Random spectral samples
                    SampledWavelengths lambda =
    SampledWavelengths::SampleUniform(rng.Uniform<Float>(), lambda_min, lambda_max);
                    Vector3f wo = sphericalDirection(Radians(theta_o), Radians(phi_o));

                    //Evaluate BSDF
                    SampledSpectrum f = bsdf.f(wo, wi, TransportMode::Radiance) *
    AbsCosTheta(wi);

                    //std::cout << lambda[0] << std::endl;
                    //std::cout << f[0] << std::endl;

                    // Save reflectance (1 entry per wavelength)
                    for (int k = 0; k < NSpectrumSamples;k++){
                        fprintf(fp, "%f %f %f %f %f %f\n", phi_i, theta_i, phi_o, theta_o,
    lambda[k], f[k]);
                    }
                }
            }
        }
    }

    */

    // Fixed incident light example
    // Float phi_i = 0.0f, theta_i = -80.0f;
    Float phi_i = 0.0f;
    for (auto theta_i : theta_i_values) {
        Vector3f wi = sphericalDirection(Radians(theta_i), Radians(phi_i));
        for (auto phi_o : phi_o_values) {
            for (auto theta_o : theta_o_values) {
                // SampledWavelengths lambda =
                // SampledWavelengths::SampleVisible(rng.Uniform<Float>());//Random
                // spectral samples
                SampledWavelengths lambda = SampledWavelengths::SampleUniform(
                    rng.Uniform<Float>(), lambda_min, lambda_max);
                // wo is the Observation Vector.
                Vector3f wo = sphericalDirection(Radians(theta_o), Radians(phi_o));

                // Evaluate BSDF
                // NB: Should we use wo instead wi?
                // SampledSpectrum f = bsdf->f(wo, wi,TransportMode::Radiance) *
                // AbsCosTheta(wi);

                SampledSpectrum f =
                    bsdf->f(wi, wo, TransportMode::Radiance) * AbsCosTheta(wo);
                // SampledSpectrum f = bsdf.f(wo, wi, TransportMode::Radiance) *
                // AbsCosTheta(wo);
                //  SampledSpectrum f = bsdf.f(wi, wo, TransportMode::Radiance) *
                //  AbsCosTheta(wo);

                // std::cout << lambda[0] << std::endl;
                // std::cout << f[0] << std::endl;

                // Save reflectance (1 entry per wavelength)
                for (int k = 0; k < NSpectrumSamples; k++) {
                    fprintf(fp, "%f %f %f %f %f %f\n", phi_i, theta_i, phi_o, theta_o,
                            lambda[k], f[k]);
                }
            }
        }
    }

    fclose(fp);

    // Estimate reflected uniform incident radiance from bsdf
    /*Float ySum = 0;
    int count = 20000;

    if (NSpectrumSamples < 4) count *= 4;

    for (int i = 0; i < count; ++i) {

        //std::cout << lambda << std::endl;

        SampledWavelengths lambda =
    SampledWavelengths::SampleVisible(rng.Uniform<Float>());
        //Vector3f wi = SampleUniformSphere({rng.Uniform<Float>(), rng.Uniform<Float>()});
        Vector3f wo = SampleUniformSphere({rng.Uniform<Float>(), rng.Uniform<Float>()});

        SampledSpectrum reflectance(1.0f);

        DiffuseBxDF bsdf = DiffuseBxDF(reflectance);

        SampledSpectrum f = bsdf.f(wo, wi, TransportMode::Radiance) * AbsCosTheta(wi);
        ySum += f.y(lambda);
    }

    Float avg = ySum / (count * UniformSpherePDF());

    printf("Average reflectance: %f\n", avg);

    //EXPECT_TRUE(avg >= .95 && avg <= 1.05) << avg;
    */

    EXPECT_TRUE(true);
}

TEST(Measurements, MeasureBRDF) {
    // Read the material configuration file: model pararameters
    std::ifstream file("./scripts/measurements/model_parameters.json");
    Json::Value root;  // Root of JSON object
    Json::Reader reader;

    bool parsingSuccessful = reader.parse(file, root);
    if (!parsingSuccessful) {
        // Report the failure and their locations in the document.
        std::cout << "Failed to parse JSON: " << reader.getFormattedErrorMessages()
                  << std::endl;
        return;
    }

    // Extract model parameters
    // printf("Extracting parameters\n");

    // Float albedo = root["albedo"].asFloat();
    // Float albedo_flake = root["albedo_flake"].asFloat();
    // Float albedo_diff = root["albedo_diff"].asFloat();

    Float surf_rough = root["surf_roughs"].asFloat();
    Float thick = root["thickness"].asFloat();
    Float refl = root["refl_levels"].asFloat();
    Float mean = root["means_or"].asFloat();
    Float std_or = root["stds_or"].asFloat();
    Float sigma = root["mean_size"].asFloat();
    Float std_size = root["stds_size"].asFloat();
    Float ior = root["ior"].asFloat();
    Float diff_particles_perc = root["diff_particles_percs"].asFloat();
    Float g1 = root["g1"].asFloat();
    Float g2 = root["g2"].asFloat();
    Float w_g = root["w_g"].asFloat();

    Float albedo_diff_r = root["albedo_diff_r"].asFloat();
    Float albedo_diff_g = root["albedo_diff_g"].asFloat();
    Float albedo_diff_b = root["albedo_diff_b"].asFloat();

    Float albedo_flake_r = root["albedo_flake_r"].asFloat();
    Float albedo_flake_g = root["albedo_flake_g"].asFloat();
    Float albedo_flake_b = root["albedo_flake_b"].asFloat();
    
    // int n_samples = root["n_samples"].asInt();
    std::vector<Float> albedo_diff {albedo_diff_r,albedo_diff_g,albedo_diff_b};
    std::vector<Float> albedo_flake {albedo_flake_r,albedo_flake_g,albedo_flake_b};
    // = json_to_vector(root["albedo_diff"]);
    // // achromatic albedo
    // if (albedo_diff.size() == 1) {
    //     albedo_diff.push_back(albedo_diff[0]);
    //     albedo_diff.push_back(albedo_diff[0]);
    // }
    // std::vector<Float> albedo_flake = json_to_vector(root["albedo_flake"]);
    // // achromatic albedo
    // if (albedo_flake.size() == 1) {
    //     albedo_flake.push_back(albedo_flake[0]);
    //     albedo_flake.push_back(albedo_flake[0]);
    // }

    std::cout<<"alb flake"<< albedo_flake[0] <<std::endl;

    std::vector<Float> n_samples = json_to_vector(root["n_samples"]);
    if (n_samples.empty()) {
        n_samples.push_back(128.f);
    }
    // printf("Surface roughness = %f\n", surf_rough);
    printf("Sigma (Mean Size) = %f, Particle Percentage = %f\n", sigma,
           diff_particles_perc);

    // std::cout<<"Measure BRDF for fitting"<< std::endl;

    // Evaluate the BSDF

    std::ostringstream out_path;
    auto refl_type = refl == 0 ? "_spec" : "_diff";
    // Float extraParams[13];
    Float extraParams[21];

    // extraParams[5] = albedo; //albedo

    extraParams[0] = mean;
    extraParams[1] = std_or;
    extraParams[2] = surf_rough;  // Interface's roughness
    extraParams[3] = thick;
    extraParams[4] = sigma;
    extraParams[5] = albedo_diff[0];   // albedo
    extraParams[6] = albedo_flake[0];  // albedo
    extraParams[7] = refl;             // Specular (0) or Diffuse (1) Reflectance?
    extraParams[8] = std_size;
    extraParams[9] = ior;
    extraParams[10] = diff_particles_perc;
    extraParams[11] = g1;
    extraParams[12] = g2;
    extraParams[13] = w_g;
    extraParams[14] = n_samples[0];

    extraParams[15] = albedo_flake[0];
    extraParams[16] = albedo_flake[1];
    extraParams[17] = albedo_flake[2];

    extraParams[18] = albedo_diff[0];
    extraParams[19] = albedo_diff[1];
    extraParams[20] = albedo_diff[2];

    //printf("Albedo flakes = %f, Albedo diffuse = %f\n", albedo_flake, albedo_diff);

    /*
    out_path << "./scripts/measurements/fitting/layered_mean_" <<std::setprecision(4) <<
    mean <<"_std_or_" \
         << std_or << "_std_sz_"<< std_size << "_t_" <<thick <<
    refl_type<<"_s_"<<sigma<<"_r_"<<surf_rough<<"_a_"<<albedo<< "_ior_" \
         <<ior <<"_conc_"<< diff_particles_perc << "_g1_" <<g1 << "_g2_"<< g2 << "_wg_" <<
    w_g << ".txt";
    */

    out_path << "./scripts/measurements/fitting/fitting_model.txt";

    EvalBsdf(
        [](const SurfaceInteraction& si, Allocator alloc, Float* extra_params) -> BSDF* {
            Float eta(0.0f);
            Float mean = extra_params[0];
            Float std_or = extra_params[1];
            Float rough = extra_params[2];
            Float thick = extra_params[3];
            Float sigma = extra_params[4];
            Float albedo_diff_a = extra_params[5];
            Float albedo_flake_a = extra_params[6];
            bool use_specular = extra_params[7] == 0 ? true : false;
            Float std_size = extra_params[8];
            Float ior = extra_params[9];
            Float diff_particles_perc = extra_params[10];
            Float g1 = extra_params[11];
            Float g2 = extra_params[12];
            Float w_g = extra_params[13];
            Float n_samples = extra_params[14];

            Float albedo_flake_R = extra_params[15];
            Float albedo_flake_G = extra_params[16];
            Float albedo_flake_B = extra_params[17];

            Float albedo_diff_R = extra_params[18];
            Float albedo_diff_G = extra_params[19];
            Float albedo_diff_B = extra_params[20];
            std::vector<Float> albedo_diff = {albedo_diff_R, albedo_diff_G,
                                              albedo_diff_B};
            std::vector<Float> albedo_flake = {albedo_flake_R, albedo_flake_G,
                                               albedo_flake_B};
            
            // for (int i=0; i < 15; i++) {
            //     std::cout << "val is " << extra_params[i] << std::endl;
            // }

            return createLayeredMaterial(si, alloc, mean, std_or, std_size, rough, thick,
                                         sigma, albedo_diff, albedo_flake, ior,
                                         diff_particles_perc, g1, g2, w_g, use_specular,
                                         int(n_samples));

            // return createLayeredMaterial(si,alloc,mean,std_or,std_size,rough,thick,
            // sigma, albedo,ior,diff_particles_perc,g1,g2,w_g,use_specular);

            // return createDielectric(si, alloc, ior,rough);
        },
        "EvalLayeredMaterial", out_path.str(), extraParams);

    // TO DO: Precompute a list of parameters

    EXPECT_TRUE(true);
}
