# ------------------------Project Installation-------------------------------------------------------------

# Clone repository and submodules (public repository)

git clone --recursive https://github.com/mmp/pbrt-v4.git

# Update submodules if you already clone the project (public repository)

git submodule init
git submodule update

# Compile project

mkdir build
cd build
cmake ..
make -j20

# Unit Tests
./pbrt_test

# Copy our plugins

cp cp -R pbrt-v4-plugin/* pbrt-v4

# Remember to add our unit tests to the CMakeLists.txt (add src/pbrt/sggx_test.cpp)

# ------------------------Running examples-------------------------------------------------------------

# Layered BSDF for sphere example
./pbrt ../../scenes/sheen_volume/sphere.pbrt

# SGGX volume for sphere example (different results as PBRT v3)
./pbrt ../../scenes/sheen_volume/cloth.pbrt

# SGGX volume for cloth example (not working)
./pbrt ../../scenes/sheen_volume/cloth.pbrt

# SGGX medium for bunny cloud example (not working)
./pbrt ../../scenes/bunny-cloud/bunny_cloud.pbrt

# SGGX normal distribution plots: Logging SGGX data first and then plot it with matplotlib
./pbrt_test --test-filter SGGX.Logging
python ./scripts/SGGX/sggx_plots.py

# -----------------------Reflectance measurement, fitting and visualization-----------------------------

# Measuring a BSDF
./pbrt_test --test-filter Measurements.Eval

# Plotting the measured BSDF
python ./scripts/measurements/reflectance_plotting2.py

# Fitting measurements

time python ./scripts/measurements/fitting_measurements_cma_rgb_patches.py -optimizer CMA-ES

# Download results from server
scp -r juanraul@10.3.27.209:/home/juanraul/cosmetic_project/scripts/measurements/fitting .