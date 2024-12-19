# ESiWACE3 Trixi.jl service

## Instructions for terrabyte cluster

You need to get an account at https://docs.terrabyte.lrz.de/services/identity/get-account/
and set up two-factor authentication.

Documentation is available here: https://docs.terrabyte.lrz.de/


### Login
```shell
ssh login.terrabyte.lrz.de
```
You have storage space at `$HOME`, `$SCRATCH` (not backed up, temporary), and `$PROJECT`
(soon to come).


### Set up t8code

**TODO: once there is $PROJECT, this step can be skipped**

1. Load modules
   ```shell
   module purge
   module load spack/23.1.0
   module load gcc/12.2.0
   module load openmpi/4.1.5-gcc11
   ```
2. Change to scratch folder
   ```shell
   cd $SCRATCH
   ```
3. Clone the repository
   ```shell
   git clone --branch 'v3.0.1' --depth 1 https://github.com/DLR-AMR/t8code.git
   cd t8code
   git submodule init
   git submodule update
   ```
4. Build using cmake:
   ```shell
   module add cmake
   mkdir build
   cd build
   cmake \
     -DCMAKE_C_COMPILER=mpicc \
     -DCMAKE_CXX_COMPILER=mpicxx \
     -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_INSTALL_PREFIX="$SCRATCH/install/t8code" \
     -DT8CODE_BUILD_TESTS=OFF \
     -DT8CODE_BUILD_TUTORIALS=OFF \
     -DT8CODE_BUILD_EXAMPLES=OFF \
     -DT8CODE_BUILD_BENCHMARKS=OFF \
     -DT8CODE_ENABLE_MPI=ON
     ..
   nice make -j8
   nice make install -j8
   ```


### Set up Julia
Julia is not available on the cluster. We need to install it manually.
1. If there is no `.bashrc` or `.bash_profile` in your `$HOME` directory, create one
   ```
   touch $HOME/.bashrc
   ```
2. Use the official Julia installer:
   ```shell 
   curl -fsSL https://install.julialang.org | sh
   ```
   Accept the defaults. Once finished you will be told to source your `.bashrc` or re-login.
3. Julia should now be available
   ```shell
   julia --version
   ```
4. Install the 1.11 branch
   ```shell
   juliaup add 1.11
   ```


### Set up Trixi.jl
1. Clone the repository
   ```shell
   git clone https://github.com/benegee/Trixi.jl.git
   git switch lc/gpu-develop
   ```
2. Go to the `esiwace` directory. We collect necessary environmental settings in
   `profile`. Edit this file as neccessary and source it:
   ```shell
   . profile
   ```
3. The Julia project is configured by several files: `Project.toml` lists dependencies,
   `Manifest.toml` lists exact version numbers for all installed packages,
   `LocalPreferences.toml` contains advanced configuration options.
   It should only be necessary to adapt `LocalPreference.toml` to reflect the t8code
   installation path.
4. Open Julia via the `$JL` command and instantiate the project:
   ```shell
   $JL --project=. -e 'using Pkg; Pkg.instantiate()'
   ```
   This will take some time! Some packages might throw errors.


### Check installation
1. Make sure that everything is precompiled by running the following:
   ```shell
   $JL --project=. -e 'using OrdinaryDiffEq, Trixi'
   ```
   If there are still some errors, they might get resolved when running on compute nodes.
2. To test CUDA, first log in to a GPU node:
   ```shell
   salloc --cluster=hpda2 --partition=hpda2_testgpu --nodes=1 --ntasks-per-node=1 --gres=gpu:1 --time=00:30:00
   ```
   Then start Julia:
   ```shell
   $JL --project=. -e 'using CUDA; CUDA.versioninfo()'
   ```
   This should print
   ```
   CUDA runtime 11.8, local installation
   ...
   ```
   <!--
   If it fails, it might help to re-set the CUDA runtime:
   ```shell
   $JL --project=. -e 'using CUDA; CUDA.set_runtime_version!(VersionNumber(11,8); local_toolkit=true)
   ```
   
   pkg = Base.PkgId(Base.UUID("76a88914-d11a-5bdc-97e0-2f5a05c973a2"), "CUDA_Runtime_jll")
   Base.compilecache(pkg)
   -->


## Launch
1. SLURM jobscripts are in `jobscripts`. Edit as necessary. At least, you have to specify
   your mail address.
2. The actual simulation is configured in `run.jl` and based on Trixi.jl files in `elixirs`.
3. Send job to queue:
   ```shell
   sbatch jobscript/single_node.sh
   ```
