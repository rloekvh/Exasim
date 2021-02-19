# specify an Exasim version to run
version = "Version0.3";

# External packages
using Revise, DelimitedFiles, SymPy

# Add Exasim to Julia search path
cdir = pwd(); ii = findlast("Exasim", cdir);
include(cdir[1:ii[end]] * "/Installation/setpath.jl");

# Exasim packages
using Preprocessing, Mesh, Gencode, Postprocessing

# create pde structure and mesh structure
pde, mesh = Preprocessing.initializeexasim(version);

# Define PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";            # ModelC, ModelD, ModelW
include("pdemodel.jl");          # include the PDE model file

# Set discretization parameters, physical parameters, and solver parameters
# NOTE: useful parameters
pde.porder = 3;          # polynomial degree
pde.torder = 3;          # time-stepping order of accuracy
pde.nstage = 3;          # time-stepping number of stages
pde.dt = 0.05*ones(100);   # time step sizes
#pde.dt = 0.05*ones(20);
pde.soltime = collect(1:length(pde.dt)); # steps at which solution are collected
pde.visdt = 0.05; # visualization timestep size

gam = 1.4;                      # specific heat ratio
Re = 1600;                      # Reynolds number
Pr = 0.71;                      # Prandtl number    
Minf = 0.2;                     # Mach number
pde.physicsparam = [gam Re Pr Minf];
pde.tau = [2.0];                # DG stabilization parameter
pde.GMRESrestart=40;
pde.linearsolvertol=0.0001;
pde.linearsolveriter=41;
pde.precMatrixType=2;
pde.NLiter=2;

# Choose computing platform and set number of processors
#pde.platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;                # number of MPI processors
pde.gpucompiler = "/home/software/cuda/9.0/bin/nvcc"
#pde.gpucompiler = "/home/software/spack/llvm/7.0.1-ggobltsozzkw3rqjnslhigxqkjxrafnh/bin/clang++";
#pde.gpuflags = " --cuda-path=/path/to/cuda/version "
pde.gpuflags = "-Wl,-rpath, -I/home/software/cuda/9.0/include -L/home/software/cuda/9.0/lib64 -lcudart -lcublas";
pde.cpuflags = pde.cpuflags * " -g -pg"

# create a linear mesh for a cube domain
#mesh.p, mesh.t = Mesh.cubemesh(16,16,16,1);
mesh.p, mesh.t = Mesh.cubemesh(8,8,8, 1);
mesh.p = 2*pi*mesh.p;
# expressions for disjoint boundaries
mesh.boundaryexpr = [p -> (p[2,:] .< 1e-3), p -> (p[1,:] .> 2*pi-1e-3), p -> (p[2,:] .> 2*pi-1e-3), p -> (p[1,:] .< 1e-3), p -> (p[3,:] .< 1e-3), p -> (p[3,:] .> 2*pi-1e-3)];
mesh.boundarycondition = [1 1 1 1 1 1]; # Set boundary condition for each disjoint boundary
mesh.periodicexpr = [2 p->p[2:3,:] 4 p->p[2:3,:]; 1 p->p[[1,3],:] 3 p->p[[1,3],:]; 5 p->p[[1,2],:] 6 p->p[[1,2],:]];

# if you just want to generate code to send to a cluster
#compilerstr, pde, mesh, master, dmd = Postprocessing.producecode(pde, mesh);

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh,~,~,~,~  = Postprocessing.exasim(pde,mesh);

print("Done!");

