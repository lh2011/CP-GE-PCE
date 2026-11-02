# AMIGA - A Magnetic Isogeometric Analysis Framework for Electric Motor Design
Welcome to AMIGA my friend! AMIGA is a Magnetic [Isogeometric Analysis](https://en.wikipedia.org/wiki/Isogeometric_analysis) Framework (AMIGA) for Electric Motor Design, written in **Matlab**. 

## Features
- **Static**, **transient** and **harmonic** simulations with IGA or lowest order FEM
- Fully parametric **NURBS** models for PMSMs / IMs / reluctance machines
- Electric network simulation with **MNA** or **MLA** (MONA in progress)
- Arbitrary coupling of fields and circuits with stranded/solid conductor models
- Non-conforming rotor/stator coupling for rotating geometries
- Importing of FEM meshes from [Gmsh](https://gmsh.info/) or [FEMM](https://www.femm.info/wiki/HomePage)
- Very **modular** structure and **intuitive** naming

## Installation
AMIGA requires MATLAB 2022b or later. Additionally, the **curve fitting toolbox** is required. All other dependencies are automatically downloaded when cloning recursively with
```
git clone --recurse-submodules git@git.rwth-aachen.de:creator/amiga.git
```
If you already cloned but forgot the recursive call, you can use
```
git submodule init 
git submodule update
```
to manually update the submodules. If you download manually, you also need to download the dependencies manually and put them in the submodules folder. These are
- [FEMM](https://github.com/cenit/FEMM)
- [GeoPDEs](https://github.com/rafavzqz/geopdes)
- [The NURBS package](https://octave.sourceforge.io/nurbs/)
- [Tensor Toolbox](https://github.com/andrewssobral/tensor_toolbox)

## Getting Started
Check out the examples folder, choose anything that catches your interest and press run!

## Authors and acknowledgment
TBD

## License
GNU General Public License v3.0 or later. See [COPYING](COPYING) to see the full text.

## Open Tasks / Next Steps
- [ ] [Plotting for lowest order FEM: speficy rectangular grid an not each triangle has own vector]
- [ ] [Evaluation Points for lowest order FEM]
- [ ] [Include Optimization again + tests]
- [ ] [Include Thermal Simulations again + tests]
- [ ] [Include Mechanical Simulations again + tests, add nonlinear material]
- [ ] [More IGA geometries - specifically one good synchronous reluctance machine]
- [ ] [Include MLA conductor models]
- [ ] [Include MONA elements and conductor models]

