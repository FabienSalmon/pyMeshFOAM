
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pyMeshFOAM

Software to ease the creation of meshes for OpenFOAM numerical simulations

pyMeshFOAM aims at performing meshes based on free open-source mesh generation software. More particularly, the available scripts are dedicated to the meshing process of 2D of 3D slender shapes within a box. The mesh of boundary layers is encompassed within the options. This numerical tool allows the simulation of fluid-structure interactions between fluid and deformable structures. Three open-source software codes are involved in such simulations: OpenFOAM, preCICE and CalculiX. This free mesh generation code is based on free tools such as blockMesh, snappyHexMesh, cfMesh and GMSH. The solid mesh interface matches automatically with the mesh of the fluid interface. Thus, there is no interpolation between the nodes of both meshes. It is however possible to use different interfaces since preCICE allows interpolations. It also enables fluid simulations, so only with OpenFOAM. The programming language of the scripts is Python.

## Getting started

There are two tutorials available as first examples: a 2D wing and a 3D wing (with the solid mesh).


