# Architecture

## Repository Structure

```
computational-structural-analysis/
├── CSA_2021_FRAME/          Frame analysis solver
│   ├── frame.f90            Source code (1026 lines)
│   ├── case-studies/        Example input/output files
│   └── Figures/             Deformed shape plots
├── CSA_2021_ISOP4/          ISOP4 element solver
│   ├── isop4.f90            Source code (1332 lines)
│   ├── case-studies/        Example input/output files
│   └── Figures/             Mesh, deformation, and stress plots
├── CSA_2021_BENCHMARKS/     Benchmark validation cases
│   ├── frame.f90            Frame solver variant for benchmarks
│   ├── case-studies/        10+ benchmark input/output files
│   └── Figures/             Benchmark result plots
└── docs/                    Documentation
```

## Frame Solver — Subroutine Call Graph

```
PROGRAM FRAME
│
├── GEOMET ─────────── Read nodes, sections, materials, elements
│
├── SCODE ──────────── Read restraints and links, number DOFs
│
├── LOADS ──────────── Read nodal loads, element loads, prestressing
│
├── ASSEMB ─────────── Assemble global stiffness matrix
│   └── MKK(NE) ───── Compute element stiffness + equivalent forces
│                      ├── Build 6×6 local K (Timoshenko beam)
│                      ├── Build 6×6 transformation matrix T
│                      ├── Transform: K_global = T·K·T^T
│                      └── Compute equivalent nodal forces:
│                          ├── Distributed loads (constant + linear)
│                          ├── Self-weight
│                          ├── Thermal loads (axial + flexural)
│                          └── Prestressing (parabolic profile)
│
├── JOINTS ─────────── Apply elastic restraints/links, imposed displacements
│
├── SOLVE ──────────── Gaussian elimination + back-substitution
│
├── STRESS ─────────── Recover element end forces: q = K·d - f_eq
│
└── PLOT ───────────── Interpolate displacements/forces along elements
    ├── F(I,AL,X) ──── Hermite shape functions
    ├── U(QSL,X,AL) ── Axial displacement interpolation
    └── V(QSL,X,AL) ── Transverse displacement interpolation
```

## ISOP4 Solver — Subroutine Call Graph

```
PROGRAM ISOP4
│
├── [if automatic mesh]
│   └── MESH ───────── Generate structured quad mesh from domain corners
│                      └── FF4(A,CSI,ETA) ── Isoparametric coordinate mapping
│
├── GEOMET ─────────── Read/write nodes, element types, connectivity
│
├── SCODE ──────────── Read restraints and links, number DOFs
│
├── LOADS ──────────── Read nodal loads, surface loads, thermal loads
│
├── ASSEMB ─────────── Assemble global stiffness matrix
│   └── MKK(NE) ───── Compute element stiffness + equivalent forces
│       ├── Build 3×3 constitutive matrix D (plane stress or strain)
│       ├── GAUSS(NGP,JGP,CC,WW) ── Gauss quadrature points/weights
│       ├── MKB(CSI,ETA,X,Y,DETJ,B) ── Strain-displacement B-matrix
│       │   ├── SHCSI(I,ETA) ── ∂N/∂ξ
│       │   ├── SHETA(I,CSI) ── ∂N/∂η
│       │   ├── FF4C(A,ETA) ── Σ aᵢ·∂Nᵢ/∂ξ (Jacobian terms)
│       │   └── FF4E(A,CSI) ── Σ aᵢ·∂Nᵢ/∂η (Jacobian terms)
│       ├── K = Σ t·B^T·D·B·det(J)·wᵢ·wⱼ  (numerical integration)
│       ├── Body force loads via SHAPEF(I,CSI,ETA)
│       └── Thermal equivalent forces via B^T·D·ε₀
│
├── SOLVE ──────────── Gaussian elimination + back-substitution
│
├── STRESS ─────────── Compute stresses/strains at:
│   │                  ├── Gauss points
│   │                  ├── Element centroid (ξ=0, η=0)
│   │                  └── Element vertices (ξ=±1, η=±1)
│   ├── Principal stress/strain (Mohr's circle)
│   ├── Deformation energy per element: U = ½·d^T·K·d
│   └── Nodal averaging (arithmetic / volumetric / energy-weighted)
│
└── PLOT ───────────── Write mesh, stresses, displacements, energy to PLOT.DAT
```

## Data Flow

```
                    ┌─────────────┐
                    │  input.txt  │
                    └──────┬──────┘
                           │
              ┌────────────▼────────────┐
              │    GEOMET + SCODE       │
              │  (mesh, BCs, DOF map)   │
              └────────────┬────────────┘
                           │
              ┌────────────▼────────────┐
              │        LOADS            │
              │  (forces → load vector) │
              └────────────┬────────────┘
                           │
              ┌────────────▼────────────┐
              │    ASSEMB + MKK         │
              │  K·d = F  (assembly)    │◄──── MKK.tmp (element matrices)
              └────────────┬────────────┘
                           │
              ┌────────────▼────────────┐
              │        SOLVE            │
              │  (Gaussian elimination) │
              └────────────┬────────────┘
                           │
              ┌────────────▼────────────┐
              │    STRESS + PLOT        │
              │  (post-processing)      │
              └────┬───────────────┬────┘
                   │               │
          ┌────────▼──────┐  ┌─────▼──────┐
          │  output.txt   │  │  PLOT.DAT  │
          │  (text results│  │  (visual-  │
          │   & tables)   │  │   ization) │
          └───────────────┘  └─────┬──────┘
                                   │
                          ┌────────▼────────┐
                          │  MATLAB viewer  │
                          └─────────────────┘
```

## Key Data Structures

### Frame (`GLOBALVAR` module)

| Variable | Type | Description |
|----------|------|-------------|
| `COORD(NNODE,2)` | Real | Node coordinates (X, Y) |
| `IN(NELE,2)` | Integer | Element connectivity (node 1, node 2) |
| `CSEC(NSEC,4)` | Real | Section properties (A, I, h, χ) |
| `CMAT(NMAT,4)` | Real | Material properties (E, ν, α, γ) |
| `IDOF(NNODE,3)` | Integer | DOF map (-1 = fixed, >0 = equation number) |
| `VK(NDOF,NDOF)` | Real | Global stiffness matrix |
| `VLOADS(NDOF)` | Real | Global load vector |
| `VDISP(NDOF)` | Real | Solution displacement vector |
| `ELOADS(NELE,5)` | Real | Element loads (Px, Py1, Py2, DTtop, DTbot) |
| `PRES(NELE,4)` | Real | Prestressing (P, e1, em, e2) |
| `ST(6,6)` | Real | Current element stiffness (local→global) |
| `EQFG(6)` | Real | Current element equivalent forces (global) |

### ISOP4 (`GLOBALVAR` module)

| Variable | Type | Description |
|----------|------|-------------|
| `COORD(NNODE,2)` | Real | Node coordinates (X, Y) |
| `IN(NELE,4)` | Integer | Element connectivity (4 nodes) |
| `CTYPE(NTYPE,5)` | Real | Element type properties (E, ν, α, γ, t) |
| `IDOF(NNODE,2)` | Integer | DOF map (-1 = fixed, >0 = equation number) |
| `VK(NDOF,NDOF)` | Real | Global stiffness matrix |
| `VLOADS(NDOF)` | Real | Global load vector |
| `VDISP(NDOF)` | Real | Solution displacement vector |
| `THERM(NELE)` | Real | Element temperature changes |
| `NGAUSS(NELE)` | Integer | Gauss points per element |
| `ST(8,8)` | Real | Current element stiffness |
| `EQF(8)` | Real | Current element equivalent forces |
| `SIGNOD(NNODE,4)` | Real | Averaged nodal stresses |
| `EPSNOD(NNODE,4)` | Real | Averaged nodal strains |
