# ISOP4 Element Analysis — Input File Format

The ISOP4 solver reads a single plain-text input file. Each section starts with a keyword header followed by data. All values are comma-separated.

The first line declares the unit system:

```
2D ISOP4 FINITE ELEMENT [m, kN]
```

---

## TYPE OF PROBLEM

Selects the 2D formulation.

```
TYPE OF PROBLEM(0 = PLANE STRESS, 1 = PLANE STRAIN)
0
```

| Value | Formulation  | Assumption                          |
|-------|-------------|--------------------------------------|
| 0     | Plane stress | Out-of-plane stress is zero (thin)  |
| 1     | Plane strain | Out-of-plane strain is zero (thick) |

---

## TYPE OF MESH

Selects manual or automatic mesh generation.

```
TYPE OF MESH <0 = STANDARD, 1 = AUTOMATIC>
1
```

| Value | Mode      | Description                                                  |
|-------|-----------|--------------------------------------------------------------|
| 0     | Standard  | Nodes and elements are defined manually                      |
| 1     | Automatic | Mesh is generated from domain corner coordinates and subdivisions |

---

## DOMAINS (Automatic mesh only)

Defines rectangular domains by their four corner coordinates and subdivision counts. The mesh generator creates a structured grid of quad elements within each domain using isoparametric mapping.

```
DOMAINS<NDom><Dom, X1, Y1, X2, Y2, X3, Y3, X4, Y4, NSubx, NSuby>
```

| Field       | Description                         |
|-------------|-------------------------------------|
| NDom        | Number of domains                   |
| Dom         | Domain ID                           |
| X1, Y1      | Corner 1 (bottom-left)             |
| X2, Y2      | Corner 2 (bottom-right)            |
| X3, Y3      | Corner 3 (top-right)               |
| X4, Y4      | Corner 4 (top-left)                |
| NSubx       | Number of element subdivisions in X |
| NSuby       | Number of element subdivisions in Y |

**Example** — trapezoidal domain, 10x10 elements:

```
DOMAINS<NDom><Dom, X1, Y1, X2, Y2, X3, Y3, X4, Y4, NSubx, NSuby>
1
1, 0, 0, 10, 0, 12, 6, 2, 4, 10, 10
```

The corners do not need to form a rectangle — the isoparametric mapping handles arbitrary quadrilateral domains.

---

## CHARACTERISTICS

Defines element type properties.

```
CHARACTERSTICS<NType><Type, E, v, alpha, self-weight, t>
```

| Field       | Description                          |
|-------------|--------------------------------------|
| NType       | Number of element types              |
| Type        | Type ID                              |
| E           | Young's modulus                      |
| v           | Poisson's ratio                      |
| alpha       | Thermal expansion coefficient        |
| self-weight | Unit weight (for body force loading) |
| t           | Element thickness                    |

**Example** — concrete plate, 1 m thick:

```
CHARACTERSTICS<NType><Type, E, v, alpha, self-weight, t>
1
1, 30000000, 0.25, 0, 0, 1
```

---

## PROPERTY OF DOMAIN

Assigns an element type and Gauss quadrature order to each domain.

```
PROPERTY OF DOMAIN<IDDom, Type, NGauss>
```

| Field  | Description                                         |
|--------|-----------------------------------------------------|
| IDDom  | Domain ID                                           |
| Type   | Element type ID (from CHARACTERISTICS)              |
| NGauss | Number of Gauss points per direction (1, 2, 3, or 4)|

**Example** — domain 1 uses type 1 with 2x2 integration:

```
PROPERTY OF DOMAIN<IDDom, Type, NGauss>
1, 1, 2
```

---

## RESTRAINTS

Defines fixed (zero-displacement) boundary conditions.

```
RESTRAINTS<NRes><Node, Dir>
```

| Field | Description                  |
|-------|------------------------------|
| NRes  | Number of restrained DOFs    |
| Node  | Node ID                      |
| Dir   | Direction: 1 = Ux, 2 = Uy    |

**Example** — pin support at node 1:

```
RESTRAINTS<NRes><Node, Dir>
2
1, 1
1, 2
```

---

## LINKS

Defines master-slave DOF links.

```
LINKS<NLinks><Master, Slave, Dir>
```

| Field  | Description                  |
|--------|------------------------------|
| NLinks | Number of links              |
| Master | Master node ID               |
| Slave  | Slave node ID                |
| Dir    | Direction: 1 = Ux, 2 = Uy    |

---

## NODAL LOADS

Applies concentrated forces at nodes.

```
NODAL LOADS<NLoads><Node, Dir, Force>
```

| Field  | Description                                        |
|--------|----------------------------------------------------|
| NLoads | Number of nodal loads                              |
| Node   | Node ID                                            |
| Dir    | Direction: 1 = Fx, 2 = Fy                          |
| Force  | Magnitude (positive follows global axis direction) |

**Example** — 100 kN in X and -100 kN in Y at node 3:

```
NODAL LOADS<NLoads><Node, Dir, Force>
2
3, 1, 100
3, 2, -100
```

---

## SURFACE LOADS

Applies distributed loads along an element edge between two nodes.

```
SURFACE LOADS<NLoads><IDSide, N1, N2, Px1, Px2, Py1, Py2>
```

| Field  | Description                                     |
|--------|-------------------------------------------------|
| NLoads | Number of surface loads                         |
| IDSide | Load ID                                         |
| N1     | Start node of the loaded edge                   |
| N2     | End node of the loaded edge                     |
| Px1    | Horizontal traction at N1                       |
| Px2    | Horizontal traction at N2                       |
| Py1    | Vertical traction at N1                         |
| Py2    | Vertical traction at N2                         |

Tractions vary linearly from N1 to N2.

**Example** — uniform vertical traction of 10 kN/m along the edge from node 3 to 4:

```
SURFACE LOADS<NLoads><IDSide, N1, N2, Px1, Px2, Py1, Py2>
1
1, 3, 4, 0, 0, 10, 10
```

---

## THERMAL LOADS

Applies uniform temperature change to elements.

```
THERMAL LOADS<Nloads><Elem, DT>
```

| Field  | Description                     |
|--------|---------------------------------|
| Nloads | Number of thermally loaded elements |
| Elem   | Element ID                      |
| DT     | Temperature change              |

---

## NODAL AVERAGING

Selects the method for averaging element stresses/strains at shared nodes.

```
NODAL AVARAGING <0 = Arithmetic, 1 = Volumetric, 2 = Energy>
0
```

| Value | Method      | Description                                       |
|-------|-------------|---------------------------------------------------|
| 0     | Arithmetic  | Simple average of contributing element values      |
| 1     | Volumetric  | Weighted by element volume                         |
| 2     | Energy      | Weighted by element deformation energy             |

---

## Complete Example

A trapezoidal plate under nodal forces and surface traction, plane stress, automatic mesh:

```
2D ISOP4 FINITE ELEMENT [m, kN]

TYPE OF PROBLEM(0 = PLANE STRESS, 1 = PLANE STRAIN)
0

TYPE OF MESH <0 = STANDARD, 1 = AUTOMATIC>
1

DOMAINS<NDom><Dom, X1, Y1, X2, Y2, X3, Y3, X4, Y4, NSubx, NSuby>
1
1, 0, 0, 10, 0, 12, 6, 2, 4, 10, 10

CHARACTERSTICS<NType><Type, E, v, alpha, self-weight, t>
1
1, 30000000, 0.25, 0, 0, 1

PROPERTY OF DOMAIN<IDDom, Type, NGauss>
1, 1, 2

RESTRAINTS<NRes><Node, Dir>
4
1, 1
1, 2
25, 1
25, 2

LINKS<NLinks><Master, Slave, Dir>
0

NODAL LOADS<NLoads><Node, Dir, Force>
2
3, 1, 100
3, 2, -100

SURFACE LOADS<NLoads><IDSide, N1, N2, Px1, Px2, Py1, Py2>
1
1, 3, 4, 0, 0, 10, 10

THERMAL LOADS<Nloads><Elem, DT>
0

NODAL AVARAGING <0 = Arithmetic, 1 = Volumetric, 2 = Energy>
0
```
