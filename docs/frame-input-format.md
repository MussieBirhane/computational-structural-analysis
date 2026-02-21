# Frame Analysis — Input File Format

The frame solver reads a single plain-text input file. Each section starts with a keyword header followed by a count and then data rows. All values are comma-separated.

The first line declares the unit system:

```
2D FRAME [m, kN]
```

---

## NODES

Defines the node coordinates.

```
NODES<NNode><IDNode,X,Y>
```

| Field  | Description              |
|--------|--------------------------|
| NNode  | Total number of nodes    |
| IDNode | Node ID                  |
| X      | X coordinate             |
| Y      | Y coordinate             |

**Example** — 2 nodes:

```
NODES<NNode><IDNode,X,Y>
2
1, 0, 0
2, 4, 0
```

---

## SECTIONS

Defines cross-section properties.

```
SECTIONS<NSec><IDSec, Area, Inertia, Depth, ShearCF>
```

| Field   | Description                                      |
|---------|--------------------------------------------------|
| NSec    | Number of section types                          |
| IDSec   | Section ID                                       |
| Area    | Cross-sectional area                             |
| Inertia | Second moment of area (I)                        |
| Depth   | Section depth (used for thermal gradient)        |
| ShearCF | Shear correction factor (0 = Euler-Bernoulli)    |

**Example** — concrete rectangular section:

```
SECTIONS<NSec><IDSec, Area, Inertia, Depth, ShearCF>
1
1, 0.12, 0.0016, 0.4, 0
```

---

## MATERIALS

Defines material properties.

```
MATERIALS<NMat><IDMat, E, v, Thermal, SelfWeight>
```

| Field      | Description                             |
|------------|-----------------------------------------|
| NMat       | Number of material types                |
| IDMat      | Material ID                             |
| E          | Young's modulus                          |
| v          | Poisson's ratio                         |
| Thermal    | Thermal expansion coefficient (alpha)   |
| SelfWeight | Unit weight (for self-weight loading)   |

**Example** — concrete and steel:

```
MATERIALS<NMat><IDMat, E, v, Thermal, SelfWeight>
2
1, 30000000, 0.2, 0, 25
2, 210000000, 0.3, 0.000012, 78.5
```

---

## ELEMENTS

Defines element connectivity and assigns section/material properties.

```
ELEMENTS<NEle><IDElem, Node1, Node2, Section, Material>
```

| Field    | Description                    |
|----------|--------------------------------|
| NEle     | Number of elements             |
| IDElem   | Element ID                     |
| Node1    | Start node                     |
| Node2    | End node                       |
| Section  | Section ID (from SECTIONS)     |
| Material | Material ID (from MATERIALS)   |

**Example** — 3 elements:

```
ELEMENTS<NEle><IDElem, Node1, Node2, Section, Material>
3
1, 1, 2, 1, 1
2, 2, 3, 2, 2
3, 4, 5, 1, 1
```

---

## RESTRAINTS

Defines fixed (zero-displacement) boundary conditions.

```
RESTRAINTS<NRes><Node, Dir>
```

| Field | Description                            |
|-------|----------------------------------------|
| NRes  | Number of restrained DOFs              |
| Node  | Node ID                                |
| Dir   | Direction: 1 = Ux, 2 = Uy, 3 = Rotz   |

**Example** — fully fixed support at node 1:

```
RESTRAINTS<NRes><Node, Dir>
3
1, 1
1, 2
1, 3
```

---

## LINKS

Defines master-slave DOF links. The slave node's DOF is tied to the master node's DOF.

```
LINKS<NLinks><Master, Slave, Dir>
```

| Field  | Description                            |
|--------|----------------------------------------|
| NLinks | Number of links                        |
| Master | Master node ID                         |
| Slave  | Slave node ID                          |
| Dir    | Direction: 1 = Ux, 2 = Uy, 3 = Rotz   |

**Example** — nodes 3 and 4 share translations (internal hinge):

```
LINKS<NLinks><Master, Slave, Dir>
2
3, 4, 1
3, 4, 2
```

This links Ux and Uy of node 4 to node 3, but leaves rotation independent — creating an internal hinge.

---

## NODAL_LOADS

Applies concentrated forces or moments at nodes.

```
NODAL_LOADS<NLoads><Node, Dir, Force>
```

| Field | Description                                        |
|-------|----------------------------------------------------|
| NLoads| Number of nodal loads                              |
| Node  | Node ID                                            |
| Dir   | Direction: 1 = Fx, 2 = Fy, 3 = Mz                 |
| Force | Magnitude (positive follows global axis direction) |

**Example** — 50 kN horizontal force at node 2:

```
NODAL_LOADS<NLoads><Node, Dir, Force>
1
2, 1, 50
```

---

## ELEMENT_LOADS

Applies distributed loads and thermal effects along elements.

```
ELEMENT_LOADS<NLoads><Elem, Px, Py1, Py2, DTtop, DTbot>
```

| Field | Description                                              |
|-------|----------------------------------------------------------|
| NLoads| Number of loaded elements                                |
| Elem  | Element ID                                               |
| Px    | Axial distributed load (constant along element)          |
| Py1   | Transverse distributed load at start node                |
| Py2   | Transverse distributed load at end node                  |
| DTtop | Temperature change at top fiber                          |
| DTbot | Temperature change at bottom fiber                       |

Transverse loads vary linearly from Py1 to Py2. Thermal loads produce axial strain (uniform component) and curvature (gradient component).

**Example** — linearly varying transverse load with uniform temperature rise:

```
ELEMENT_LOADS<NLoads><Elem, Px, Py1, Py2, DTtop, DTbot>
1
2, 0, -25, -45, 20, 20
```

---

## PRESTRESSING ELEMENTS

Defines prestressing tendons with parabolic eccentricity profiles.

```
PRESTRESSING ELEMENTS<NEp><IDElep, NTend><IDTend, P, e1, em, e2>
```

| Field  | Description                                     |
|--------|-------------------------------------------------|
| NEp    | Number of elements with prestressing            |
| IDElep | Element ID                                      |
| NTend  | Number of tendons in this element               |
| IDTend | Tendon ID                                       |
| P      | Prestressing force                              |
| e1     | Eccentricity at start node (+ above centroid)   |
| em     | Eccentricity at midspan                         |
| e2     | Eccentricity at end node                        |

**Example** — element 1 with two parabolic tendons:

```
PRESTRESSING ELEMENTS<NEp><IDElep, NTend><IDTend, P, e1, em, e2>
1
1, 2
1, 10, 0.15, -0.2, 0.15
2, 10, -0.05, -0.4, -0.05
```

---

## ELASTIC RESTRAINTS

Applies spring supports at nodes.

```
ELASTIC RESTRAINTS<Number><Node, Dir, Stiffness>
```

| Field     | Description                            |
|-----------|----------------------------------------|
| Number    | Number of elastic restraints           |
| Node      | Node ID                                |
| Dir       | Direction: 1 = Ux, 2 = Uy, 3 = Rotz   |
| Stiffness | Spring stiffness                       |

---

## ELASTIC LINKS

Applies spring connections between two nodes.

```
ELASTIC LINKS<Number><N1, N2, Dir, Stiffness>
```

| Field     | Description                            |
|-----------|----------------------------------------|
| Number    | Number of elastic links                |
| N1        | First node ID                          |
| N2        | Second node ID                         |
| Dir       | Direction: 1 = Ux, 2 = Uy, 3 = Rotz   |
| Stiffness | Spring stiffness                       |

---

## IMPOSED DISPLACEMENT

Prescribes an absolute displacement at a node.

```
IMPOSED DISPLACEMENT<Number><Node, Dir, Displacement>
```

| Field        | Description                            |
|--------------|----------------------------------------|
| Number       | Number of imposed displacements        |
| Node         | Node ID                                |
| Dir          | Direction: 1 = Ux, 2 = Uy, 3 = Rotz   |
| Displacement | Prescribed displacement value          |

---

## IMPOSED RELATIVE DISPLACEMENT

Prescribes a relative displacement between two nodes.

```
IMPOSED RELATIVE DISPLACEMENT<Number><N1, N2, Dir, Displacement>
```

| Field        | Description                            |
|--------------|----------------------------------------|
| Number       | Number of imposed relative displacements|
| N1           | First node ID                          |
| N2           | Second node ID                         |
| Dir          | Direction: 1 = Ux, 2 = Uy, 3 = Rotz   |
| Displacement | Prescribed relative displacement       |

---

## Complete Example

A portal frame with a horizontal beam, one hinge, distributed load, prestressing, and a horizontal point load:

```
2D FRAME [m, kN]

NODES<NNode><IDNode,X,Y>
5
1, 0, 0
2, 0, 4
3, 6, 4
4, 6, 4
5, 6, 0

SECTIONS<NSec><IDSec, Area, Inertia, Depth, ShearCF>
2
1, 0.12, 0.0016, 0.4, 0
2, 0.00538, 0.00008356, 0.3, 0

MATERIALS<NMat><IDMat, E, v, Thermal, SelfWeight>
2
1, 30000000, 0.2, 0, 25
2, 210000000, 0.3, 0.000012, 78.5

ELEMENTS<NEle><IDElem, Node1, Node2, Section, Material>
3
1, 1, 2, 1, 1
2, 2, 3, 2, 2
3, 4, 5, 1, 1

RESTRAINTS<NRes><Node, Dir>
5
1, 1
1, 2
5, 1
5, 2
5, 3

LINKS<NLinks><Master, Slave, Dir>
2
3, 4, 1
3, 4, 2

NODAL_LOADS<NLoads><Node, Dir, Force>
1
2, 1, 50

ELEMENT_LOADS<NLoads><Elem, Px, Py1, Py2, DTtop, DTbot>
1
2, 0, -25, -45, 20, 20

PRESTRESSING ELEMENTS<NEp><IDElep, NTend><IDTend, P, e1, em, e2>
1
1, 2
1, 10, 0.15, -0.2, 0.15
2, 10, -0.05, -0.4, -0.05

ELASTIC RESTRAINTS<Number><Node, Dir, Stiffness>
0
ELASTIC LINKS<Number><N1, N2, Dir, Stiffness>
0
IMPOSED DISPLACEMENT<Number><Node, Dir, Displacement>
0
IMPOSED RELATIVE DISPLACEMENT<Number><N1, N2, Dir, Displacement>
0
```
