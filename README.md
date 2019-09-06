# Triangulation
C++ library implementing triangulation methods described by Hartley &amp; Sturm.

## Getting started
The shared library can be built and installed using CMake.
```
mkdir build && cd build
cmake ..
make && make install
```
## Usage
There are 6 triangulation methods in total:
- Linear-LS
- Linear-Eigen
- Iterative-LS
- Iterative-Eigen
- Poly
- PolyAbs

For each method there is a separate class in `Triangulation` namespace with constructor that expects two camera projection matrices **P0**, **P1**. The points can be triangulated using the `triangulate` method which expect either two points or two vectors of points.
```
#include <Triangulation/LinearLS.h>
#include <Triangulation/IterativeLS.h>
#include <Triangulation/LinearEigen.h>
#include <Triangulation/IterativeEigen.h>
#include <Triangulation/Poly.h>
#include <Triangulation/PolyAbs.h>
...
cv::Mat P0, P1;             // 3x4 camera projection matrices
cv::Point2d point0, point1; // triangulated image points
cv::Point3d point;          // triangulation result
...
point = Triangulation::LinearLS(P0, P1).triangulate(point0, point1);
point = Triangulation::IterativeLS(P0, P1).triangulate(point0, point1);
point = Triangulation::LinearEigen(P0, P1).triangulate(point0, point1);
point = Triangulation::IterativeEigen(P0, P1).triangulate(point0, point1);
point = Triangulation::Poly(P0, P1).triangulate(point0, point1);
point = Triangulation::PolyAbs(P0, P1).triangulate(point0, point1);
...
```
## Citation
If you find this library useful, please cite:

HARTLEY, Richard I.; STURM, Peter. Triangulation. Computer vision and image understanding, 1997, 68.2: 146-157.
