#ifndef RINGVERTEX_H
#define RINGVERTEX_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"

using namespace CGL;

// Forward declarations
class Halfedge;

struct RingVertex {
    RingVertex(Vector3D position, double thickness)
      : position(position), start_position(position), last_position(position),
        temp_position(position), thickness(thickness) {}

  // static values
  Vector3D start_position;
  double circulation = 2; // for now, static

  // dynamic values
  Vector3D position;
  Vector3D temp_position;
  Vector3D last_position;
  Vector3D velocity;

  // RK4 integration
  Vector3D K1;
  Vector3D K2;
  Vector3D K3;
  Vector3D K4;

  double thickness;
  // mesh reference
  Halfedge *halfedge;
};

#endif /* RINGVERTEX_H */
