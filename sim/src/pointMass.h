#ifndef POINTMASS_H
#define POINTMASS_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"

using namespace CGL;

// Forward declarations
class Halfedge;

struct PointMass {
  PointMass(Vector3D position, bool pinned)
      : pinned(pinned), start_position(position), position(position),
        last_position(position) {}

  PointMass(Vector3D position, bool pinned, double thickness)
      : pinned(pinned), start_position(position), position(position),
        last_position(position), thickness(thickness) {}

  Vector3D normal();
  Vector3D velocity(double delta_t) {
    return (position - last_position) / delta_t;
  }

  // static values
  bool pinned;
  Vector3D start_position;
  double circulation = 2; // for now, static

  // dynamic values
  Vector3D position;
  Vector3D temp_position;
  Vector3D last_position;
  Vector3D forces;
  Vector3D point_velocity;

  // Other physics
  double grav;
  double length;

  // RK4 integration
  Vector3D K1;
  Vector3D K2;
  Vector3D K3;
  Vector3D K4;

  double thickness;
  // mesh reference
  Halfedge *halfedge;
};

#endif /* POINTMASS_H */
