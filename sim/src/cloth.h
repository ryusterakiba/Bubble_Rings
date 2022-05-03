#ifndef CLOTH_H
#define CLOTH_H

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "clothMesh.h"
#include "collision/collisionObject.h"
#include "spring.h"

using namespace CGL;
using namespace std;

enum e_orientation { HORIZONTAL = 0, VERTICAL = 1 };

struct ClothParameters {
  ClothParameters() {}
  ClothParameters(bool enable_structural_constraints,
                  bool enable_shearing_constraints,
                  bool enable_bending_constraints, double damping,
                  double density, double ks)
      : enable_structural_constraints(enable_structural_constraints),
        enable_shearing_constraints(enable_shearing_constraints),
        enable_bending_constraints(enable_bending_constraints),
        damping(damping), density(density), ks(ks) {}
  ~ClothParameters() {}

  // Global simulation parameters

  bool enable_structural_constraints;
  bool enable_shearing_constraints;
  bool enable_bending_constraints;

  double damping;

  // Mass-spring parameters
  double density;
  double ks;
};

struct Cloth {
  Cloth() {}
  Cloth(int num_vertices, double initial_ring_radius);
  ~Cloth();

  void buildGrid();

  void simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                vector<Vector3D> external_accelerations,
                vector<CollisionObject *> *collision_objects);

  void reset();
  void deleteExtraPoints();
  void createExtraPoints();
  void buildClothMesh();

  void build_spatial_map();
  void self_collide(PointMass &pm, double simulation_steps);
  float hash_position(Vector3D pos);

  // Bubble ring physics functions
  void velocity();
  Vector3D biotsavart_edge(PointMass p, PointMass v0, PointMass v1);
  Vector3D induction(PointMass v0, PointMass v1, PointMass v2);
  double boussinesq(PointMass pm0, PointMass pm1);
  void modify_thickness();
  void resample();
  double volume();
  void burgers_flow(double delta_t);

  // Cloth properties
  double width;
  double height;
  int num_width_points;
  int num_height_points;
  double initial_ring_radius;
  int num_vertices;
  double circulation;
  double thickness;
  e_orientation orientation;
  double vol0;
  double min_dist;

  // Cloth components
  vector<PointMass> point_masses;
  vector<vector<int>> pinned;
  vector<Spring> springs;
  ClothMesh *clothMesh;

  // Spatial hashing
  unordered_map<float, vector<PointMass *> *> map;
};

#endif /* CLOTH_H */
