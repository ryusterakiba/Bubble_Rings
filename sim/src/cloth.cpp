#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(int num_vertices, double initial_ring_radius) {
  this->num_vertices = num_vertices;
  this->initial_ring_radius = initial_ring_radius;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
  // TODO (Part 1): Build a grid of masses and springs.
    // create num_width_points * num_height_points point masses, first starts at (0, 0) and last at (width, height)
    
    // point masses
    //for (double a = 0.; a <= this->width; a+= this->width / (double(num_width_points) - 1)) {
      //  for (double b = 0.; b <= this->height; b+= this->height / (double(num_height_points) - 1)) {
    double angle = 2 * PI / num_vertices;
    for (int i = 0; i < num_vertices; i++) {
        Vector3D position;
        position.x = initial_ring_radius * cos(angle * i);
        position.z = initial_ring_radius * sin(angle * i);
        position.y = 0;
        PointMass pm = PointMass(position, false, 0.2);
        point_masses.emplace_back(pm);
        
        // make a circle around the pointmass
        double new_angle = 2 * PI / num_vertices;
        for (int j = 0; j < num_vertices; j++) {
            Vector3D position;
            position.y = pm.thickness / 2 * sin(new_angle * j);
            position.x = pm.thickness / 2 * sin(new_angle * j) * cos(angle * i);
            position.y = pm.thickness / 2 * sin(new_angle * j) * sin(angle * i);
            PointMass new_pm = PointMass(position, false, 0.0);
            point_masses.emplace_back(new_pm);
        }
        int point_index = 0;
        Spring left_spring = Spring(&(this->point_masses[point_masses.size() - 1]), &(this->point_masses[point_masses.size() - num_vertices]), STRUCTURAL);
        springs.emplace_back(left_spring);
        while (point_index < num_vertices - 1) {
            left_spring = Spring(&(this->point_masses[point_masses.size() - num_vertices + point_index]), &(this->point_masses[point_masses.size() - num_vertices + point_index + 1]), STRUCTURAL);
            springs.emplace_back(left_spring);
            point_index += 1;
        }
    }
    
    int point_index = 0;
    
    Spring left_spring = Spring(&(this->point_masses[point_masses.size() - num_vertices - 1]), &(this->point_masses[0]), STRUCTURAL);
    springs.emplace_back(left_spring);

    while (point_index < num_vertices - 1) {
        left_spring = Spring(&(this->point_masses[point_index * (num_vertices + 1)]), &(this->point_masses[(point_index + 1) * (num_vertices + 1)]), STRUCTURAL);
        springs.emplace_back(left_spring);
        point_index += 1;
    }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  // This function performs 1 timestep of Runge Kutta
  double delta_t = 1.0f / frames_per_sec / simulation_steps;
  // TODO: fill this in
}

/**
 * Update the point velocities.
 */
void Cloth::velocity() {

}

/**
 * Biosavart term for velocity calculations.
 */
double Cloth::biotsavart_edge() {

}

/**
 * Induction term for velocity calculations.
 */
double Cloth::induction() {

}

/**
 * Boussinesq term for velocity calculations.
 */
double Cloth::boussinesq(PointMass pm0, PointMass pm1) {

}

/**
 * Updates the thickness at each edge (stored on PointMasses).
 * Uses each PointMass's position and last_position.
 * Sets values for each PointMass's thickness.
 */
void Cloth::modify_thickness() {

}

/**
 * Resamples the number of the points on the bubble ring.
 * Modifies positions, thicknesses, C values.
 */
void Cloth::resample(double min_dist) {

}

/**
 * Calculates the volume of the bubble ring.
 * Uses the positions and thicknesses stored in each PointMass.
 */
double Cloth::volume() {

}

/**
 * Modifies thicknesses according Burgers flow.
 */
void Cloth::burgers_flow(double delta_t) {

}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.
    for (PointMass& pm : this->point_masses) {
        float hash = hash_position(pm.position);
        if (map.count(hash)) {
            // key exists
            map[hash]->emplace_back(&pm);
        } else {
            vector<PointMass*> *bucket = new vector<PointMass*>;
            bucket->emplace_back(&pm);
            map[hash] = bucket;
        }
    }

}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.
    float hash = hash_position(pm.position);
    Vector3D average_correction = Vector3D(0., 0., 0.);
    int count = 0;
    vector<PointMass*> *bucket = map[hash];
    
    for (auto candidate_pm : *bucket) {
        if (candidate_pm != &pm) {
            double distance = (pm.position - candidate_pm->position).norm();
            if (distance < (2 * thickness)) {
                count += 1;
                Vector3D correction_vector = 2 * thickness * (pm.position - candidate_pm->position).unit() - (pm.position - candidate_pm->position);
//                printf("%f %f %f \n", pm.position.x, pm.position.y, pm.position.z);
                average_correction += correction_vector;
            }
        }
    }
    if (count != 0) {
        average_correction = average_correction / count;
        pm.position = pm.position + average_correction / simulation_steps;
    }

}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
    float w = 3 * width / num_width_points;
    float h = 3 * height / num_height_points;
    float t = max(w, h);
    
    int w_coord = int(pos.x / w);
    int h_coord = int(pos.y / h);
    int t_coord = int(pos.z / t);
    
    return float(61379 * w_coord + 52609 * h_coord + 84347 * t_coord);
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;
    
  //

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B, 
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D, 
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
