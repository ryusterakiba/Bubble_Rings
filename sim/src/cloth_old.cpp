#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

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
    int count = 0;
    for (int y = 0; y < this->num_height_points; y += 1) {
        for (int x = 0; x < this->num_width_points; x += 1) {
            Vector3D position;
            double a = double(x) / double(this->num_width_points - 1) * double(this->width);
            double b = double(y) / double(this->num_height_points - 1) * double(this->height);
            if (orientation == HORIZONTAL) {
                position.y = 1.;
                position.x = a;
                position.z = b;
            } else {
                float random = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2.0)));
                position.z = random / 1000.;
                position.x = a;
                position.y = b;
            }
            bool found = false;
            for (auto vec : pinned) {
                if (vec[0] == x && vec[1] == y) {
                    PointMass pm = PointMass(position, true);
                    point_masses.emplace_back(pm);
                    found = true;
                    break;
                }
            }
            if (found == false) {
                PointMass pm = PointMass(position, false);
                point_masses.emplace_back(pm);
            }
            count += 1;
           
        }
    }
    
    // constraints
    int point_index = 0;
    for (PointMass& pm : this->point_masses) {
        // make structural constraints
        
        if (point_index % num_width_points != 0) {
            Spring left_spring = Spring(&(this->point_masses[point_index]), &(this->point_masses[point_index - 1]), STRUCTURAL);
            springs.emplace_back(left_spring);
        }
        
        if (point_index >= num_width_points) {
            Spring up_spring = Spring(&(this->point_masses[point_index]), &(this->point_masses[point_index - num_width_points]), STRUCTURAL);
            springs.emplace_back(up_spring);
        }

            
        // make shearing constraints
            
        if (point_index % num_width_points > 0 && point_index >= num_width_points) {
            Spring diag_left_spring = Spring(&(this->point_masses[point_index]), &(this->point_masses[point_index - num_width_points - 1]), SHEARING);
            springs.emplace_back(diag_left_spring);
        }
            
        if (point_index % num_width_points < (num_width_points - 1) && point_index >= num_width_points) {
            Spring diag_right_spring = Spring(&(this->point_masses[point_index]), &(this->point_masses[point_index - num_width_points + 1]), SHEARING);
            springs.emplace_back(diag_right_spring);
        }
        
        
        // make bending constraints
        
        if (point_index % num_width_points >= 2) {
            Spring two_left_spring = Spring(&(this->point_masses[point_index]), &(this->point_masses[point_index - 2]), BENDING);
            springs.emplace_back(two_left_spring);
        }
        
        if (point_index >= 2 * num_width_points) {
            Spring two_up_spring = Spring(&(this->point_masses[point_index]), &(this->point_masses[point_index - 2 * num_width_points]), BENDING);
            springs.emplace_back(two_up_spring);
        }
        point_index += 1;
    }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2): Compute total force acting on each point mass.
    int i = 0;
    Vector3D total_force = Vector3D(0., 0., 0.);
    for (Vector3D accel : external_accelerations) {
        total_force += accel * mass;
    }
    for (PointMass &pm : this->point_masses) {
        pm.forces = total_force;
    }
    
    for (Spring &spring : this->springs) {
        double spring_force;
        if (spring.spring_type == STRUCTURAL) {
            if (cp->enable_structural_constraints) {
                double position_diff_norm = (spring.pm_a->position - spring.pm_b->position).norm();
                spring_force = cp->ks * (position_diff_norm - spring.rest_length);
            }
        }
        if (spring.spring_type == SHEARING) {
            if (cp->enable_shearing_constraints) {
                double position_diff_norm = (spring.pm_a->position - spring.pm_b->position).norm();
                spring_force = cp->ks * (position_diff_norm - spring.rest_length);
            }
        }
        if (spring.spring_type == BENDING) {
            if (cp->enable_bending_constraints) {
                double position_diff_norm = (spring.pm_a->position - spring.pm_b->position).norm();
                spring_force = 0.2 * cp->ks * (position_diff_norm - spring.rest_length);
            }
        }
        // apply force to one and equal and opposite to other
        spring.pm_a->forces -= spring_force * (spring.pm_a->position - spring.pm_b->position).unit();
        spring.pm_b->forces += spring_force * (spring.pm_a->position - spring.pm_b->position).unit();
    }
    

  // TODO (Part 2): Use Verlet integration to compute new point mass positions
    for (PointMass &pm : this->point_masses) {
        if (!pm.pinned) {
            Vector3D old_position = pm.position;
            pm.position = pm.position + (1. - (cp->damping / 100.)) * (pm.position - pm.last_position) + (pm.forces / mass) * pow(delta_t, 2.);
            pm.last_position = old_position;
        }
    }


  // TODO (Part 4): Handle self-collisions.


  // TODO (Part 3): Handle collisions with other primitives.
    build_spatial_map();
    for (PointMass &pm : point_masses) {
        self_collide(pm, simulation_steps);
    }
    for (auto c : *collision_objects) {
        for (PointMass &pm : this->point_masses) {
            c->collide(pm);
        }
    }

  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
    for (Spring &spring : this->springs) {
        double spring_length = (spring.pm_a->position - spring.pm_b->position).norm();
        // have to correct
        if (spring_length > (spring.rest_length * 1.1)) {
            // a pinned
            if (spring.pm_a->pinned) {
                Vector3D direction = (spring.pm_b->position - spring.pm_a->position).unit();
                spring.pm_b->position = direction * 1.1 * spring.rest_length + spring.pm_a->position;
            }
            // b pinned
            else if (spring.pm_b->pinned) {
                Vector3D direction = (spring.pm_a->position - spring.pm_b->position).unit();
                spring.pm_a->position = direction * 1.1 * spring.rest_length + spring.pm_b->position;
            }
            // neither pinned
            else {
                Vector3D midpoint = (spring.pm_a->position + spring.pm_b->position) / 2.;
                spring.pm_a->position = (spring.pm_a->position - midpoint).unit() * 1.1 * spring.rest_length / 2. + midpoint;
                spring.pm_b->position = (spring.pm_b->position - midpoint).unit() * 1.1 * spring.rest_length / 2. + midpoint;
            }
        }
    }

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
