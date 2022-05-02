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
        printf("%f\n", initial_ring_radius * cos(angle * i));
        position.z = initial_ring_radius * sin(angle * i);
        position.y = 0;
        PointMass pm = PointMass(position, false, 0.2);
        point_masses.emplace_back(pm);
        
        // make a circle around the pointmass
        double new_angle = 2 * PI / num_vertices;
        for (int j = 0; j < num_vertices; j++) {
            Vector3D cs_position;
            cs_position.y = pm.position.y + pm.thickness * sin(new_angle * j);
            cs_position.x = pm.position.x + pm.thickness * cos(new_angle * j) * cos(angle * i);
            cs_position.z = pm.position.z + pm.thickness * cos(new_angle * j) * sin(angle * i);
            PointMass new_pm = PointMass(cs_position, false, 0.0);
            point_masses.emplace_back(new_pm);
        }
        int point_index = 0;
        Spring left_spring = Spring(&(this->point_masses[point_masses.size() - num_vertices]), &(this->point_masses[point_masses.size() - 1]), STRUCTURAL);
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

void Cloth::deleteExtraPoints() {
    for (int i = point_masses.size() - num_vertices; i >= 0; i-=num_vertices+1) {
        point_masses.erase(point_masses.begin()+i);
    }
    
    // delete all springs, and just add back the outer ring ones
    
    springs.clear();
    
    int point_index = 0;
    
    Spring left_spring = Spring(&(this->point_masses[point_masses.size() - 1]), &(this->point_masses[0]), STRUCTURAL);
    springs.emplace_back(left_spring);

    while (point_index < num_vertices - 1) {
        left_spring = Spring(&(this->point_masses[point_index]), &(this->point_masses[(point_index + 1)]), STRUCTURAL);
        springs.emplace_back(left_spring);
        point_index += 1;
    }
}

void Cloth::createExtraPoints() {
    double angle = 2 * PI / num_vertices;
    for (int i = 0; i < num_vertices; i++) {
        double new_angle = 2 * PI / num_vertices;
        for (int j = 0; j < num_vertices; j++) {
            Vector3D position;
            PointMass pm = point_masses[i];
            position.y = pm.position.y + pm.thickness * sin(new_angle * j);
            position.x = pm.position.x + pm.thickness * cos(new_angle * j) * cos(angle * i);
            position.z = pm.position.z + pm.thickness * cos(new_angle * j) * sin(angle * i);
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

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;
    
    // generate triangles
    
    for (int i = 0; i < point_masses.size(); i++) {
        // if point is a ring center or last in ring
        if (i % (num_vertices + 1) == 0) {
            continue;
        }
        
        // pointmass is on the cross sectional ring; connect
        // it to adjacent neighbor and two neighbors on the next ring
        
        PointMass *pm = &point_masses[i];
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
        
        //not a big deal, just for textures
        float u_min = i;
        u_min /= num_vertices - 1;
        float u_max = i + 1;
        u_max /= num_vertices - 1;
        float v_min = i;
        v_min /= num_vertices - 1;
        float v_max = i + 1;
        v_max /= num_vertices - 1;
        
        PointMass *pm_A = pm;
        PointMass *pm_B = pm + 1;
        if (i + 1 >= point_masses.size()) {
            pm_B = pm + 1 - num_vertices;
        }
        
        PointMass *pm_C = pm + num_vertices + 1;
        if (i + num_vertices + 1 >= point_masses.size()) {
            pm_C = pm + num_vertices + 1 - point_masses.size();
        }
        
        PointMass *pm_D = pm_D = pm + num_vertices + 2;
        if (i + num_vertices + 2 >= point_masses.size()) {
            pm_D = pm + num_vertices + 2 - point_masses.size();
        }
        
        // specially account for last vertex in each cross-sectional ring, and need to also account for if last ring overall too
        if (i % (num_vertices + 1) == num_vertices) {
            pm_A = pm;
            pm_B = pm - num_vertices + 1;
            pm_C = pm + num_vertices + 1;
            if (i + num_vertices + 1 >= point_masses.size()) {
                pm_C = pm + num_vertices + 1 - point_masses.size();
            }
            pm_D = pm + 2;
            if (i + 2 >= point_masses.size()) {
                pm_D = pm + 2 - point_masses.size();
            }
        }
        
        
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
    
    int num_edges_in_cross_sectional_ring = num_vertices - 1;
  int num_tris_in_cross_sectional_ring = 2 * num_edges_in_cross_sectional_ring;


    bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, always exists
        int left_index = i - 1;
        if (left_index < 0) {
            left_index += num_tris_in_cross_sectional_ring;
        }
        Triangle *temp = triangles[left_index];
        t->pm1->halfedge->twin = temp->pm3->halfedge;

      // Get triangle above, always exists
        int above_index = i - num_tris_in_cross_sectional_ring + 1;
        if (above_index < 0) {
            above_index += triangles.size();
        }
        temp = triangles[above_index];
        t->pm3->halfedge->twin = temp->pm2->halfedge;

      // Get triangle to bottom right; guaranteed to exist
      temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
        int right_index = i + 1;
        if (i + 1 >= triangles.size()) {
            right_index -= num_tris_in_cross_sectional_ring;
        }
        Triangle *temp = triangles[right_index];
        t->pm3->halfedge->twin = temp->pm1->halfedge;

      // Get triangle below, always exists
        int below_index = i + num_tris_in_cross_sectional_ring - 1;
        if (below_index >= triangles.size()) {
            below_index -= triangles.size();
        }
        temp = triangles[below_index];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
        
        // Get triangle to top left; guaranteed to exist
        temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm2->halfedge;
      }
    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
