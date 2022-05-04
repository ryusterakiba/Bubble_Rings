#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

const double muRM = exp(-0.75);
const double dRM = 0.5 * exp(0.25);
const double nu = 1e-6;
const Vector3D g = Vector3D(0, 9.8, 0);
const double At = -1;

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
    this->min_dist = initial_ring_radius * angle;
    cout << "initial_ring_radius" << initial_ring_radius << endl;
    cout << "angle" << angle << endl;
    cout << "min_dist" << min_dist << endl;

    for (int i = 0; i < num_vertices; i++) {
        Vector3D position;
        position.x = initial_ring_radius * cos(angle * i);
        // printf("%f\n", initial_ring_radius * cos(angle * i));
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
  // This function performs 1 timestep of Runge Kutta
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // Remove thickness-rending points
  deleteExtraPoints();
  if (num_vertices != point_masses.size()) {
    cout << "deleteExtraPoints() did not work!" << endl;
    cout << "point_masses.size(): " << point_masses.size() << endl;
    cout << "num_vertices: " << num_vertices << endl << endl;
    throw 11; // throw error, arbitrary number idk
  }

  // Set vol0 if we haven't done so already
  if (this->vol0 == -1) {
    this->vol0 = volume();
  }

  // Use temp_position for Runge Kutta calculations
  for (int i = 0; i < point_masses.size(); i++) {
    PointMass& pm = point_masses.at(i);
    pm.temp_position = pm.position;
  }

  // Step 1
  velocity();
  vector<Vector3D> K1(point_masses.size());
  for (int i = 0; i < point_masses.size(); i++) {
    Vector3D& k1 = K1.at(i);
    Vector3D vel = point_masses.at(i).point_velocity;
    k1 = delta_t * vel;
  }

  // Step 2
  for (int i = 0; i < point_masses.size(); i++) {
    PointMass& pm = point_masses.at(i);
    pm.last_position = pm.temp_position;
    pm.temp_position += 0.5 * K1.at(i);
  }
  modify_thickness();
  velocity();

  vector<Vector3D> K2(point_masses.size());
  for (int i = 0; i < point_masses.size(); i++) {
    Vector3D& k2 = K2.at(i);
    Vector3D vel = point_masses.at(i).point_velocity;
    k2 = delta_t * vel;
  }

  // Step 3
  for (int i = 0; i < point_masses.size(); i++) {
    PointMass& pm = point_masses.at(i);
    pm.last_position = pm.temp_position;
    pm.temp_position += 0.5 * K2.at(i);
  }
  modify_thickness();
  velocity();
  
  vector<Vector3D> K3(point_masses.size());
  for (int i = 0; i < point_masses.size(); i++) {
    Vector3D& k3 = K3.at(i);
    Vector3D vel = point_masses.at(i).point_velocity;
    k3 = delta_t * vel;
  }

  // Step 4
  for (int i = 0; i < point_masses.size(); i++) {
    PointMass& pm = point_masses.at(i);
    pm.last_position = pm.temp_position;
    pm.temp_position += K3.at(i);
  }
  modify_thickness();
  velocity();
  
  vector<Vector3D> K4(point_masses.size());
  for (int i = 0; i < point_masses.size(); i++) {
    Vector3D& k4 = K4.at(i);
    Vector3D vel = point_masses.at(i).point_velocity;
    k4 = delta_t * vel;
  }

  // RK4
  for (int i = 0; i < point_masses.size(); i++) {
    PointMass& pm = point_masses.at(i);
    pm.position += (K1.at(i) + 2 * K2.at(i) + 2 * K3.at(i) + K4.at(i)) / 6;
  }

  // Volume conservation
  double vol1 = volume();
  for (auto pm_iter = point_masses.begin(); pm_iter != point_masses.end(); pm_iter++) {
    pm_iter->thickness = pm_iter->thickness * sqrt(vol0 / vol1);
  }

  // Resample number of points
  // resample();

  // Burgers thickness flow
  //burgers_flow(delta_t);

  // Add thickness-rendering points back in
  createExtraPoints();
}

/**
 * Update the point velocities.
 */
void Cloth::velocity() {
  // Initialize point velocities to 0 before calculations
  for (auto pm_iter = point_masses.begin(); pm_iter != point_masses.end(); pm_iter++) {
    pm_iter->point_velocity = Vector3D(0, 0, 0);
  }

  for (int i = 0; i < num_vertices; i++) {
    // Over all vertices of filament
    PointMass& pm = point_masses.at(i);

    for (int j = 0; j < num_vertices; j++) {
      // Biot Savart sum over edges (TERM 1 uBSdisc)
      int j_plus1 = (j + 1) % point_masses.size();
      pm.point_velocity += biotsavart_edge(pm, point_masses.at(j), point_masses.at(j_plus1));
    }

    // Apply induction term  (TERM 2 uLIA)
    int i_minus1 = (i + point_masses.size() - 1) % point_masses.size();
    int i_plus1 = (i + 1) % point_masses.size();
    pm.point_velocity += induction(point_masses.at(i_minus1), pm, point_masses.at(i_plus1));
  }

  // Boussinesq term (TERM 3) evaluated at edges then interpolated to vertices
  vector<Vector3D> edge_vel(num_vertices);
  for (int i = 0; i < num_vertices; i++) {
    // Boussinesq over edges
    int i_plus1 = (i + 1) % point_masses.size();
    Vector3D& ev = edge_vel.at(i);
    ev = boussinesq(point_masses.at(i), point_masses.at(i_plus1));
  }
  for (int i = 0; i < num_vertices; i++) {
    // interpolate to vertices
    PointMass &pm = point_masses.at(i);
    int i_minus1 = (i + point_masses.size() - 1) % point_masses.size();
    pm.point_velocity += (edge_vel.at(i) + edge_vel.at(i_minus1)) / 2;
  }
}

/**
 * Biosavart term for velocity calculations.
 */
Vector3D Cloth::biotsavart_edge(PointMass p, PointMass v0, PointMass v1) {
  double a = v0.thickness;
  double C = v0.circulation;

  Vector3D r0 = v0.temp_position - p.temp_position;
  Vector3D r1 = v1.temp_position - p.temp_position;
  Vector3D T = r1 - r0;
  double a2mu = pow(a * muRM, 2);
  Vector3D crossr01 = cross(r0, r1);

  double term1 = dot(r1, T) / (sqrt(a2mu + r1.norm2()) * (T.norm2() * a2mu + crossr01.norm2()));
  double term2 = dot(r0, T) / (sqrt(a2mu + r0.norm2()) * (T.norm2() * a2mu + crossr01.norm2()));
  return C / (4 * PI) * (term1 - term2) * crossr01;
}

/**
 * Induction term for velocity calculations.
 */
Vector3D Cloth::induction(PointMass v0, PointMass v1, PointMass v2) {
  double c0 = v0.circulation;
  double c1 = v1.circulation;
  double a0 = v0.thickness;
  double a1 = v1.thickness;

  double s0 = (v1.temp_position - v0.temp_position).norm();
  double s1 = (v2.temp_position - v1.temp_position).norm();
  Vector3D T0 = (v1.temp_position - v0.temp_position) / s0;
  Vector3D T1 = (v2.temp_position - v1.temp_position) / s1;
  double C = (c0 + c1) / 2;

  Vector3D kB = 2 * cross(T0, T1) / (s0 + s1);
  double logterm = log(s0 * s1 / (a0 * a1 * dRM * dRM));
  return C / (4 * PI) * .5 * logterm * kB;
}

/**
 * Boussinesq term for velocity calculations.
 */
Vector3D Cloth::boussinesq(PointMass v0, PointMass v1) {
  double C = v0.circulation;
  double a = v0.thickness;

  Vector3D edge = v1.temp_position - v0.temp_position;
  double ds = edge.norm();
  Vector3D T = edge / ds;
  Vector3D Atg_n = At * g - dot(At * g, T) * T;

  double coeff1 = 16 * PI*PI * nu * a*a / (256 * PI*PI * nu*nu + C*C);
  double coeff2 = PI * a*a * C / (256 * PI*PI * nu*nu + C*C);
  return coeff1 * Atg_n + coeff2 * cross(T, Atg_n);
}

/**
 * Updates the thickness at each edge (stored on PointMasses).
 * Uses each PointMass's temp_position and last_position.
 * Sets values for each PointMass's thickness.
 */
void Cloth::modify_thickness() {
  for (int i = 0; i < num_vertices; i++) {
    int i_plus1 = (i + 1) % num_vertices;
    PointMass& pm0 = point_masses.at(i);
    PointMass& pm1 = point_masses.at(i_plus1);

    double l_old = (pm1.last_position - pm0.last_position).norm();
    double l_new = (pm1.temp_position - pm0.temp_position).norm();
    pm0.thickness = pm0.thickness * sqrt(l_old / l_new);
  }
}

static double lerp(double x, double y0, double y1) {
  return y0 + x * (y1 - y0);
}

static Vector3D lerp(double x, Vector3D y0, Vector3D y1) {
  return y0 + x * (y1 - y0);
}

/**
 * Resamples the number of the points on the bubble ring.
 * Modifies positions, thicknesses, C values.
 * Also changes num_vertices of the Cloth.
 */
void Cloth::resample() {
  // Create list of cumulative sums
  vector<double> d;
  double d_sum = 0;
  d.push_back(d_sum);
  for (int i = 0; i < num_vertices; i++) {
    int i_plus1 = (i + 1) % num_vertices;
    d_sum += (point_masses.at(i_plus1).position - point_masses.at(i).position).norm();
    d.push_back(d_sum);
  }

  // Determine new number of vertices and distance between them
  int new_num_vertices = floor(d_sum / min_dist);
  if (new_num_vertices <= num_vertices)
    return;
  num_vertices = new_num_vertices;
  double space = d_sum / num_vertices;
  cout << "Resampling!" << endl;
  cout << "new_num_vertices=" << new_num_vertices << endl;
  cout << "d_sum=" << d_sum << endl;
  cout << "d:" << endl;
  for (auto it = d.begin(); it != d.end(); it++) {
    cout << *it << " | ";
  }
  cout << endl;

  // Create new points
  vector<PointMass> new_pms;
  double d_new = 0;
  int j = 0;
  for (int i = 0; i < num_vertices; i++) {
    while (d_new >= d.at(j + 1)) {
      j++;
    }
    int j_plus1 = (j + 1) % point_masses.size();
    PointMass& pm0 = point_masses.at(j);
    PointMass& pm1 = point_masses.at(j_plus1);
    /*cout << "j=" << j << " | ";
    cout << "d_new=" << d_new << " | " << endl;*/

    // Linearly interpolate positions, squared thicknesses, and circulation
    double x = (d_new - d.at(j)) / (d.at(j + 1) - d.at(j));
    Vector3D pos_new = lerp(x, pm0.position, pm1.position);
    //cout << "pos_new=" << pos_new << endl;
    double a2_new = lerp(x, pow(pm0.thickness, 2), pow(pm1.thickness, 2));
    double a_new = sqrt(a2_new);
    double C_new = lerp(x, pm0.circulation, pm1.circulation);

    // Create new point mass
    new_pms.emplace_back(pos_new, false, a_new, C_new);
    
    d_new += space;
  }
  cout << "point_masses.size()=" << point_masses.size() << endl;
  cout << "point_masses:" << endl;
  for (auto it = point_masses.begin(); it != point_masses.end(); it++) {
    cout << it->position << " | " << it->thickness << " | " << it->circulation << endl;
  }
  point_masses = new_pms;
  cout << "new point_masses.size()=" << point_masses.size() << endl;
  cout << "new point_masses:" << endl;
  for (auto it = point_masses.begin(); it != point_masses.end(); it++) {
    cout << it->position << " | " << it->thickness << " | " << it->circulation << endl;
  }
  cout << "Leaving resample()" << endl;
  //throw 11;
}

/**
 * Calculates the volume of the bubble ring.
 * Uses the positions and thicknesses stored in each PointMass.
 */
double Cloth::volume() {
  double vol = 0;

  for (int i = 0; i < num_vertices; i++) {
    int i_plus1 = (i + 1) % num_vertices;
    PointMass& pm0 = point_masses.at(i);
    PointMass& pm1 = point_masses.at(i_plus1);

    double d = (pm1.position - pm0.position).norm();
    vol += d * pm0.thickness * pm0.thickness;
  }

  vol *= PI;
  return vol;
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

void Cloth::deleteExtraPoints() {
  // Delete all point masses except the main outer ring
  for (int i = point_masses.size() - num_vertices - 1; i >= 0; i -= num_vertices + 1) {
    vector<PointMass>::iterator first = point_masses.begin() + i + 1;
    vector<PointMass>::iterator last = point_masses.begin() + i + num_vertices + 1;
    point_masses.erase(first, last);
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
    int main_ring_idx = i * (num_vertices + 1);

    double new_angle = 2 * PI / num_vertices;
    for (int j = 0; j < num_vertices; j++) {
        Vector3D position;
        PointMass pm = point_masses[main_ring_idx];
        position.y = pm.position.y + pm.thickness * sin(new_angle * j);
        position.x = pm.position.x + pm.thickness * cos(new_angle * j) * cos(angle * i);
        position.z = pm.position.z + pm.thickness * cos(new_angle * j) * sin(angle * i);
        PointMass new_pm = PointMass(position, false, 0.0);

        vector<PointMass>::iterator it = point_masses.begin() + main_ring_idx + 1 + j;
        point_masses.insert(it, new_pm);
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
