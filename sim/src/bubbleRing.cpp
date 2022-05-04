#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "bubbleRing.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

const double muRM = exp(-0.75);
const double dRM = 0.5 * exp(0.25);
const double nu = 1e-6;
const Vector3D g = Vector3D(0, 0, -9.8);
const double At = -1;

BubbleRing::BubbleRing(int num_vertices, double initial_ring_radius, Vector3D center, double thickness) {
    this->num_vertices = num_vertices;
    this->initial_ring_radius = initial_ring_radius;
    this->center = center;
    this->num_render_circle = 20; //Resolution of mesh circles

    this->width = 2.0 * initial_ring_radius;
    this->height = 2.0 * thickness;

    buildRing();
    buildRingRender();
    buildClothMesh();

    this->vol0 = volume();
    double angle = 2 * PI / num_vertices;
    this->min_dist = initial_ring_radius * angle;
}

BubbleRing::~BubbleRing() {
    point_masses.clear();
    vertices.clear();

    if (clothMesh) {
        delete clothMesh;
    }
}

void BubbleRing::buildRing() {
    // Build the bubble ring initial
    double angle = 2 * PI / num_vertices;
    for (int i = 0; i < num_vertices; i++) {
        Vector3D position;
        position.x = initial_ring_radius * cos(angle * i) + center.x;
        position.y = initial_ring_radius * sin(angle * i) + center.y;
        position.z = center.z;
        vertices.emplace_back(position, thickness);
    }
}

void BubbleRing::buildRingRender() {
    // Build point masses for rendering
    point_masses.clear();
        double angle = 2 * PI / num_vertices;
    for (int i = 0; i < num_vertices; i++) {
        RingVertex v = vertices[i];
        point_masses.emplace_back(v.position, false, thickness);

        double new_angle = 2 * PI / num_render_circle;
        for (int j = 0; j < num_render_circle; j++) {
            Vector3D position;
            position.x = v.position.x + v.thickness * cos(new_angle * j) * cos(angle * i);
            position.y = v.position.y + v.thickness * cos(new_angle * j) * sin(angle * i);
            position.z = v.position.x + v.thickness * sin(new_angle * j);
            point_masses.emplace_back(position, false, 0.0);
        }
    }
}

/**
 * Main Simulation Function
 */
void BubbleRing::simulate(double frames_per_sec, double simulation_steps) {
    // Bubble Ring integrator using Runge Kutta 4th order
    // Use temp_position for Runge Kutta calculations
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        it->temp_position = it->position;
    }

    // Step 1
    velocity();
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        it->K1 = delta_t * it->velocity;
        it->last_position = it->temp_position;
        it->temp_position += 0.5 * it->K1;
    }
    modify_thickness();

    // Step 2
    velocity();
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        it->K2 = delta_t * it->velocity;
        it->last_position = it->temp_position;
        it->temp_position += 0.5 * it->K2;
    }
    modify_thickness();

    // Step 3
    velocity();
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        it->K3 = delta_t * it->velocity;
        it->last_position = it->temp_position;
        it->temp_position += it->K3;
    }
    modify_thickness();

    // Step 4
    velocity();
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        it->K4 = delta_t * it->velocity;
    }

    // RK4 position update
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        it->position += (it->K1 + 2.0 * it->K2 + 2.0 * it->K3 + it->K4) / 6.0;
    }

    // Volume conservation
    double vol1 = volume();
    double vol_scale = sqrt(vol0 / vol1);
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        it->thickness = it->thickness * vol_scale;
    }

    // Resample number of points
    //resample();

    // Burgers thickness flow
    //burgers_flow(delta_t);

    // Add in the point masses
    buildRingRender();
}


/**
 * Update the point velocities.
 */
void BubbleRing::velocity() {
    // Initialize point velocities to 0 before calculations
    for (auto it = vertices.begin(); it != vertices.end(); it++) {
        it->velocity = Vector3D(0, 0, 0);
    }

    for (int i = 0; i < num_vertices; i++) {
        // Over all vertices of filament
        RingVertex& vi = vertices.at(i);

        for (int j = 0; j < num_vertices; j++) {
            // Biot Savart sum over edges (TERM 1 uBSdisc)
            int j_plus1 = (j + 1) % vertices.size();
            vi.velocity += biotsavart_edge(vi, vertices.at(j), vertices.at(j_plus1));
        }

        // Apply induction term  (TERM 2 uLIA)
        int i_minus1 = (i + vertices.size() - 1) % vertices.size();
        int i_plus1 = (i + 1) % vertices.size();
        vi.velocity += induction(vertices.at(i_minus1), vi, vertices.at(i_plus1));
    }

    // Boussinesq term (TERM 3) evaluated at edges then interpolated to vertices
    vector<Vector3D> edge_vel(num_vertices);
    for (int i = 0; i < num_vertices; i++) {
        // Boussinesq over edges
        int i_plus1 = (i + 1) % vertices.size();
        Vector3D& ev = edge_vel.at(i);
        ev = boussinesq(vertices.at(i), vertices.at(i_plus1));
    }
    for (int i = 0; i < num_vertices; i++) {
        // interpolate to vertices
        RingVertex& vi = vertices.at(i);
        PointMass& pm = point_masses.at(i);
        int i_minus1 = (i + vertices.size() - 1) % vertices.size();
        vi.velocity += (edge_vel.at(i) + edge_vel.at(i_minus1)) / 2;
    }
}

/**
 * Biosavart term for velocity calculations.
 */
Vector3D BubbleRing::biotsavart_edge(RingVertex p, RingVertex v0, RingVertex v1) {
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
Vector3D BubbleRing::induction(RingVertex v0, RingVertex v1, RingVertex v2) {
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
Vector3D BubbleRing::boussinesq(RingVertex v0, RingVertex v1) {
    double C = v0.circulation;
    double a = v0.thickness;

    Vector3D edge = v1.temp_position - v0.temp_position;
    double ds = edge.norm();
    Vector3D T = edge / ds;
    Vector3D Atg_n = At * g - dot(At * g, T) * T;

    double coeff1 = 16 * PI * PI * nu * a * a / (256 * PI * PI * nu * nu + C * C);
    double coeff2 = PI * a * a * C / (256 * PI * PI * nu * nu + C * C);
    return coeff1 * Atg_n + coeff2 * cross(T, Atg_n);
}

/**
 * Updates the thickness at each edge (stored on PointMasses).
 * Uses each PointMass's temp_position and last_position.
 * Sets values for each PointMass's thickness.
 */
void BubbleRing::modify_thickness() {
    for (int i = 0; i < num_vertices; i++) {
        int i_plus1 = (i + 1) % num_vertices;
        RingVertex& v0 = vertices.at(i);
        RingVertex& v1 = vertices.at(i_plus1);

        double l_old = (v1.last_position - v0.last_position).norm();
        double l_new = (v1.temp_position - v0.temp_position).norm();
        v0.thickness = v0.thickness * sqrt(l_old / l_new);
    }
}

/**
 * Resamples the number of the points on the bubble ring.
 * Modifies positions, thicknesses, C values.
 * Also changes num_vertices of the Cloth.
 */
void BubbleRing::resample() {

}

/**
 * Calculates the volume of the bubble ring.
 * Uses the positions and thicknesses stored in each RingVertex.
 */
double BubbleRing::volume() {
    double vol = 0;
    for (int i = 0; i < num_vertices; i++) {
        int i_plus1 = (i + 1) % num_vertices;
        RingVertex& v0 = vertices.at(i);
        RingVertex& v1 = vertices.at(i_plus1);

        double d = (v1.position - v0.position).norm();
        vol += d * v0.thickness * v0.thickness;
    }
    vol *= PI;
    return vol;
}

/**
 * Modifies thicknesses according Burgers flow.
 */
void BubbleRing::burgers_flow(double delta_t) {

}

void BubbleRing::reset() {
    RingVertex* v = &vertices[0];
    PointMass* pm = &point_masses[0];
    for (int i = 0; i < point_masses.size(); i++) {
        pm->position = pm->start_position;
        pm->last_position = pm->start_position;
        pm++;
    }
    for (int i = 0; i < vertices.size(); i++) {
        v->position = v->start_position;
        v->last_position = v->start_position;
        v++;
    }
}

void BubbleRing::buildClothMesh() {
    if (point_masses.size() == 0) return;

    ClothMesh* clothMesh = new ClothMesh();
    vector<Triangle*> triangles;

    // generate triangles

    for (int i = 0; i < point_masses.size(); i++) {
        // if point is a ring center or last in ring
        if (i % (num_vertices + 1) == 0) {
            continue;
        }

        // pointmass is on the cross sectional ring; connect
        // it to adjacent neighbor and two neighbors on the next ring

        PointMass* pm = &point_masses[i];
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

        PointMass* pm_A = pm;
        PointMass* pm_B = pm + 1;
        if (i + 1 >= point_masses.size()) {
            pm_B = pm + 1 - num_vertices;
        }

        PointMass* pm_C = pm + num_vertices + 1;
        if (i + num_vertices + 1 >= point_masses.size()) {
            pm_C = pm + num_vertices + 1 - point_masses.size();
        }

        PointMass* pm_D = pm_D = pm + num_vertices + 2;
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
        Triangle* t = triangles[i];

        // Allocate new halfedges on heap
        Halfedge* h1 = new Halfedge();
        Halfedge* h2 = new Halfedge();
        Halfedge* h3 = new Halfedge();

        // Allocate new edges on heap
        Edge* e1 = new Edge();
        Edge* e2 = new Edge();
        Edge* e3 = new Edge();

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
        Triangle* t = triangles[i];

        if (topLeft) {
            // Get left triangle, always exists
            int left_index = i - 1;
            if (left_index < 0) {
                left_index += num_tris_in_cross_sectional_ring;
            }
            Triangle* temp = triangles[left_index];
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
        }
        else {
            // Get right triangle, if it exists
            int right_index = i + 1;
            if (i + 1 >= triangles.size()) {
                right_index -= num_tris_in_cross_sectional_ring;
            }
            Triangle* temp = triangles[right_index];
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


