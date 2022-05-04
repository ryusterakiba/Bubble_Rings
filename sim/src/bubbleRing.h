#ifndef BUBBLERING_H
#define BUBBLERING_H

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "clothMesh.h"
#include "spring.h"
#include "ringVertex.h"

using namespace CGL;
using namespace std;

struct BubbleRing {
	BubbleRing() {}
	BubbleRing(int num_vertices, double initial_ring_radius, Vector3D center, double thickness);
	~BubbleRing();

	void buildRing();
	void buildRingRender();
	void simulate(double frames_per_sec, double simulation_steps);
	void reset();

	void buildClothMesh();


	// Bubble ring physics functions
	void velocity();
	Vector3D biotsavart_edge(RingVertex p, RingVertex v0, RingVertex v1);
	Vector3D induction(RingVertex v0, RingVertex v1, RingVertex v2);
	Vector3D boussinesq(RingVertex v0, RingVertex v1);
	void modify_thickness();
	void resample();
	double volume();
	void burgers_flow(double delta_t);

	// Bubble ring properties
	double initial_ring_radius;
	int num_vertices;
	int num_render_circle;
	Vector3D center;
	double circulation;
	double thickness;

	// Properties for physics
	double vol0 = -1; // FIXME? Set this value in simulate()
	double min_dist; // set in constructor

	// Bubble ring components
	vector<RingVertex> vertices;

	// For rendering
	double width;
	double height;
	vector<PointMass> point_masses;
	ClothMesh* clothMesh;
};


#endif /* BUBBLERING_H */
