/*************************************************
*	Gyrification.h
*
*	Release: February 2016
*	Update: April 2016
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <vector>
#include "Mesh.h"

using std::vector;

class Gyrification
{
private:
	struct prop
	{
		double weight;
		double value;
		bool operator<(const prop& rhs) const { return value < rhs.value; }
	};
	struct tuple
	{
		int id;
		double value;
		bool operator<(const tuple& rhs) const { return value < rhs.value; }
	};
	// sulcal/gyral curves
	struct point
	{
		struct cpt
		{
			int id;		// vertex id
			int group;	// group
			int eid;	// id in group
			int type;	// curve type: 1 - sulcus, 2 - gyrus
			double dist;
			bool operator < (const point::cpt &p) const
			{
				return (dist < p.dist);
			}
		};	// struct for cand point
		int type;
		int id;	// vertex id
		int eid; // point id in curve
		int group;	// curve id
		int source;	// propagation source (from where?)
		vector<cpt> cands;	// candidate points
		double tan[3];		// tangent vector
		vector<double> GI;			// gyrification index
		vector<double> kArea1;		// kernel on input surface
		vector<double> kArea2;		// kernel on hull surface
		vector<int> boundary;		// vertex indecies on the kernel boundary
		double adjDist;		// adjusted distance
		double intpWeight;
		double intpSum;
		vector<prop> distrib;
		int state;		// is contained in a kernel?
	};	// struct for sulcal/gyral point
	vector< vector<point> > m_sulcus, m_gyrus;	// sulcal/gyral points for each curve
	vector<int> m_lookup;		// lookup table for correspondence
	point *m_vertex;
	vector<double> m_areamap1, m_areamap2;
	bool *m_check_work;		// work space for outer surface flag

	Mesh *m_mesh;
	Mesh *m_outer;
	
	double **T, *pT;
	double **m_u;
	double *m_pu;
	double *m_lam1, *m_lam2;
	double m_intv, m_maxArea, m_speed;
	double *m_medial_point;
	double m_populationArea;
	double m_adjRatio;
	
	char m_outfile[1024];
	
public:
	Gyrification(void);
	~Gyrification(void);
	void run(double rad);
	void run(const char *map);
	void open(const char *mesh, const char *sulcus, const char *gyrus, const char *outer, const char *corr);
	void open(const char *mesh, const char *outer, const char *corr);
	void open(const char *mesh, const char *gpoint, const char *outer, const char *corr);
	void setKernelInterval(double intv);
	void setMaxKernelSize(double area);
	void setSpeed(double speed);
	void saveGI(const char *filename);
	void saveCorrespondence(const char *filename);
	void saveGMap(const char *filename);
	void loadGMap(const char *filename);
	void saveKMap1(const char *filename);
	void saveKMap2(const char *filename);
	void saveKernelBoundary(const char *filename);
	void saveKernelBoundary(FILE *fp, int index, const double *area, double maxArea);
	void loadSpoint(const char *filename);
	void loadScurve(const char *filename);
	void loadGcurve(const char *filename);
	void setOutputFileName(const char *filename);
	void setPopulationArea(double populationArea);
private:
	void initVertex(void);
	void setupSDist(double distance);	// compute a geodesic distnace map for sulcus
	void setupGDist(double distance);	// compute a geodesic distnace map for gyrus
	void setupGDist2(double distance);	// compute a geodesic distnace map for gyrus
	void setupGDist2_old(double distance);	// compute a geodesic distnace map for gyrus
	void setupTensor(void);				// distance-based tensor
	void setupTensorEqualArea(void);	// distance-based tensor - constant eigen area
	void computeGyrification(void);		// gyrification index
	void precomputeArea(void);
	void smoothingTensor(int n);
	void smoothingScalar(double *scalar, int n);
	void medialPoints(void);
	void getKernelBoundary(int index, const double *area, double maxArea);
	double kernelArea(point *p, const double *dist, const int *state, const double *area = NULL, const double size = 0);			// kernel area
	double vertexArea(const Mesh *mesh, int id);	// vertex area
	double totalArea(void);
};
