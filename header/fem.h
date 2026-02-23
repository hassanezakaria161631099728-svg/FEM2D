#ifndef FEM_H
#define FEM_H
#include <vector>
#include <string>
using namespace std;

struct Node { int id; double x; double y; };
struct Element { int id; int n1; int n2;double E;double A;double I;double q;string load_type;double L;double c;double s;};

void computeGeometry(Element &e, const vector<Node> &nodes);
void localStiffness(const Element &e, double k[6][6]);
void transformationMatrix(const Element &e, double T[6][6]);
void multiply6(double A[6][6], double B[6][6], double C[6][6]);
void transpose6(double A[6][6], double AT[6][6]);
void globalStiffness(const Element &e, double kg[6][6]);
void equivalentLoad(const Element &e, double fe[6]);
void assembleGlobal(vector<vector<double>>& K,double kg[6][6],int n1, int n2);
void assembleLoad(vector<double>& F,double fe[6],int n1, int n2);
void applyBC(vector<vector<double>>& K,vector<double>& F,const vector<int>& bc);
void solveSystem(vector<vector<double>> K,vector<double> F,vector<double>& U);
void elementForces(const Element& e,const vector<double>& U,ofstream& out);
#endif
