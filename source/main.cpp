// 2D Frame FEM Solver - STEP 4: Post Processing // Pure C++ implementation without external libraries // 3 DOF per node: ux, uy, rotz
#include "fem.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

int main() {
ofstream out("results2.txt");
 if(!out) return 1;

vector<Node> nodes = {
    {1,0,0},
    {2,0,7},
    {3,9,9},
    {4,18,7},
    {5,18,0}
};

double E=2.1e8;
double A1=131.4e-4;
double I1=19270e-8;
double A2=84.5e-4;
double I2=23130e-8;
double q=-195;

vector<Element> elements;

elements.push_back({1,1,2,E,A1,I1,0,"local"});
elements.push_back({2,2,3,E,A2,I2,q,"global"});
elements.push_back({3,3,4,E,A2,I2,q,"global"});
elements.push_back({4,4,5,E,A1,I1,0,"local"});

for(size_t i=0;i<elements.size();i++)
    computeGeometry(elements[i],nodes);

int totalDOF = 15;
vector<vector<double>> K(totalDOF,vector<double>(totalDOF, 0.0));
vector<double> F(totalDOF, 0.0);

for(size_t i=0;i<elements.size();i++){

    double kg[6][6];
    double fe[6];

    globalStiffness(elements[i],kg);
    equivalentLoad(elements[i],fe);

    assembleGlobal(K,kg,elements[i].n1,elements[i].n2);
    assembleLoad(F,fe,elements[i].n1,elements[i].n2);
}
vector<vector<double>> K_original(totalDOF,vector<double>(totalDOF, 0.0));
vector<double> F_original(totalDOF, 0.0);


for(int i=0;i<15;i++){
    F_original[i] = F[i];
    for(int j=0;j<15;j++)
        K_original[i][j] = K[i][j];
}
// nodal loads
F[4]+= -1387;
F[10]+= -1387;

vector<int> bc = {1,2,3,13,14,15};

applyBC(K,F,bc);

vector<double> U(totalDOF, 0.0);

solveSystem(K,F,U);

vector<double> R(totalDOF, 0.0);

for(int i=0;i<15;i++){
    for(int j=0;j<15;j++)
        R[i] += K_original[i][j] * U[j];

    R[i] -= F_original[i];
}

out << "\n=========================================\n";
out << "        NODAL DISPLACEMENTS [cm]\n";
out << "=========================================\n\n";

out << setw(9)  << "Node"
    << setw(18) << "ux"
    << setw(18) << "uy"
    << setw(18) << "rotz" << "\n";

out << "--------------------------------------------------------------\n";
int numNodes=5;
for (int i = 0; i < numNodes; i++)
{
    double ux   = U[3*i + 0];
    double uy   = U[3*i + 1];
    double rotz = U[3*i + 2];

    out << setw(9)  << i+1
        << setw(18) << ux
        << setw(18) << uy
        << setw(18) << rotz
        << "\n";
}
out << "\n=========================================\n";
out << "        SUPPORT REACTIONS [kg & m]\n";
out << "=========================================\n\n";

out << setw(9)  << "Node"
    << setw(18) << "Rx"
    << setw(18) << "Ry"
    << setw(18) << "Mz" << "\n";

out << "--------------------------------------------------------------\n";

vector<int> fixedNodes = {1,5};

for(int node : fixedNodes){

    int base = 3*(node-1);

    out << setw(9)  << node
        << setw(18) << R[base+0]
        << setw(18) << R[base+1]
        << setw(18) << R[base+2]
        << "\n";
}
// -------- POST PROCESSING --------
out << "\nELEMENT FORCES (LOCAL AXIS) [kg & m]\n";
out << "--------------------------------------------------------------\n";
out << setw(9)  << "Element"
    << setw(9)  << "Node"
    << setw(14) << "N"
    << setw(14) << "V"
    << setw(14) << "M" << "\n";
out << "--------------------------------------------------------------\n";
for(size_t i=0;i<elements.size();i++)
    elementForces(elements[i],U,out);
 out.close();
return 0;

}


