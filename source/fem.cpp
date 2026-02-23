#include "fem.h"
#include <cmath>
#include <fstream>
#include <iomanip>

// ===================================================== //                 BASIC DATA STRUCTURES // =====================================================
// nodes elements
// ===================================================== //            PREPROCESSING // =====================================================
void computeGeometry(Element &e, const vector<Node> &nodes) {

double x1 = nodes[e.n1-1].x;
double y1 = nodes[e.n1-1].y;

double x2 = nodes[e.n2-1].x;
double y2 = nodes[e.n2-1].y;

double dx = x2 - x1;
double dy = y2 - y1;

e.L = sqrt(dx*dx + dy*dy);

e.c = dx / e.L;
e.s = dy / e.L;

}
// =====================================================
//ELEMENT ROUTINES // =====================================================
void localStiffness(const Element &e, double k[6][6]) {

double E = e.E;
double A = e.A;
double I = e.I;
double L = e.L;

double EA = E*A / L;
double EI = E*I;

double a = EA;
double b = 12*EI/(L*L*L);
double c = 6*EI/(L*L);
double d = 4*EI/L;
double f = 2*EI/L;

for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
        k[i][j]=0;

k[0][0]= a;   k[0][3]= -a;
k[3][0]= -a;  k[3][3]= a;

k[1][1]= b;   k[1][2]= c;   k[1][4]= -b;  k[1][5]= c;

k[2][1]= c;   k[2][2]= d;   k[2][4]= -c;  k[2][5]= f;

k[4][1]= -b;  k[4][2]= -c;  k[4][4]= b;   k[4][5]= -c;

k[5][1]= c;   k[5][2]= f;   k[5][4]= -c;  k[5][5]= d;

}

void transformationMatrix(const Element &e, double T[6][6]) {

double c = e.c;
double s = e.s;

for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
        T[i][j]=0;

T[0][0]= c;  T[0][1]= s;
T[1][0]= -s; T[1][1]= c;
T[2][2]= 1;

T[3][3]= c;  T[3][4]= s;
T[4][3]= -s; T[4][4]= c;
T[5][5]= 1;

}
void multiply6(double A[6][6], double B[6][6], double C[6][6]) {

for(int i=0;i<6;i++){
    for(int j=0;j<6;j++){

        C[i][j]=0;

        for(int k=0;k<6;k++)
            C[i][j]+=A[i][k]*B[k][j];
    }
}

}
void transpose6(double A[6][6], double AT[6][6]) {

for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
        AT[j][i]=A[i][j];

}

void globalStiffness(const Element &e, double kg[6][6]) {

double k[6][6];
double T[6][6];
double TT[6][6];
double temp[6][6];

localStiffness(e,k);
transformationMatrix(e,T);
transpose6(T,TT);

multiply6(TT,k,temp);
multiply6(temp,T,kg);

}

void equivalentLoad(const Element &e, double fe[6])
{
    double L = e.L;
    double q = e.q;

    double c = e.c;   // cos(theta)
    double s = e.s;   // sin(theta)

    double qx_local = 0.0;
    double qy_local = 0.0;

    if (e.load_type == "local") {
        // load already in local y direction
        qy_local = q;
    }
    else if (e.load_type == "global") {
        // load is vertical global load (0 , q)
        // convert to local coordinates
        qx_local =  s * q;
        qy_local =  c * q;
    }

    // --- Equivalent nodal loads in LOCAL system ---
    double fl[6] = {0};

    // Axial contribution from qx_local
    fl[0] = qx_local * L / 2.0;
    fl[3] = qx_local * L / 2.0;

    // Transverse contribution from qy_local
    fl[1] = qy_local * L / 2.0;
    fl[2] = qy_local * L * L / 12.0;

    fl[4] = qy_local * L / 2.0;
    fl[5] = -qy_local * L * L / 12.0;

    // --- Transform to GLOBAL coordinates ---
    double T[6][6];
    transformationMatrix(e, T);

    for(int i = 0; i < 6; i++){
        fe[i] = 0.0;
        for(int j = 0; j < 6; j++)
            fe[i] += T[j][i] * fl[j];   // T^T * fl
    }
}

// =====================================================
//GLOBAL SYSTEM ROUTINES // =====================================================

void assembleGlobal(vector<vector<double>>& K,double kg[6][6],int n1, int n2){

int map[6] = {
    3*(n1-1), 3*(n1-1)+1, 3*(n1-1)+2,
    3*(n2-1), 3*(n2-1)+1, 3*(n2-1)+2
};

for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
        K[ map[i] ][ map[j] ] += kg[i][j];

}
void assembleLoad(vector<double>& F,double fe[6],int n1, int n2) {

int map[6] = {
    3*(n1-1), 3*(n1-1)+1, 3*(n1-1)+2,
    3*(n2-1), 3*(n2-1)+1, 3*(n2-1)+2
};

for(int i=0;i<6;i++)
    F[ map[i] ] += fe[i];

}

void applyBC(vector<vector<double>>& K,vector<double>& F,const vector<int>& bc){

double big = 1e20;

for(size_t i=0;i<bc.size();i++){

    int d = bc[i]-1;

    K[d][d] += big;
    F[d] = 0;
}

}

void solveSystem(vector<vector<double>> K,vector<double> F,vector<double>& U){

int n=15;

double A[15][16];

for(int i=0;i<n;i++){
    for(int j=0;j<n;j++)
        A[i][j]=K[i][j];

    A[i][n]=F[i];
}

for(int k=0;k<n;k++){

    double piv=A[k][k];

    for(int j=k;j<=n;j++)
        A[k][j]/=piv;

    for(int i=0;i<n;i++){
        if(i==k) continue;

        double f=A[i][k];

        for(int j=k;j<=n;j++)
            A[i][j]-=f*A[k][j];
    }
}

for(int i=0;i<n;i++)
    U[i]=A[i][n];

}

// =====================================================
//STEP 4 : POST PROCESSING // =====================================================

void elementForces(const Element& e,const vector<double>& U,ofstream& out){

// 1. Extract element displacement vector in global

int map[6] = {
    3*(e.n1-1), 3*(e.n1-1)+1, 3*(e.n1-1)+2,
    3*(e.n2-1), 3*(e.n2-1)+1, 3*(e.n2-1)+2
};

double ug[6];

for(int i=0;i<6;i++)
    ug[i]=U[ map[i] ];

// 2. Transform to local

double T[6][6];
transformationMatrix(e,T);

double ul[6]={0};

for(int i=0;i<6;i++)
    for(int j=0;j<6;j++)
        ul[i]+=T[i][j]*ug[j];

// 3. Local stiffness

double k[6][6];
localStiffness(e,k);

// 4. Fixed end load

double fe[6];
equivalentLoad(e,fe);

// 5. Internal forces

double f_int[6]={0};

for(int i=0;i<6;i++){
    for(int j=0;j<6;j++)
        f_int[i]+=k[i][j]*ul[j];

    f_int[i]-=fe[i];
}

double  qx_local =  e.s * e.q;

double N1 = f_int[0]- (qx_local * e.L)/2;

double N2 = f_int[3]- (qx_local * e.L)/2;

// Print two rows per element (start and end node)

out << setw(9)  << e.id << " "
    << setw(9)  << e.n1 << " "
    << setw(14) << -N1 << " "
    << setw(14) << -f_int[1] << " "
    << setw(14) << -f_int[2] << "\n";

out << setw(9)  << e.id << " "
    << setw(9)  << e.n2 << " "
    << setw(14) << N2 << " "
    << setw(14) << f_int[4] << " "
    << setw(14) << f_int[5] << "\n";}


