#include <iostream>
#include <CDT.h>
#include <complex>
#include <cmath>
#include <random>
#include <fstream>

using namespace std;

using cmp = complex<double>;
using cdtv = CDT::V2d<double>;

const double CIRCLE_RAD = .5e-3, SQUARE_RAD = 12e-3;
const int M = 16, N = 80, R = 50, TOTAL = 2000;

cdtv polarToDecart(double theta, double radius) {
	cmp t = exp(cmp(0, theta)) * radius;
    return {t.real(), t.imag()};
}

double dist2(cdtv v1, cdtv v2) {
    return (v1.x - v2.x)*(v1.x - v2.x) + (v1.y - v2.y)*(v1.y - v2.y);
}

void printVtk(ostream &out, CDT::Triangulation<double> &cdt) {
    out << "MFEM mesh v1.0\ndimension\n2\n";

    out << "elements\n" << cdt.triangles.size() << "\n";
    for(auto t : cdt.triangles) 
        out << "1 2 " << t.vertices[0] << ' ' << t.vertices[1] << ' ' << t.vertices[2] << '\n';
    
    out << "boundary\n" << (M+8*N) << "\n";
    for(int i = 0; i < M-1; i++) out << "1 1 " << i << ' ' << (i+1) << '\n';  // circle edges
    out << "1 1 1 " << M-1 << '\n';

    for(int i = 0; i < N*8-1; i++) out << (i/(2*N)+2) << " 1 " << (i+M) << ' ' << (i+M+1) << '\n';  // square edges
    out << "5 1 " << M << ' ' << (M+8*N+1) << '\n';

    out << "vertices\n" << cdt.vertices.size() << "\n2\n";
    for(auto t : cdt.vertices) 
        out << t.x << ' ' << t.y << '\n';
}

int main(int argv, char *argc[]) {
    // if (argv !=)

    vector<cdtv> verts;
    for(int i = 0; i < M; i++) 
        verts.push_back(polarToDecart(i*2*M_PI/M, CIRCLE_RAD));  // inner circle

    for(int i = -N; i < N; i++) verts.push_back({ i*SQUARE_RAD/N,  SQUARE_RAD});
    for(int i = -N; i < N; i++) verts.push_back({ SQUARE_RAD, -i*SQUARE_RAD/N});
    for(int i = -N; i < N; i++) verts.push_back({-i*SQUARE_RAD/N, -SQUARE_RAD});
    for(int i = -N; i < N; i++) verts.push_back({-SQUARE_RAD,  i*SQUARE_RAD/N});


    vector<CDT::Edge> edges;
    for(int i = 0; i < M-1; i++) edges.push_back({i, i+1});
    edges.push_back({0, M-1});  // circle edges

    for(int i = 0; i < N*8-1; i++) edges.push_back({i+M, i+1+M});
    edges.push_back({M, M+8*N-1});  // square edges

    
    uniform_real_distribution<double> uniform(-SQUARE_RAD, SQUARE_RAD);
    mt19937_64 reng;
    for(int i = 0; i < TOTAL; i++) {  // generate uniform distribution of points in plate
        vector<cdtv> generated(R);
        double maxd2 = 1;
        int maxdi = 0;
        for(int j = 0; j < R; j++) {
            do {
                generated[j] = {uniform(reng), uniform(reng)};
            } while(dist2(generated[j], {0, 0}) < CIRCLE_RAD * CIRCLE_RAD * (1 + 0.2/M));
            double mind2 = 1e12;
            for(auto v2 : verts) 
                mind2 = min(mind2, dist2(generated[j], v2));
            if (maxd2 < mind2) {
                maxd2 = mind2;
                maxdi = j;
            }
        }
        verts.push_back(generated[maxdi]);
    }

    CDT::Triangulation<double> cdt;
    cdt.insertVertices(verts);
    cdt.insertEdges(edges);
    cdt.eraseOuterTrianglesAndHoles();
    
    ofstream out("../meshes/mesh_cpp_16_80.mesh");
    printVtk(out, cdt);
    out.close();
}
