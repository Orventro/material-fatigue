#include <iostream>
#include <complex>
#include <cmath>
#include <random>
#include <fstream>
#include <set>

using namespace std;

using cmp = complex<double>;

const double CIRCLE_RAD = .5e-3, SQUARE_RAD = 12e-3;
const int M = 16, N = 80, R = 50, TOTAL = 2000;


void readVtk(istream &in, vector<cmp> &verts, vector<set<int>> &edges, vector<bool> &boundary) {
    string line;
    line.resize(50);
    do {
        getline(in, line);
    } while(line.substr(0, 4) != "elem");
    int n, m;
    in >> n;
    for(int i = 0; i < n; i++) {
        int a, b, el[4];
        in >> a >> b >> el[0] >> el[1] >> el[2] >> el[3];
        for(int j = 0; j < 4; j++) {
            if (edges.size() <= el[j]+1) 
                edges.resize(el[j]+1);
            edges[el[j]].insert(el[(j+1)%4]);
            edges[el[j]].insert(el[(j+3)%4]);
        }
    }
    do getline(in, line); while(line.substr(0, 4) != "boun");
    in >> n;
    boundary.resize(edges.size(), 0);
    for(int i = 0; i < n; i++) {
        int a, b, u, v;
        in >> a >> b >> u >> v;
        boundary[u] = 1;
        boundary[v] = 1;
    }
    do getline(in, line); while(line.substr(0, 4) != "vert");
    in >> n >> m;
    verts.resize(n);
    for(int i = 0; i < n; i++) {
        double x, y;
        in >> x >> y;
        verts[i] = {x, y};
    }
}

void adjustMesh(vector<cmp> &verts, vector<set<int>> &edges, vector<bool> &boundary) {
    vector<cmp> adj(verts.size(), 0);
    double alpha = 1e-1, threshold = .1e-7, oldmaxadj = 1;
    int lastk = 0, lastjump = 1;
    for(int k = 1; k < 100000; k++){
        for(int i = 0; i < verts.size(); i++) {
            if (boundary[i]) continue;
            adj[i] = 0;
            for (auto j : edges[i]) 
                adj[i] += (verts[j] - verts[i])*abs((verts[j] - verts[i]))*3e2;//*(k>1000?0.3:1.0);
        }
        double maxadj = 0;
        for(int i = 0; i < verts.size(); i++) {
            if (boundary[i]) continue;
            verts[i] += adj[i];
            maxadj = max(maxadj, abs(adj[i]));
        }
        if (maxadj*lastjump*9 <= threshold) {
            cout << maxadj*lastjump*9 << ' ' << lastjump << endl;
            break;
        }
        if (maxadj <= oldmaxadj*0.9) {
            cout << maxadj*lastjump*9 << ' ' << lastjump << endl;
            oldmaxadj = maxadj;
            lastjump = k - lastk;
            lastk = k;
        }
    }
}

void adjustVtk(istream &in, vector<cmp> &verts, ostream &out) {
    string line;
    do {
        getline(in, line);
        out << line << endl;
    } while(line.substr(0, 4) != "vert");
    out << verts.size() << "\n2\n";
    for(auto v : verts)
        out << v.real() << ' ' << v.imag() << endl;
}

int main(int argv, char *argc[]) {
    if (argv != 3) {
        cout << "Arguments:\n";
        cout << "\t1. Input mesh\n";
        cout << "\t2. Output mesh (!=input)\n";
        exit(1);
    }
    ifstream in(argc[1]);
    ofstream out(argc[2]);
    vector<cmp> verts;
    vector<bool> boundary;
    vector<set<int>> edges;
    readVtk(in, verts, edges, boundary);
    adjustMesh(verts, edges, boundary);
    in.seekg(0);
    adjustVtk(in, verts, out);
}
