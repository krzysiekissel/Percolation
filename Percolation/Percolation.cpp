#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <map>
#include <ctime>
using namespace std;

struct Params {
    int L;
    int T;
    float p0;
    float pk;
    float dp;
};
Params getParams(string inputFile) {
    Params result = Params();
    ifstream file = ifstream();
    file.open(inputFile);
    if (file.is_open()) {
        string buffer;
        getline(file, buffer);
        buffer = buffer.substr(0, buffer.find("%", 0));
        result.L = stoi(buffer);
        getline(file, buffer);
        buffer = buffer.substr(0, buffer.find("%", 0));
        result.T = stoi(buffer);
        getline(file, buffer);
        buffer = buffer.substr(0, buffer.find("%", 0));
        result.p0 = stof(buffer);
        getline(file, buffer);
        buffer = buffer.substr(0, buffer.find("%", 0));
        result.pk = stof(buffer);
        getline(file, buffer);
        buffer = buffer.substr(0, buffer.find("%", 0));
        result.dp = stof(buffer);
    }
    file.close();

    return result;
}
void printParams(Params p) {
    cout << "# L: " << p.L << endl;
    cout << "# T: " << p.T << endl;
    cout << "# p0: " << p.p0 << endl;
    cout << "# pk: " << p.pk << endl;
    cout << "# dp: " << p.dp << endl;
}
void printConfig(int** a, int L) {
    for (int row = 0; row < L; row++) {
        for (int col = 0; col < L; col++) {
            cout << a[row][col] << "  ";
        }
        cout << endl;
    }
    cout << endl;
}
bool isOccupied(int** a,int L,int row,int col) {
    if (row < 0 || row >= L) {
        return false;
    }
    if (col < 0 || col >= L) {
        return false;
    }
    return a[row][col] != 0;
}
void unionClusters(int** a, int L, int sm, int lg) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (a[i][j] == sm) {
                a[i][j] = lg;
            }
        }
    }
}
map<int,int> findClusters(int** a, int L) {
    int k = 2;
    map<int, int> clusters;
    for (int row = 0; row < L; row++) {
        for (int col = 0; col < L; col++) {
            if (a[row][col]) {
                //top and left empty
                if (!isOccupied(a, L, row - 1, col) && !isOccupied(a, L, row, col - 1)) {
                    k++;
                    a[row][col] = k;
                    clusters.insert({ k, 1 });
                }
                //one is occupied with k0
                if ((isOccupied(a, L, row - 1, col) && !isOccupied(a, L, row, col - 1)) || (!isOccupied(a, L, row - 1, col) && isOccupied(a, L, row, col - 1))) {
                    if ((isOccupied(a, L, row - 1, col))) {
                        a[row][col] = a[row-1][col];
                        if (clusters[a[row - 1][col]] >= 0) {
                            clusters[a[row - 1][col]]++;
                        }
                        
                    }
                    else {
                        a[row][col] = a[row][col-1];
                        if (clusters[a[row][col - 1]] >= 0) {
                            clusters[a[row][col - 1]]++;
                        }
                        
                    }
                    
                }
                //top and left occupied 
                if (isOccupied(a, L, row - 1, col) && isOccupied(a, L, row, col - 1)) {
                    if (a[row - 1][col] == a[row][col - 1]) {
                        a[row][col] = a[row - 1][col];
                        if (clusters[a[row - 1][col]] >= 0) {
                            clusters[a[row - 1][col]]++;
                        }
                    }
                    //top and left occupied but k1!=k2
                    else {
                        int sm = 0;
                        int lg = 0;
                        if (clusters[a[row - 1][col]] > clusters[a[row][col - 1]]) {
                            sm = a[row][col - 1];
                            lg = a[row - 1][col];
                        }
                        else {
                            sm = a[row-1][col];
                            lg = a[row][col-1];
                        }
                        a[row][col] = lg;
                        clusters[lg] = clusters[lg]+ clusters[sm]+1;
                        clusters.erase(sm);
                        unionClusters(a, L, sm, lg);
                    }
                }
            }
        }
    }
    return clusters;


}


bool checkPath(int** a,int L) {
    //init burn side
    for (int i = 0; i < L; i++) {
        if (a[0][i] != 0) {
            a[0][i] = 2;
        }
    }

    //burn process
    bool burnFlag = true;
    bool endFlag = false;
    int burnLevel = 2;
    while (burnFlag) {
        burnFlag = false;
        for (int row = 0; row < L; row++) {
            for (int col = 0; col < L; col++) {
                if (a[row][col] == burnLevel) {
                    burnFlag = true;
                    //n1 row-1 col
                    if (row - 1 >= 0) {
                        if (a[row - 1][col] == 1) {
                            a[row - 1][col] = burnLevel + 1;
                        }
                    }
                    //n2 row col-1
                    if (col - 1 >= 0) {
                        if (a[row][col - 1] == 1) {
                            a[row][col - 1] = burnLevel + 1;
                        }
                    }
                    //n3 row col+1
                    if (col + 1 < L) {
                        if (a[row][col + 1] == 1) {
                            a[row][col + 1] = burnLevel + 1;
                        }
                    }
                    //n4 row+1 col
                    if (row + 1 < L) {
                        if (a[row + 1][col] == 1) {
                            a[row + 1][col] = burnLevel + 1;
                            if (row + 1 == L - 1) {
                                goto stop;
                            }
                        }
                    }
                }
            }
        }
        burnLevel++;
    }
    return false;
    stop:
    return true;
}

int max(map<int, int> clusters) {
    int max = 0;
    for (auto const& x : clusters) {
        if (x.second > max) {
            max = x.second;
        }
    }
    return max;
}
int main(int argc, char **argv)
{
    auto params = getParams(string(argv[1]));
    printParams(params);
    const int L = params.L;
    const int T = params.T;
    const float p0 = params.p0;
    const float pk = params.pk;
    const float dp = params.dp;

    //latice
    int** a = new int* [L];
    int** b = new int* [L];
    for (int i = 0; i < L; i++) {
        a[i] = new int[L];
        b[i] = new int[L];
    }

    //random
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    string fileName = "Ave_L_" + to_string(L) + "_T_" + to_string(T) + "_tmstmp_" + to_string(time(nullptr)) + ".txt";
    fstream mainOutput;
    mainOutput.open(fileName,ios::out);
    if (!mainOutput.is_open()) {
        cout << "file error" << endl;
        return -1;
    }

    //p loop
    float p = p0;
    float epsilon = 0.0003;
    while (p <= (pk+epsilon)) {
        //mc loop
        float pflow = 0;
        float s = 0;
        map<int, float> clustersDistibution;

        for (int mcs = 0; mcs < T; mcs++) {
            //lattice loop
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    float rand = distribution(generator);
                    if (rand < p) {
                        a[i][j] = 1;
                        b[i][j] = 1;
                    }
                    else {
                        a[i][j] = 0;
                        b[i][j] = 0;
                    }
                }
            }
            //printConfig(a, L);
            auto clusters = findClusters(a, L);
            int max = 0;
            for (auto const& x : clusters) {
                clustersDistibution[x.second] += 1;
                if (x.second > max) {
                    max = x.second;
                }
            }
            s += max;
            bool path = checkPath(b, L);
            if (path) {
                pflow += 1;
            }

        }
        for (auto const& x : clustersDistibution) {
            clustersDistibution[x.first] /= (float)T;
        }
        pflow = pflow / T;
        s = s / T;
        string fileName = "Dist_p_" + to_string(p) + "_L_" + to_string(L)+"_T_"+to_string(T)+ "_tmstmp_" + to_string(time(nullptr)) + ".txt";
        fstream clustersOutput;
        clustersOutput.open(fileName, ios::out);
        for (auto const& x : clustersDistibution) {
            cout << x.first << "\t" << x.second << endl;
            clustersOutput << x.first << "\t" << x.second << endl;
        }
        clustersOutput.close();
        cout << p << "\t" << pflow << "\t" << s << endl;
        mainOutput << p << "\t"<<pflow<<"\t"<< s << endl;
        p += dp;
    }
    mainOutput.close();

    for (int i = 0; i < L; i++) {
        delete[] a[i];
        delete[] b[i];
    }
    delete[] a;
    delete[] b;
}





