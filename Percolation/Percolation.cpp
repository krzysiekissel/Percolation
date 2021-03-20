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
enum NN {
    TOPLEFT,
    TOP,
    LEFT,
    RIGHT,
    BOTTOM,
    BOTTOMRIGHT
};
struct Node {
    int val;
    Node* nn[6];
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
void printConfig(Node* a, int L) {
    int LL = L * L;
    for (int i = 0; i < LL; i++) {

        cout << a[i].val << "  ";

        if (i % L == L-1) {
            cout << endl;
        }
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
void unionClusters(Node* a, int L, int sm, int lg) {
    int LL = L * L;
    for (int i = 0; i < LL; i++) {
        if (a[i].val == sm) {
            a[i].val = lg;
        }
    }
}
map<int,int> findClusters(Node* a, int L) {
    int k = 2;
    int LL = L * L;
    map<int, int> clusters;
    for (int i = 0; i < LL; i++) {
        if (a[i].val) {
            if (!a[i].nn[NN::TOPLEFT]->val && !a[i].nn[NN::TOP]->val && !a[i].nn[NN::LEFT]->val) {
                k++;
                a[i].val = k;
                clusters[k] = 1;
                continue;
            }
            if ((a[i].nn[NN::TOP]->val && !a[i].nn[NN::LEFT]->val) || (!a[i].nn[NN::TOP]->val && a[i].nn[NN::LEFT]->val)) {
                if (a[i].nn[NN::TOP]->val) {
                    a[i].val = a[i].nn[NN::TOP]->val;
                    clusters[a[i].val]++;
                }
                else {
                    a[i].val = a[i].nn[NN::LEFT]->val;
                    clusters[a[i].val]++;
                }
                continue;
            }
            
            if (a[i].nn[NN::TOP]->val && a[i].nn[NN::LEFT]->val) {
                if (a[i].nn[NN::TOP]->val == a[i].nn[NN::LEFT]->val) {
                    a[i].val = a[i].nn[NN::LEFT]->val;
                    clusters[a[i].val]++;
                }
                else {
                    int sm = 0;
                    int lg = 0;
                    if (a[i].nn[NN::TOP]->val > a[i].nn[NN::LEFT]->val) {
                        sm = a[i].nn[NN::LEFT]->val;
                        lg = a[i].nn[NN::TOP]->val;
                    }
                    else {
                        sm = a[i].nn[NN::TOP]->val;
                        lg = a[i].nn[NN::LEFT]->val;
                    }
                    a[i].val = lg;
                    clusters[lg] += clusters[sm] + 1;
                    clusters.erase(sm);
                    unionClusters(a, L, sm, lg);
                    
                }
            }
            if (a[i].nn[NN::TOPLEFT]->val) {
                a[i].val = a[i].nn[NN::TOPLEFT]->val;
                clusters[a[i].val]++;
            }
        }
    }
    return clusters;
}


bool checkPath(Node* a,int L) {
    //init burn side
    int LL = L * L;
    for (int i = 0; i < L; i++) {
        if (a[i].val) {
            a[i].val = 2;
        }
    }
    //burn process
    bool burnFlag = true;
    bool endFlag = false;
    int burnLevel = 2;
    while (burnFlag) {
        burnFlag = false;
        for (int i = 0; i < LL; i++) {
            if (a[i].val == burnLevel) {
                burnFlag = true;
                if (a[i].nn[NN::TOPLEFT]->val == 1) {
                    a[i].nn[NN::TOPLEFT]->val = burnLevel + 1;
                }
                if (a[i].nn[NN::TOP]->val==1) {
                    a[i].nn[NN::TOP]->val = burnLevel + 1;
                }
                if (a[i].nn[NN::LEFT]->val==1) {
                    a[i].nn[NN::LEFT]->val = burnLevel + 1;
                }
                if (a[i].nn[NN::RIGHT]->val==1) {
                    a[i].nn[NN::RIGHT]->val = burnLevel + 1;
                }
                if (a[i].nn[NN::BOTTOM]->val==1) {
                    a[i].nn[NN::BOTTOM]->val = burnLevel + 1;
                    if (i >= LL - (2 * L)) {
                        goto stop;
                    }
                }
                if (a[i].nn[NN::BOTTOMRIGHT]->val == 1) {
                    a[i].nn[NN::BOTTOMRIGHT]->val = burnLevel + 1;
                    if (i >= LL - (2 * L)) {
                        goto stop;
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
    const int LL = L * L;
    const int T = params.T;
    const float p0 = params.p0;
    const float pk = params.pk;
    const float dp = params.dp;

    //lattice
    Node* a = new Node [LL];
    Node* b = new Node [LL];
    Node boundryNode = Node();
    boundryNode.val = 0;
    //nn init
    for (int i = 0; i < LL; i++) {
        //n0 i-L-1
        if (i>= L&& i % L != 0) {
            a[i].nn[NN::TOPLEFT] = &a[i - L-1];
            b[i].nn[NN::TOPLEFT] = &b[i - L-1];
        }
        else {
            a[i].nn[NN::TOPLEFT] = &boundryNode;
            b[i].nn[NN::TOPLEFT] = &boundryNode;
        }
        //n1 i-L
        if (i >= L) {
            a[i].nn[NN::TOP] = &a[i - L];
            b[i].nn[NN::TOP] = &b[i - L];
        }
        else {
            a[i].nn[NN::TOP] = &boundryNode;
            b[i].nn[NN::TOP] = &boundryNode;
        }
        //n2 i-1
        if (i % L != 0) {
            a[i].nn[NN::LEFT] = &a[i - 1];
            b[i].nn[NN::LEFT] = &b[i - 1];
        }
        else {
            a[i].nn[NN::LEFT] = &boundryNode;
            b[i].nn[NN::LEFT] = &boundryNode;
        }
        //n3 i+1
        if (i % L != L-1) {
            a[i].nn[NN::RIGHT] = &a[i + 1];
            b[i].nn[NN::RIGHT] = &b[i + 1];
        }
        else {
            a[i].nn[NN::RIGHT] = &boundryNode;
            b[i].nn[NN::RIGHT] = &boundryNode;
        }
        //n4 i+L
        if (i < LL-L) {
            a[i].nn[NN::BOTTOM] = &a[i + L];
            b[i].nn[NN::BOTTOM] = &b[i + L];
        }
        else {
            a[i].nn[NN::BOTTOM] = &boundryNode;
            b[i].nn[NN::BOTTOM] = &boundryNode;
        }
        //n5 i+L+1
        if (i < LL - L&& i % L != L - 1) {
            a[i].nn[NN::BOTTOMRIGHT] = &a[i + L+1];
            b[i].nn[NN::BOTTOMRIGHT] = &b[i + L+1];
        }
        else {
            a[i].nn[NN::BOTTOMRIGHT] = &boundryNode;
            b[i].nn[NN::BOTTOMRIGHT] = &boundryNode;
        }
    }

    //random
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    /*string fileName = "Ave_L_" + to_string(L) + "_T_" + to_string(T) + "_tmstmp_" + to_string(time(nullptr)) + ".txt";
    fstream mainOutput;
    mainOutput.open(fileName,ios::out);
    if (!mainOutput.is_open()) {
        cout << "file error" << endl;
        return -1;
    }*/

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
            for (int i = 0; i < LL; i++) {
                float rand = distribution(generator);
                if (rand < p) {
                    a[i].val = 1;
                    b[i].val = 1;
                }
                else {
                    a[i].val = 0;
                    b[i].val = 0;
                }
            }
            printConfig(b, L);
            /*auto clusters = findClusters(a, L);
            int max = 0;
            for (auto const& x : clusters) {
                clustersDistibution[x.second] += 1;
                if (x.second > max) {
                    max = x.second;
                }
            }
            s += max;*/
            bool path = checkPath(b, L);
            if (path) {
                pflow += 1;
            }
            printConfig(b, L);

        }
        for (auto const& x : clustersDistibution) {
            clustersDistibution[x.first] /= (float)T;
        }
        pflow = pflow / T;
        s = s / T;
        /*string fileName = "Dist_p_" + to_string(p) + "_L_" + to_string(L)+"_T_"+to_string(T)+ "_tmstmp_" + to_string(time(nullptr)) + ".txt";
        fstream clustersOutput;
        clustersOutput.open(fileName, ios::out);
        for (auto const& x : clustersDistibution) {
            cout << x.first << "\t" << x.second << endl;
            clustersOutput << x.first << "\t" << x.second << endl;
        }
        clustersOutput.close();*/
        cout << p << "\t" << pflow << "\t" << s << endl;
        //mainOutput << p << "\t"<<pflow<<"\t"<< s << endl;
        p += dp;
    }
    //mainOutput.close();

    delete[] a;
    delete[] b;
}





