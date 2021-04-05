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
    _X,
    X,
    _Y,
    Y,
    _Z,
    Z
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
    int L3 = L * L * L;
    for (int i = 0; i < L3; i++) {
        if (a[i].val == sm) {
            a[i].val = lg;
        }
    }
}
void unionClusters(Node* a, int L, int sm1,int sm2, int lg) {
    int L3 = L * L*L;
    for (int i = 0; i < L3; i++) {
        if (a[i].val == sm1||a[i].val==sm2) {
            a[i].val = lg;
        }
    }
}
map<int, int> findClusters(Node* a, int L) {
    int k = 2;
    int L3 = L * L * L;
    map<int, int> clusters;
    for (int i = 0; i < L3; i++) {
        if (a[i].val) {
            if (!a[i].nn[NN::_X]->val && !a[i].nn[NN::_Y]->val && !a[i].nn[NN::_Z]->val) {
                k++;
                a[i].val = k;
                clusters[k] = 1;
                continue;
            }
            if ((a[i].nn[NN::_X]->val && !a[i].nn[NN::_Y]->val && !a[i].nn[NN::_Z]->val)
                || (!a[i].nn[NN::_X]->val && a[i].nn[NN::_Y]->val && !a[i].nn[NN::_Z]->val)
                || (!a[i].nn[NN::_X]->val && !a[i].nn[NN::_Y]->val && a[i].nn[NN::_Z]->val)) {
                if (a[i].nn[NN::_X]->val) {
                    a[i].val = a[i].nn[NN::_X]->val;
                    clusters[a[i].val]++;
                }
                else if (a[i].nn[NN::_Y]->val) {
                    a[i].val = a[i].nn[NN::_Y]->val;
                    clusters[a[i].val]++;
                }
                else {
                    a[i].val = a[i].nn[NN::_Z]->val;
                    clusters[a[i].val]++;
                }
                continue;
            }
            if ((a[i].nn[NN::_X]->val && a[i].nn[NN::_Y]->val && !a[i].nn[NN::_Z]->val)
                || (a[i].nn[NN::_X]->val && !a[i].nn[NN::_Y]->val && a[i].nn[NN::_Z]->val)
                || (!a[i].nn[NN::_X]->val && a[i].nn[NN::_Y]->val && a[i].nn[NN::_Z]->val))
            {
                if (a[i].nn[NN::_X]->val && a[i].nn[NN::_Y]->val) {
                    if (a[i].nn[NN::_X]->val == a[i].nn[NN::_Y]->val) {
                        a[i].val = a[i].nn[NN::_X]->val;
                        clusters[a[i].val]++;
                    }
                    else {
                        int sm = 0;
                        int lg = 0;
                        if (a[i].nn[NN::_X]->val > a[i].nn[NN::_Y]->val) {
                            lg = a[i].nn[NN::_X]->val;
                            sm = a[i].nn[NN::_Y]->val;
                        }
                        else {
                            sm = a[i].nn[NN::_X]->val;
                            lg = a[i].nn[NN::_Y]->val;
                        }
                        a[i].val = lg;
                        clusters[lg] += clusters[sm] + 1;
                        clusters.erase(sm);
                        unionClusters(a, L, sm, lg);
                    }
                    continue;
                }
                if (a[i].nn[NN::_X]->val && a[i].nn[NN::_Z]->val) {
                    if (a[i].nn[NN::_X]->val == a[i].nn[NN::_Z]->val) {
                        a[i].val = a[i].nn[NN::_X]->val;
                        clusters[a[i].val]++;
                    }
                    else {
                        int sm = 0;
                        int lg = 0;
                        if (a[i].nn[NN::_X]->val > a[i].nn[NN::_Z]->val) {
                            lg = a[i].nn[NN::_X]->val;
                            sm = a[i].nn[NN::_Z]->val;
                        }
                        else {
                            sm = a[i].nn[NN::_X]->val;
                            lg = a[i].nn[NN::_Z]->val;
                        }
                        a[i].val = lg;
                        clusters[lg] += clusters[sm] + 1;
                        clusters.erase(sm);
                        unionClusters(a, L, sm, lg);
                    }
                    continue;
                }
                if (a[i].nn[NN::_Z]->val && a[i].nn[NN::_Y]->val) {
                    if (a[i].nn[NN::_Z]->val == a[i].nn[NN::_Y]->val) {
                        a[i].val = a[i].nn[NN::_Z]->val;
                        clusters[a[i].val]++;
                    }
                    else {
                        int sm = 0;
                        int lg = 0;
                        if (a[i].nn[NN::_Z]->val > a[i].nn[NN::_Y]->val) {
                            lg = a[i].nn[NN::_Z]->val;
                            sm = a[i].nn[NN::_Y]->val;
                        }
                        else {
                            sm = a[i].nn[NN::_Z]->val;
                            lg = a[i].nn[NN::_Y]->val;
                        }
                        a[i].val = lg;
                        clusters[lg] += clusters[sm] + 1;
                        clusters.erase(sm);
                        unionClusters(a, L, sm, lg);
                    }
                }
            }

            if (a[i].nn[NN::_X]->val && a[i].nn[NN::_Y]->val && a[i].nn[NN::_Z]->val) {
                if (a[i].nn[NN::_X]->val == a[i].nn[NN::_Y]->val && a[i].nn[NN::_X]->val == a[i].nn[NN::_Z]->val && a[i].nn[NN::_Y]->val == a[i].nn[NN::_Z]->val) {
                    a[i].val = a[i].nn[NN::_X]->val;
                    clusters[a[i].val]++;
                }
                else if (a[i].nn[NN::_X]->val == a[i].nn[NN::_Y]->val) {
                    int sm = 0;
                    int lg = 0;
                    if (a[i].nn[NN::_X]->val > a[i].nn[NN::_Z]->val) {
                        lg = a[i].nn[NN::_X]->val;
                        sm = a[i].nn[NN::_Z]->val;
                    }
                    else {
                        sm = a[i].nn[NN::_X]->val;
                        lg = a[i].nn[NN::_Z]->val;
                    }
                    a[i].val = lg;
                    clusters[lg] += clusters[sm] + 1;
                    clusters.erase(sm);
                    unionClusters(a, L, sm, lg);

                }
                else if (a[i].nn[NN::_X]->val == a[i].nn[NN::_Z]->val) {
                    int sm = 0;
                    int lg = 0;
                    if (a[i].nn[NN::_X]->val > a[i].nn[NN::_Y]->val) {
                        lg = a[i].nn[NN::_X]->val;
                        sm = a[i].nn[NN::_Y]->val;
                    }
                    else {
                        sm = a[i].nn[NN::_X]->val;
                        lg = a[i].nn[NN::_Y]->val;
                    }
                    a[i].val = lg;
                    clusters[lg] += clusters[sm] + 1;
                    clusters.erase(sm);
                    unionClusters(a, L, sm, lg);
                }
                else if (a[i].nn[NN::_Y]->val == a[i].nn[NN::_Z]->val)
                {
                    int sm = 0;
                    int lg = 0;
                    if (a[i].nn[NN::_Y]->val > a[i].nn[NN::_X]->val) {
                        lg = a[i].nn[NN::_Y]->val;
                        sm = a[i].nn[NN::_X]->val;
                    }
                    else {
                        sm = a[i].nn[NN::_Y]->val;
                        lg = a[i].nn[NN::_X]->val;
                    }
                    a[i].val = lg;
                    clusters[lg] += clusters[sm] + 1;
                    clusters.erase(sm);
                    unionClusters(a, L, sm, lg);
                }
                else {
                    int sm1 = 0;
                    int sm2 = 0;
                    int lg = 0;
                    if (a[i].nn[NN::_X]->val > a[i].nn[NN::_Y]->val) {
                        if (a[i].nn[NN::_X]->val > a[i].nn[NN::_Z]->val) {
                            sm1 = a[i].nn[NN::_Z]->val;
                            sm2 = a[i].nn[NN::_Y]->val;
                            lg = a[i].nn[NN::_X]->val;
                        }
                        else {
                            sm1 = a[i].nn[NN::_X]->val;
                            sm2 = a[i].nn[NN::_Y]->val;
                            lg = a[i].nn[NN::_Z]->val;
                        }

                    }
                    else {
                        if (a[i].nn[NN::_Y]->val > a[i].nn[NN::_Z]->val) {
                            sm1 = a[i].nn[NN::_Z]->val;
                            sm2 = a[i].nn[NN::_X]->val;
                            lg = a[i].nn[NN::_Y]->val;
                        }
                        else {
                            sm1 = a[i].nn[NN::_X]->val;
                            sm2 = a[i].nn[NN::_Y]->val;
                            lg = a[i].nn[NN::_Z]->val;
                        }
                    }
                    a[i].val = lg;
                    clusters[lg] += clusters[sm1] + clusters[sm2] + 1;
                    clusters.erase(sm1);
                    clusters.erase(sm2);
                    unionClusters(a, L, sm1, sm2, lg);

                }
            }

        }
    }
    return clusters;
}


bool checkPath(Node* a,int L) {
    //init burn side
    int L3 = L * L*L;
    int L2 = L * L;
    for (int i = 0; i < L2; i++) {
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
        for (int i = 0; i < L3; i++) {
            if (a[i].val == burnLevel) {
                burnFlag = true;
                if (a[i].nn[NN::_X]->val == 1) {
                    a[i].nn[NN::_X]->val = burnLevel + 1;
                }
                if (a[i].nn[NN::X]->val==1) {
                    a[i].nn[NN::X]->val = burnLevel + 1;
                }
                if (a[i].nn[NN::_Y]->val==1) {
                    a[i].nn[NN::_Y]->val = burnLevel + 1;
                }
                if (a[i].nn[NN::Y]->val==1) {
                    a[i].nn[NN::Y]->val = burnLevel + 1;
                    if (i >= L3 - (2*L2)) {
                        goto stop;
                    }
                }
                if (a[i].nn[NN::_Z]->val==1) {
                    a[i].nn[NN::_Z]->val = burnLevel + 1;
                }
                if (a[i].nn[NN::Z]->val == 1) {
                    a[i].nn[NN::Z]->val = burnLevel + 1;
                    
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
    const int L2 = L*L;
    const int L3 = L * L*L;
    const int T = params.T;
    const float p0 = params.p0;
    const float pk = params.pk;
    const float dp = params.dp;

    //lattice
    Node* a = new Node [L3];
    Node* b = new Node [L3];
    Node boundryNode = Node();
    boundryNode.val = 0;
    //nn init
    for (int i = 0; i < L3; i++) {
        //n0 _X i-1
        if ( i % L != 0) {
            a[i].nn[NN::_X] = &a[i - 1];
            b[i].nn[NN::_X] = &b[i - 1];
        }
        else {
            a[i].nn[NN::_X] = &boundryNode;
            b[i].nn[NN::_X] = &boundryNode;
        }
        //n1 X i+1
        if (i % L != L - 1) {
            a[i].nn[NN::X] = &a[i + 1];
            b[i].nn[NN::X] = &b[i + 1];
        }
        else {
            a[i].nn[NN::X] = &boundryNode;
            b[i].nn[NN::X] = &boundryNode;
        }
        //n2 _Y i-2*L
        if (i>=L2) {
            a[i].nn[NN::_Y] = &a[i - L2];
            b[i].nn[NN::_Y] = &b[i - L2];
        }
        else {
            a[i].nn[NN::_Y] = &boundryNode;
            b[i].nn[NN::_Y] = &boundryNode;
        }
        //n3 Y i+2*L
        if (i<L3-L2) {
            a[i].nn[NN::Y] = &a[i + L2];
            b[i].nn[NN::Y] = &b[i + L2];
        }
        else {
            a[i].nn[NN::Y] = &boundryNode;
            b[i].nn[NN::Y] = &boundryNode;
        }
        //n4 _Z i-L
        if (i >= L) {
            a[i].nn[NN::_Z] = &a[i - L];
            b[i].nn[NN::_Z] = &b[i - L];
        }
        else {
            a[i].nn[NN::_Z] = &boundryNode;
            b[i].nn[NN::_Z] = &boundryNode;
        }
        //n5 Z i+L
        if (i<L3-L) {
            a[i].nn[NN::Z] = &a[i + L];
            b[i].nn[NN::Z] = &b[i + L];
        }
        else {
            a[i].nn[NN::Z] = &boundryNode;
            b[i].nn[NN::Z] = &boundryNode;
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
            for (int i = 0; i < L3; i++) {
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
            //printConfig(a, L);
            auto clusters = findClusters(a, L);
            int max = 0;
            for (auto const& x : clusters) {
                clustersDistibution[x.second] += 1;
                if (x.second > max) {
                    max = x.second;
                }
            }
            /*s += max;
            bool path = checkPath(b, L);
            if (path) {
                pflow += 1;
            }*/
            //printConfig(a, L);

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
            //cout << x.first << "\t" << x.second << endl;
            clustersOutput << x.first << "\t" << x.second << endl;
        }
        clustersOutput.close();
        //cout << p << "\t" << pflow << "\t" << s << endl;
        //mainOutput << p << "\t"<<pflow<<"\t"<< s << endl;
        p += dp;
    }
    //mainOutput.close();

    delete[] a;
    delete[] b;
}





