#include <queue>
#include <mutex>
#include <atomic>

#define INT_MAX 2147483647

using namespace std;

struct partialMatching {
    int newMatch;
    int newMatchCost;
    int lb;
    bool operator<(const partialMatching& a) const
    {
        return lb > a.lb; // smaller at top
    }
    partialMatching() {
        newMatch = -1;
        newMatchCost = 0;
        lb = -1;
    }
    partialMatching(int _newMatch, int _newMatchCost, int _lb) {
        newMatch = _newMatch;
        newMatchCost = _newMatchCost;
        lb = _lb;
    }
};

struct node
{
    // only the top one is unstabled
    int* stabledMatch;
    int* stabledMatch2;
    // cost of stabled match
    int stabledCost;
    // length of stabled + unstabled match
    int lf; 
    priority_queue<partialMatching> unusedData;
    partialMatching data;
    node* left;
    node* right;
    node() {
        left = nullptr;
        right = nullptr;
        stabledMatch = nullptr;
        stabledMatch2 = nullptr;
        stabledCost = 0;
        lf = 0;
        //mark = false;
    }
    node(const node& x) {
        data = x.data;
        unusedData = x.unusedData;
        stabledMatch = x.stabledMatch;
        stabledMatch2 = x.stabledMatch2;
        stabledCost = x.stabledCost;
        left = nullptr;
        right = nullptr;
    }
    ~node() {
        delete[] stabledMatch;
        delete[] stabledMatch2;
    }
    node &operator=(const node& x) {
        data = x.data;
        unusedData = x.unusedData;
        stabledMatch = x.stabledMatch;
        stabledMatch2 = x.stabledMatch2;
        stabledCost = x.stabledCost;
        left = nullptr;
        right = nullptr;
        return *this;
    }
    bool operator>(const node& a) const
    {
        return data.lb > a.data.lb;
    }
};

struct dataBlock {
    int* f;
    int l;
    int lb;
    // -1 means unmatched, -2 means matched by other branches
    int* C;

    dataBlock() {
        f = nullptr;
        l = 0;
        lb = 0;
        C = nullptr;
    }

    dataBlock(int* _f, int _l, int _lb, int* _C) {
        f = _f;
        l = _l;
        lb = _lb;
        C = _C;
    }

    dataBlock& operator=(const dataBlock& x) {
        f = x.f;
        l = x.l;
        lb = x.lb;
        C = x.C;
        return *this;
    }

    bool operator> (const dataBlock& a) const {
        return lb < a.lb;
    }
    bool operator< (const dataBlock& a) const {
        return lb > a.lb;
    } 
    bool operator== (const dataBlock& a) const {
        return f == a.f && C == a.C;
    }
};

struct nodeHGED {
    dataBlock data;
    int children;
    
    bool isLeft;
    nodeHGED* parent;
    nodeHGED* left;
    nodeHGED* right;
    
    nodeHGED() {
        parent = nullptr;
        left = nullptr;
        right = nullptr;
        children = 0;
    }

    nodeHGED(dataBlock x) {
        data = x;
        parent = nullptr;
        left = nullptr;
        right = nullptr;
        children = 0;
    }
};

class minHeap{
public:
    node* top;
    int size;
    
    minHeap();
    void check();
    int checkHelper(node* x);
    void clear(node* x, int ub);
    void push(node* newNode, int ub);
    node* pop(int ub);
};

class heapHGED{
public:
    nodeHGED* top;

    heapHGED();
    void check(nodeHGED* x);
    void pushHelper(nodeHGED* x, nodeHGED* root);
    void push(dataBlock data);
    nodeHGED* popHelper(nodeHGED* x);
    dataBlock pop();
    void replaceHelper2(nodeHGED* x);
    nodeHGED* replaceHelper1(nodeHGED* x, dataBlock data);
    dataBlock replace(dataBlock data);
};

struct jobBlock {
    int wsize;
    int requester;
    heapHGED* requesterGlobal;
    priority_queue<dataBlock>* requesterQ;

    jobBlock() {
        wsize = INT_MAX;
        requester = -1;
        requesterGlobal = nullptr;
        requesterQ = nullptr;
    }
};

class GED{
public:
	priority_queue<dataBlock>* serialQ;
    	int ub;
	bool done;

    	int qSize;
    	int gSize;
    	int* V1, * V2, * E1, * E2;
    	int maxLabelEdge;
    	int maxLabelVertex;
	
	GED(int l1, int l2, int* _V1, int* _V2, int* _E1, int* _E2, int _maxE, int _maxV);
	int lb_Sa_StartOver(int* f, int lf, int* match2, int newMatch2);
	int start();
	void genNext(int* f, int l, int* C, priority_queue<dataBlock>* Q);
};

class GED2 {
public:
    minHeap* heapList;
    bool* avaliable;
    mutex lock;
    mutex waitLock;
    int nHeap;
    int totalThread;
    mutex lockUb;
    int upperBound;
    int waiting;
    int** _lvq;
    int** _lvg;

    int l1;
    int l2;
    int* V1;
    int* V2;
    int* E1;
    int* E2;

    atomic<long> totalWaiting;

    int maxLabelEdge;
    int maxLabelVertex;

    GED2(int _nHeap, int nThread, int _l1, int _l2, int* _V1, int* _V2, int* _E1, int* _E2, int _maxEdge, int _maxVertex);
    node* findBest(int waitState);
    void insert(node* x);
    int lb_Sa_StartOver(int id, int* f, int lf, int* match2, int newMatch2);
    void aStar(int id);
    void aStar2(int id);
};

class HGED {
public:
    mutex lock;
    int ub;
    atomic<int> done;
    mutex lockSharedArray;
    jobBlock* sharedArray;
    priority_queue<dataBlock>* serialQ;
    bool* waiting;
    int nThread;
    int beta;

    int qSize;
    int gSize;
    int* V1, * V2, * E1, * E2;
    int maxLabelEdge;
    int maxLabelVertex;

    HGED(int l1, int l2, int* _V1, int* _V2, int* _E1, int* _E2, int _maxE, int _maxV);
    HGED(int _nThread, int _beta, int l1, int l2, int* _V1, int* _V2, int* _E1, int* _E2, int _maxE, int _maxV);
    int lb_Sa_StartOver(int* f, int lf, int* match2, int newMatch2);
    void push(dataBlock m, priority_queue<dataBlock>* heap, heapHGED* globalHeap);
    dataBlock pop(int lv, priority_queue<dataBlock>* heap, heapHGED* globalHeap);
    void TwoPhaseSearch(int id, dataBlock* W);
    dataBlock** Start();
};
