#include <iostream>
#include <queue>
#include <memory>
#include <mutex>
#include <thread>
#include <atomic>
#include <string>
#include <vector>
#include <limits>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <cstring>
#include <fstream>
#include <chrono>
#include <tbb/concurrent_priority_queue.h>
#include <tbb/tick_count.h>
#include <ctime>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>

#define preAllocSize 60000000

using namespace std;
using namespace boost;

typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double>> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef property_map<Graph, edge_weight_t>::type WeightMap;
    

void set_affinity(int core_id) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    pthread_t current_thread = pthread_self();
    pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}

struct partialMatching {
    int newMatch;
    int newMatchCost;
    double lb;
    bool operator<(const partialMatching& a) const
    {
        return lb > a.lb; // smaller at top
    }
    bool operator>(const partialMatching& a) const
    {
        return lb < a.lb; // smaller at top
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
    int addr;
    int* stabledMatch;
    int* stabledMatch2;
    // cost of stabled match
    double stabledCost;
    // length of stabled + unstabled match
    int lf; 
    priority_queue<partialMatching> unusedData;
    partialMatching data;
    node() {
        stabledMatch = nullptr;
        stabledMatch2 = nullptr;
        stabledCost = 0;
        lf = 0;
        addr = -1;
        //mark = false;
    }
    node(const node& x) {
        addr = x.addr;
        data = x.data;
        lf = x.lf;
        unusedData = x.unusedData;
        stabledMatch = x.stabledMatch;
        stabledMatch2 = x.stabledMatch2;
        stabledCost = x.stabledCost;
    }
    node &operator=(const node& x) {
        data = x.data;
        unusedData = x.unusedData;
        stabledMatch = x.stabledMatch;
        stabledMatch2 = x.stabledMatch2;
        stabledCost = x.stabledCost;
        addr = x.addr;
        lf = x.lf;
        return *this;
    }node(node&& x) noexcept {
	addr = x.addr;
	data = std::move(x.data);
	lf = x.lf;
	unusedData = std::move(x.unusedData);
        stabledMatch = x.stabledMatch;
        stabledMatch2 = x.stabledMatch2;
	stabledCost = x.stabledCost;
	x.stabledMatch = nullptr;
	x.stabledMatch2 = nullptr;
    }
    node& operator=(node&& x) noexcept {
	if (this != &x) {
	    addr = x.addr;
	    data = std::move(x.data);
	    lf = x.lf;
	    unusedData = std::move(x.unusedData);
	    stabledMatch = x.stabledMatch;
	    stabledMatch2 = x.stabledMatch2;
	    stabledCost = x.stabledCost;
            x.stabledMatch = nullptr;
	    x.stabledMatch2 = nullptr;
	}
	return *this;
    }
    // for min-heap
    bool operator<(const node& a) const
    {
	if(data.lb != a.data.lb)
            return data.lb > a.data.lb;
	else
	    return lf < a.lf;
    }
    bool operator>(const node& a) const
    {
	if(data.lb != a.data.lb)
            return data.lb < a.data.lb;
	else
	    return lf > a.lf;
    }
    ~node() {}
};

class MPGED {
public:
    priority_queue<node>* heapList[256];
    atomic<bool>* avaliable;
    mutex lock;
    int nHeap;
    int totalThread;
    mutex lockUb;
    double upperBound;
    atomic<int> waiting;
    int** _lvq;
    int** _lvg;

    node emptyNode;

    bool* memUsed;
    int* memG1;
    int* memG2;
    
    bool* memUsed_2;
    int* memG1_2;
    int* memG2_2;

    bool socket2Ready;

    int l1;
    int l2;
    int* V1;
    int* V2;
    int* E1;
    int* E2;

    int maxLabelEdge;
    int maxLabelVertex;

    int* countCAS;
    int* countHS;

    MPGED(int _nHeap, int nThread, int _l1, int _l2, int* _V1, int* _V2, int* _E1, int* _E2, int _maxEdge, int _maxVertex) {
    	memUsed = new bool[preAllocSize]();
    	for (int i = 0; i < preAllocSize; i++) { memUsed[i] = 0; }
        memG1 = new int[preAllocSize*_l1]();
        memG2 = new int[preAllocSize*_l2]();

        avaliable = new atomic<bool>[_nHeap];
        for (int i = 0; i < _nHeap; i++) {
            avaliable[i].store(true);
        }
        nHeap = _nHeap;
        upperBound = INT_MAX;
        totalThread = nThread;
        waiting.store(0);

        l1 = _l1;
        l2 = _l2;
        V1 = _V1;
        V2 = _V2;
        E1 = _E1;
        E2 = _E2;

        maxLabelEdge = _maxEdge;
        maxLabelVertex = _maxVertex;

        int tmp = max(maxLabelVertex, maxLabelEdge);

        _lvq = new int*[nThread];
        _lvg = new int*[nThread];
        for (int i = 0; i < nThread; i++) {
            _lvq[i] = new int[tmp]();
            _lvg[i] = new int[tmp]();
        }

        emptyNode.stabledMatch = nullptr;
        emptyNode.stabledMatch2 = nullptr;

	countCAS = new int[nThread]();
	countHS = new int[nThread]();
    }

    node findBest(int waitState, int id) {
	    int N = 2;
	    int lastSearchLimit = 0;
	    int searchLimit = (nHeap+N-1)/N;
	    while(true){
		double minValue = INT_MAX;
            	int minIndex = -1;
            	for (int i2 = id+lastSearchLimit; i2 < id+min(nHeap, searchLimit); i2++) {
		    int i = i2%nHeap;
            	    if (!avaliable[i] || heapList[i]->empty())continue;
		    if (heapList[i]->top().data.lb >= upperBound) continue;
		    countHS[id]++;
		    if (minIndex < 0) {
			bool tmpBool = true;
			if(avaliable[i].compare_exchange_strong(tmpBool, false)){
			    countCAS[id]++;
            	            minIndex = i;
            	            minValue = heapList[i]->top().data.lb;
			}else{
			    countCAS[id]++;
			}
            	        continue;
            	    }
		    else if (minValue > heapList[i]->top().data.lb) {
			bool tmpBool = true;
            	        if(avaliable[i].compare_exchange_strong(tmpBool, false)){
			    countCAS[id]++;
			    avaliable[minIndex].store(true);
			    minIndex = i;
			    minValue = heapList[i]->top().data.lb;
			}else{
			    countCAS[id]++;
			}
			continue;
		    }
            	}
            	if (minIndex != -1 && minValue >= upperBound) {
		    avaliable[minIndex].store(true);
            	    minIndex = -1;
            	}
            	if (minIndex >= 0) {
		    if(heapList[minIndex]->empty() || heapList[minIndex]->top().data.lb >= upperBound){
			avaliable[minIndex].store(true);
			continue;
		    }
            	    node ret = heapList[minIndex]->top();
		    heapList[minIndex]->pop();
            	    avaliable[minIndex].store(true);
		    return ret;
            	}
            	else if (searchLimit >= nHeap){
		    return emptyNode;
            	}
		lastSearchLimit = searchLimit;
		searchLimit *= 2;
	    }
    }

    void insert(node x, int id) {
	    int minIndex = -1;
	    if(id<64){
        	while (minIndex<0) {
        	    minIndex = -1;
		    int tmpSize = min(totalThread, 64);
        	    for (int i = id; i < id+tmpSize; i++) {
		    	if(i>=min(totalThread, 64))i-=tmpSize;
	    	    	bool tmpBool = true;
                    	if(!avaliable[i].compare_exchange_strong(tmpBool, false)) continue;
                    	minIndex = i;
            	    	break;
	            }
            	}
	    }
	    else if(id>=64 && id<128){
        	while (minIndex<0) {
		    int tmpSize = min(totalThread, 128)-64;
        	    minIndex = -1;
        	    for (int i = id; i < id+tmpSize; i++) {
		    	if(i>=min(totalThread, 128))i-=tmpSize;
	    	    	bool tmpBool = true;
                    	if(!avaliable[i].compare_exchange_strong(tmpBool, false)) continue;
                    	minIndex = i;
            	    	break;
	            }
            	}
	    }
	    else if(id>=128 && id<192){
        	while (minIndex<0) {
		    int tmpSize = min(totalThread, 192)-128;
        	    minIndex = -1;
        	    for (int i = id; i < id+tmpSize; i++) {
		    	if(i>=min(totalThread, 192))i-=tmpSize;
	    	    	bool tmpBool = true;
                    	if(!avaliable[i].compare_exchange_strong(tmpBool, false)) continue;
                    	minIndex = i;
            	    	break;
	            }
            	}
	    }
	    else if(id>=192 && id<256){
        	while (minIndex<0) {
		    int tmpSize = min(totalThread, 256)-192;
        	    minIndex = -1;
        	    for (int i = id; i < id+tmpSize; i++) {
		    	if(i>=min(totalThread, 256))i-=tmpSize;
	    	    	bool tmpBool = true;
                    	if(!avaliable[i].compare_exchange_strong(tmpBool, false)) continue;
                    	minIndex = i;
            	    	break;
	            }
            	}
	    }
            heapList[minIndex]->push(x);
            avaliable[minIndex].store(true);
    }

    int lb_Sa_StartOver(int id, int* f, int lf, int* match2) {
        int LSa = 0;
        int i, j, tmp;
        // if f is a complete match
        if (lf == l1) {
            LSa = l2 - lf;
            int minus = 0;
            for (i = 0; i < l2; i++) {
                // find unmatched node
                if (match2[i] >= 0)continue;
                for (j = 0; j < l2; j++) {
                    // egde exist between unmatched node and other nodes
                    if (E2[i * l2 + j] >= 0) {
                        LSa++;
                        // edge between 2 unmatched nodes will be calculated twice
                        // so minus a half at the end
                        if (match2[j] < 0)minus += 1;
                    }
                }
            }
            return LSa - minus / 2;
        }

        int* lvq = _lvq[id];
        int* lvg = _lvg[id];
        int sum1;
        int sum2;

        // First part, Lv
        sum1 = 0;
        sum2 = 0;
        for (i = 0; i < maxLabelVertex; i++) {
            lvq[i] = 0;
            lvg[i] = 0;
        }
        for (i = lf; i < l1; i++) {
            lvq[V1[i]]++;
            sum1++;
        }
        for (i = 0; i < l2; i++) {
            if (match2[i]>=0)continue;
            lvg[V2[i]]++;
            sum2++;
        }
        LSa += max(sum1, sum2);
        for (i = 0; i < maxLabelVertex; i++) {
            LSa -= min(lvq[i], lvg[i]);
        }
        // Second part, LEi, egdes between unmatched nodes
        sum1 = 0;
        sum2 = 0;
        for (i = 0; i < maxLabelEdge; i++) {
            lvq[i] = 0;
            lvg[i] = 0;
        }
        for (i = lf; i < l1; i++) {
            for (j = i + 1; j < l1; j++) {
                tmp = E1[i * l1 + j];
                if (tmp >= 0) {
                    lvq[tmp]++;
                    sum1++;
                }
            }
        }
        for (i = 0; i < l2; i++) {
            if (match2[i]>=0)continue;
            for (j = i + 1; j < l2; j++) {
                if (match2[j]>=0)continue;
                tmp = E2[i * l2 + j];
                if (tmp >= 0) {
                    lvg[tmp]++;
                    sum2++;
                }
            }
        }
        LSa += max(sum1, sum2);
        for (i = 0; i < maxLabelEdge; i++) {
            LSa -= min(lvq[i], lvg[i]);
        }

        // Third part, LEc, edges between matched and unmatched nodes
        for (i = 0; i < lf; i++) {
            sum1 = 0;
            sum2 = 0;
            for (j = 0; j < maxLabelEdge; j++) {
                lvq[j] = 0;
                lvg[j] = 0;
            }
            for (j = lf; j < l1; j++) {
                tmp = E1[i * l1 + j];
                if (tmp >= 0) {
                    lvq[tmp]++;
                    sum1++;
                }
            }
            for (j = 0; j < l2; j++) {
                if (match2[j] >= 0)continue;
                else if (j == f[i]) continue;
                
                tmp = E2[f[i] * l2 + j];
                if (tmp >= 0) {
                    sum2++;
                    lvg[tmp]++;
                }
            }
            LSa += max(sum1, sum2);
            for (j = 0; j < maxLabelEdge; j++) {
                LSa -= min(lvq[j], lvg[j]);
            }
        }
        return LSa;
    }

    void aStar(int id) {
        node current = emptyNode;
        node current2 = emptyNode;
	int memTop = preAllocSize/totalThread*id;
	int memTopLimit = preAllocSize/totalThread*(id+1);
	if(id==0){
	    memTop+=l2;
	}
        int waitState = 0;

        int lf = -1;
        int* match = nullptr;
        int* match2 = nullptr;
        int* _match = nullptr;
        int* _match2 = nullptr;
        double stabledCost;
        while (true)
        {
	    
            // find a new job
            if (current.stabledMatch == nullptr) {
		if(waitState==1){waiting.fetch_add(1);}
                current = findBest(waitState, id);
                waitState += 1;
                // no thread have job to do, end
                if (int(waiting) >= totalThread) {
                    return;
                }
                continue;
            }
            // now I have a job
            if (waitState>1){waiting.fetch_add(-1);}
            waitState = 0;

            // new node found
            lf = current.lf;
            int addr=-1;
            while (addr<0) {
                if (!memUsed[memTop]) { addr = memTop; memUsed[addr] = true; }
                memTop = (memTop + 1);
	   	if(memTop>=memTopLimit)memTop=preAllocSize/totalThread*id;
            }
            match = &memG1[addr*l1];
            match2 = &memG2[addr*l2];
            stabledCost = current.stabledCost + current.data.newMatchCost;

            for (int i = 0; i < l1; i++) {
                match[i] = current.stabledMatch[i];
            }
            for (int i = 0; i < l2; i++) {
                match2[i] = current.stabledMatch2[i];
            }
            // lf after the new matching at this round
            lf++;
	    memUsed[addr] = 0;
	    
            // clean current
            priority_queue<partialMatching>().swap(current.unusedData);
            current.stabledMatch = match;
            current.stabledMatch2 = match2;
            current.stabledCost = stabledCost;
            current.lf = lf;
            current.addr = addr;
	    
            // find an unmatched node "i" in V2 to match the node "lf-1" in V1
            partialMatching tmpMatch;
            for (int i = 0; i < l2; i++) {
                // already matched node
                if (match2[i] >= 0)continue;

                match[lf - 1] = i;
                match2[i] = lf - 1;
                double lb = lb_Sa_StartOver(id, match, lf, match2);

                // flush back
                match2[i] = -1;
                match[lf - 1] = -1;

                // update matched cost
                int delta = 0;
                if (V1[lf-1] != V2[i]) delta++;
                for (int j = 0; j < lf-1; j++) {
                    if (E1[j * l1 + lf-1] != E2[match[j] * l2 + i])delta++;
                }
                lb = lb + delta + stabledCost;
                if (lb >= upperBound) {
                    continue;
                }
                tmpMatch.lb = lb;
                tmpMatch.newMatch = i;
                tmpMatch.newMatchCost = delta;
                current.unusedData.push(tmpMatch);
            };
            // prune
            if (current.unusedData.empty()) {
                current = emptyNode;
		memUsed[addr] = 0;
                lf = -1;
                continue;
            }
            // at leaf
            double xxx = current.unusedData.top().lb;
            if (lf == l1) {
                // better ub
                if (xxx < upperBound) {
                    lockUb.lock();
                    if (xxx < upperBound) {
                        upperBound = xxx;
                    }
                    lockUb.unlock();
                }
                current = emptyNode;
		memUsed[addr] = 0;
                lf = -1;
            }
            // all new lb not better than ub
            else if (xxx >= upperBound){
                current = emptyNode;
		memUsed[addr] = 0;
                lf = -1;
            }
            else {
		current.data = current.unusedData.top();
                current.unusedData.pop();
		while(!current.unusedData.empty()){
		    if(current.unusedData.top().lb >= upperBound){current.unusedData.pop();continue;}
		    addr = -1;
                	while (addr<0) {
                    	    if (!memUsed[memTop]) { addr = memTop; memUsed[addr] = true; }
                    	    memTop = (memTop + 1);
			    if(memTop>=memTopLimit)memTop=preAllocSize/totalThread*id;
            		}
            		_match = &memG1[addr*l1];
            		_match2 = &memG2[addr*l2];
            	    for (int i = 0; i < l1; i++) {
                	_match[i] = match[i];
        	    }
        	    for (int i = 0; i < l2; i++) {
        	        _match2[i] = match2[i];
        	    }
		    current2.data = current.unusedData.top();
		    current.unusedData.pop();
		    current2.stabledCost = current.stabledCost;
        	    current2.stabledMatch = _match;
        	    current2.stabledMatch2 = _match2;
        	    current2.lf = current.lf;
        	    current2.addr = addr;
            	    _match[lf - 1] = current2.data.newMatch;
            	    _match2[current2.data.newMatch] = lf - 1;
		    insert(current2, id);
		}
                match[lf - 1] = current.data.newMatch;
		match2[current.data.newMatch] = lf - 1;
            }
        }
    }
};


void middle(int id, MPGED* multi) {
    set_affinity(id+1);
    multi->aStar(id);
}

void readOneGraph(int l, int* V, int* E, ifstream* data){
	char line[4096];
	data->getline(line, 4096);
	while(line[0] != 'v'){
		data->getline(line, 4096);
	}
	int a = 0;
	while(line[0] == 'v'){
		int i = 2;
		while(line[i++] != ' ');
		string x(&line[i]);
		V[a] = stoi(x) - 1;
		data->getline(line, 4096);
		a++;
	}
	while(line[0] == 'e'){
		int i = 1;
		while(line[i++] != ' ');
		int j = i;
		while(line[j++] != ' ');
		int k = j;
		while(line[k++] != ' ');
		line[j-1] = '\0';
		line[k-1] = '\0';

		string x(&line[i]);
		string y(&line[j]);
		string z(&line[k]);
		i = stoi(x);
		j = stoi(y);
		k = stoi(z) - 1;
		E[i*l+j] = k;
		E[j*l+i] = k;
                data->getline(line, 4096);
	}
}

int main(int argc, char* argv[]) {
    set_affinity(0);
    
    int l1, l2, maxLabelEdge, maxLabelVertex;
    int nThread = 4;
    ifstream* data = new ifstream();
	
    if(argc != 4){
    	printf("your input should be: <data file name> <# vertices> <# thread>\n");
	return 0;
    }else{
        data->open(argv[1]);
        string x(argv[2]);
        l1 = stoi(x);
        l2 = stoi(x);
	string y(argv[3]);
	nThread = stoi(y);
    }
    int* V1 = new int[l1]();
    int* E1 = new int[l1*l1]();
    for(int i = 0; i < l1; i++){
        for(int j = 0; j < l1; j++){
            E1[i*l1+j] = -1;
        }
    }
    readOneGraph(l1, V1, E1, data);

    int* V2 = new int[l2]();
    int* E2 = new int[l2*l2]();
    for(int i = 0; i < l2; i++){
        for(int j = 0; j < l2; j++){
            E2[i*l2+j] = -1;
        }
    }
    readOneGraph(l2, V2, E2, data);

    maxLabelEdge = 0;
    maxLabelVertex = 0;
    for(int i = 0; i < l1; i++){
        if(maxLabelVertex < V1[i]) maxLabelVertex = V1[i];
    }
    for(int i = 0; i < l1; i++){
        for(int j = 0; j < l1; j++){
            if(maxLabelEdge < E1[i*l1+j]) maxLabelEdge = E1[i*l1+j];
        }
    }

    for(int i = 0; i < l2; i++){
        if(maxLabelVertex < V2[i]) maxLabelVertex = V2[i];
    }
    for(int i = 0; i < l2; i++){
        for(int j = 0; j < l2; j++){
            if(maxLabelEdge < E2[i*l2+j]) maxLabelEdge = E2[i*l2+j];
        }
    }
    maxLabelVertex += 1;
    maxLabelEdge += 1;
    
	if(nThread==64 || nThread==128 || nThread==192 || nThread==256)nThread-=1;
        MPGED* multi = new MPGED(nThread, nThread, l1, l2, V1, V2, E1, E2, maxLabelEdge, maxLabelVertex);
    	for(int i = 0; i<min(64, nThread); i++){
		multi->heapList[i] = new priority_queue<node>();
	}
	if(nThread>64){
	        set_affinity(64);
    		for(int i = 64; i<min(nThread, 128); i++){
			multi->heapList[i] = new priority_queue<node>();
		}
	    	set_affinity(0);
	}
	if(nThread>128){
	        set_affinity(128);
    		for(int i = 128; i<min(nThread, 192); i++){
			multi->heapList[i] = new priority_queue<node>();
		}
	    	set_affinity(0);
	}
	if(nThread>192){
	        set_affinity(192);
    		for(int i = 192; i<min(nThread, 256); i++){
			multi->heapList[i] = new priority_queue<node>();
		}
	    	set_affinity(0);
	}
	chrono::time_point<std::chrono::steady_clock> startTime, endTime;
        std::chrono::duration<double> totalDuration;
	
	startTime = std::chrono::steady_clock::now();
        for (int i = 0; i < l2; i++) {
            node tmp;
            tmp.addr = i;
            multi->memUsed[i] = 1;
            tmp.stabledMatch = &multi->memG1[i * l1];
            for (int j = 0; j < l1; j++) {
                tmp.stabledMatch[j] = -1;
            }
            tmp.stabledMatch2 = &multi->memG2[i * l2];
            for (int j = 0; j < l2; j++) {
                tmp.stabledMatch2[j] = -1;
            }
            tmp.lf = 1;
            tmp.stabledCost = 0;

            partialMatching data(i, V1[0] != V2[i], V1[0] != V2[i]);
            tmp.data = data;
            tmp.stabledMatch[0] = tmp.data.newMatch;
            tmp.stabledMatch2[tmp.data.newMatch] = 0;
	    multi->heapList[i % nThread]->push(tmp);
        }
        
	thread* x;
        for (int i = 0; i < nThread; i++) {
            x = new thread(middle, i, multi);
            x->detach();
        }
	double last = 0;
        while (int(multi->waiting) < nThread) {
            endTime = std::chrono::steady_clock::now();
	    totalDuration = endTime - startTime;
	    if(totalDuration.count() >= 3600){printf("unfinished\n");return 1;}
            if(totalDuration.count() - last < 0.01)continue;
	    last = totalDuration.count();
	    printf("%f\t", totalDuration.count());
	    printf("%d\n", int(multi->upperBound));
        }
}
