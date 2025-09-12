#include "road_network.h"
#include "radix_heap.h"
#include "util.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <thread>
#include <sys/resource.h>
#include <queue>
#include <unordered_map>
#include <array>
#include <omp.h>
#include <math.h>
#include <random>
#include <cmath>
#include <stack>

#include <vector>
#include <future>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>

#define pb push_back
using namespace std;
using namespace road_network;

#define REMOVE_REDUNDANT
#define CONTRACT
#include "road_network.h"
const size_t repeats = 10;
const size_t nr_queries = 10000000;
const size_t nr_query_tests = 10;
const size_t nr_buckets = 10;
const size_t bucket_size = 10000;
const distance_t bucket_min = 1000;

const size_t MB = 1024 * 1024;
class ThreadPool {
public:
    explicit ThreadPool(size_t n)
    : stop_(false) {
        workers_.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            workers_.emplace_back([this] {
                for (;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lk(mtx_);
                        cv_.wait(lk, [this] { return stop_ || !tasks_.empty(); });
                        if (stop_ && tasks_.empty()) return;
                        task = std::move(tasks_.front());
                        tasks_.pop();
                    }
                    task();
                }
            });
        }
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    template <class F, class... Args>
    auto submit(F&& f, Args&&... args)
    -> std::future<typename std::invoke_result<F, Args...>::type> {
        using R = typename std::invoke_result<F, Args...>::type;
        auto task = std::make_shared<std::packaged_task<R()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        std::future<R> fut = task->get_future();
        {
            std::lock_guard<std::mutex> lk(mtx_);
            tasks_.emplace([task]{ (*task)(); });
        }
        cv_.notify_one();
        return fut;
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lk(mtx_);
            stop_ = true;
        }
        cv_.notify_all();
        for (auto& w : workers_) w.join();
    }

private:
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    std::mutex mtx_;
    std::condition_variable cv_;
    bool stop_;
};

// Reusable global pool (31 workers to match original code’s default)
static ThreadPool& global_pool() {
    static ThreadPool pool(31);
    return pool;
}
struct ResultData
{
    size_t label_count;
    size_t max_label_count;
    size_t index_size;
    size_t index_height;
    double index_time;
    double avg_cut_size;
    size_t max_cut_size;
    size_t pruning_2hop;
    size_t pruning_3hop;
    size_t pruning_tail;
    double random_query_time;
    double random_hoplinks;
    vector<double> bucket_query_times;
    vector<double> bucket_hoplinks;
};

struct FileResults
{
    string filename;
    vector<ResultData> results;
    FileResults(string filename, vector<ResultData> results) : filename(filename), results(results) {}
};

ostream& operator<<(ostream& os, const FileResults& fr)
{
    if (fr.results.empty())
        return os;
    os << endl << "Summary for " << fr.filename << ":" << endl << setprecision(5);
    os << "Index size (MB): " << util::summarize(fr.results, [](ResultData r) -> double { return r.index_size; }) << endl;
    os << "Index time (s): " << util::summarize(fr.results, [](ResultData r) -> double { return r.index_time; }) << endl;
    os << "Index height: " << util::summarize(fr.results, [](ResultData r) -> double { return r.index_height; }) << endl;
    os << "Avg cut size: " << util::summarize(fr.results, [](ResultData r) -> double { return r.avg_cut_size; }) << endl;
    os << "Max cut size: " << util::summarize(fr.results, [](ResultData r) -> double { return r.max_cut_size; }) << endl;
    os << "Query time (s): " << util::summarize(fr.results, [](ResultData r) -> double { return r.random_query_time; }) << endl;
    os << "Avg Hoplinks: " << util::summarize(fr.results, [](ResultData r) -> double { return r.random_hoplinks; }) << endl;
    if (!fr.results[0].bucket_query_times.empty())
        for (size_t bucket = 0; bucket < nr_buckets; bucket++)
        {
            os << "Bucket " << bucket << ": time = " << util::summarize(fr.results, [bucket](ResultData r) -> double { return r.bucket_query_times[bucket]; }) * (nr_queries / bucket_size);
            os << ", hoplinks = " << util::summarize(fr.results, [bucket](ResultData r) -> double { return r.bucket_hoplinks[bucket]; }) << endl;
        }
        return os;
}

#ifdef PRUNING
size_t get_2hop_pruning(const vector<CutIndex> &ci)
{
    size_t total = 0;
    for (NodeID node = 1; node < ci.size(); node++)
        total += ci[node].pruning_2hop;
    return total;
}

size_t get_3hop_pruning(const vector<CutIndex> &ci)
{
    size_t total = 0;
    for (NodeID node = 1; node < ci.size(); node++)
        total += ci[node].pruning_3hop;
    return total;
}

size_t get_tail_pruning(const vector<CutIndex> &ci)
{
    size_t total = 0;
    for (NodeID node = 1; node < ci.size(); node++)
        total += ci[node].pruning_tail;
    return total;
}
#endif
// helper struct to enque nodes by distance
struct SearchNode
{
    distance_t distance;
    NodeID node;
    // reversed for min-heap ordering
    bool operator<(const SearchNode& other) const { return distance > other.distance; }
    SearchNode(distance_t distance, NodeID node) : distance(distance), node(node) {}
};

/*
 v *oid run_dijkstra_llsub_par_upd(const vector<CutIndex> &ci,const std::vector<NodeID>& vertices)
 {
 CHECK_CONSISTENT;
 vector<thread> threads;
 auto dijkstra = [this](NodeID v, size_t distance_id) {
 assert(contains(v));
 assert(distance_id < MULTI_THREAD_DISTANCES);
 const uint16_t pruning_level = node_data[v].landmark_level;
 // init distances
 for (NodeID node : nodes)
     node_data[node].distances[distance_id] = infinity;
 node_data[v].distances[distance_id] = 0;
 // init queue
 priority_queue<SearchNode> q;
 q.push(SearchNode(0, v));
 // dijkstra
 while (!q.empty())
 {
 SearchNode next = q.top();
 q.pop();

 for (Neighbor n : node_data[next.node].neighbors)
 {
 Node& n_data = node_data[n.node];
 // filter neighbors nodes not belonging to subgraph or having higher landmark level
 if (!contains(n.node) || n_data.landmark_level >= pruning_level)
     continue;
 // update distance and enque
 distance_t new_dist = next.distance + n.distance;
 if (new_dist < n_data.distances[distance_id])
 {
 n_data.distances[distance_id] = new_dist;
 q.push(SearchNode(new_dist, n.node));
 }
 }
 }
 };
 for (size_t i = 0; i < vertices.size(); i++)
     threads.push_back(thread(dijkstra, vertices[i], i));
 for (size_t i = 0; i < vertices.size(); i++)
     threads[i].join();
 }
 */
//cut->a->b->border
inline static uint16_t get_offset(const uint16_t* dist_index, size_t cut_level)
{
    return cut_level ? dist_index[cut_level - 1] : 0;
}
void Graph::run_dij_upd_down(ContractionIndex& con_index,const vector<CutIndex>& ci, const vector<NodeID>& bor,NodeID a,NodeID b,distance_t w,NodeID cut, vector<pair<pair<NodeID, NodeID>, distance_t>> &change)
{

    CHECK_CONSISTENT;

    // init distances
    unordered_map<NodeID, distance_t>dist;
    priority_queue<SearchNode> q;return;
    for (NodeID x : bor) {
        dist[x] = con_index.get_distance(cut, a) + w + (b == x ? 0 : con_index.get_distance(b, x));
        q.push(SearchNode(dist[x], x));
    }
    distance_t mb = 0;

    //return;
    // init queue
    // dijkstra
    unordered_map<NodeID, bool>vis;

    ContractionLabel cv = con_index.labels[cut];
    FlatCutIndex cu = cv.cut_index;
    uint16_t cut_level = cu.cut_level();
    //if (cut_level ) { cout << "?"; exit(0); }
    //return;
    int tot = 0;
    uint16_t cof = node_data[cut].cut_offset;
    //cout << cof << endl;
    while (!q.empty())
    {
        SearchNode next = q.top();
        q.pop();
        if (vis[next.node])continue;
        vis[next.node] = true;
        if (con_index.is_contracted(next.node))continue;
        if (!con_index.isAncestor(cut, next.node))continue;
        if (con_index.get_distance(cut,next.node,cof) <= node_data[next.node].distance)continue;
        tot++;
        change.pb({ {cut,next.node},dist[next.node]});
        {
            /*
             *           //upd.pb({ {cut,next.node},next.distance });
             *           FlatCutIndex x = con_index.labels[next.node].cut_index;
             *           uint16_t x_offset = get_offset(x.dist_index(), cut_level);
             *           //if (cof + x_offset >= x.dist_index()[cut_level])
             *           // if (cut_level && x.dist_index()[cut_level] < x.dist_index()[cut_level - 1])continue;
             *           //cout << x.dist_index()[cut_level]<<" "<< x_offset << " " << cof << " " << x.dist_index()[cut_level] - x_offset << " !" << cut << " " << next.node << " " << x.distances() + x_offset + cof << " " << a << "  " << b << " " << cut_level << endl;
             *           if (cof + x_offset < x.dist_index()[cut_level] && cof < 10000) {
             *               //x.distances();
             *               distance_t* x_ptr = x.distances() + x_offset + cof;
             *x_ptr = next.distance;
        }
        */
        }
        //upd.pb({ {cut,next.node},node_data[next.node].distance});
        for (Neighbor n : node_data[next.node].neighbors)
        {
            // filter neighbors nodes not belonging to subgraph or having higher landmark level
            if (con_index.labels[n.node].cut_index.cut_level() >cut_level)
                continue;
            // update distance and enque
            distance_t new_dist = next.distance + n.distance;
            if (!dist.count(n.node)||new_dist < dist[n.node])
            {
                dist[n.node] = new_dist;
                q.push(SearchNode(new_dist, n.node));
            }
        }
    }
    cout << tot << endl;

}
//down == true if reduce w
void Graph::find_affected_point(const ContractionIndex &con_index,distance_t &get_dis, NodeID u, NodeID v,distance_t w,distance_t bfw, NodeID etime) {
    CHECK_CONSISTENT;
    assert(contains(u));
    assert(contains(v));
    // init distances
    get_dis=0;
    bool down = w < bfw;
    vector<NodeID>fgNode;
    node_data[u].time[32] = etime;
    node_data[v].time[32] = etime;
    node_data[u].distance = 0;
    // init queue
    node_data[v].distance = 0;
    priority_queue<SearchNode> q;
    q.push(SearchNode(0, u));
    q.push(SearchNode(0, v));
    // dijkstra
    int p = 0;
    //util::start_timer();
    int mk = nodes.size()<=300000?5000:10000;
    distance_t mxdis = 0;
    while (!q.empty())
    {
        SearchNode next = q.top();
        q.pop();
        get_dis = max(get_dis, next.distance);
        p++;
        if (p>=mk)break;
        for (Neighbor n : node_data[next.node].neighbors)
        {
            // filter neighbors nodes not belonging to subgraph
            // update distance and enque
            distance_t new_dist = next.distance + n.distance;
            if (node_data[n.node].time[32]!=etime||new_dist < node_data[n.node].distance)
            {
                node_data[n.node].time[32] = etime;
                node_data[n.node].distance = new_dist;
                q.push(SearchNode(new_dist, n.node));
            }
        }
    }
    //cout << util::stop_timer() << " "<<mxdis<<endl;
    return;
}


void Graph::update_down_par(size_t start, size_t end, const std::vector<NodeID>& nodes, const ContractionIndex& con_index, const vector<pair<pair<NodeID, NodeID>, pair<distance_t, distance_t>>>& edge)
{
    vector<NodeID> affected_cut;
    for (size_t i = start; i < end; ++i) {
        NodeID cut = nodes[i];
        distance_t max_cut_bd_distance = 0;
        int bdsz = node_data[cut].border.size();
        unordered_map<NodeID, bool>isbd;
        for (NodeID bd : node_data[cut].border) {
            max_cut_bd_distance = max(max_cut_bd_distance, con_index.get_distance(cut, bd));
            isbd[bd] = true;
        }
        priority_queue<SearchNode>q;
        unordered_map<NodeID, distance_t>dist;
        unordered_map<NodeID, bool>vis,vis2,vis3;
        q.push(SearchNode(0, cut));
        size_t totbd=0;
        int p = 0,t=0;
        while (!q.empty()) {
            NodeID x = q.top().node;
            distance_t ww = q.top().distance;
            q.pop();if (!con_index.isAncestor(cut, x))continue;
            //add cut ->v
            if (vis[x])continue;
            vis[x] = 1;
            if (isbd[x])totbd++;
            if (vis2.size() == edge.size())break;//used all edges
            if (ww >= max_cut_bd_distance)break;
            if (totbd == bdsz)break;
            for (auto e : node_data[x].neighbors) {
                distance_t disw;
                if (e.distance >= (distance_t)200000000) {
                    disw = edge[e.distance - (distance_t)200000000].second.second;
                    vis2[e.distance - (distance_t)200000000] = true;
                }
                else disw = e.distance;
                if (!dist.count(e.node)||dist[e.node]>ww+disw) {
                    dist[e.node] = ww + disw;
                    q.push({ ww+disw, e.node});
                }
            }
        }
        //cout << p <<"  "<<t<< endl;

        //cout << q2.size() << endl;
    }

}

vector<NodeID>edge_fg;
void Graph::change_remain(size_t start, size_t end,
                          std::vector<NodeID>& cutt,
                          ContractionIndex& con_index,
                          const std::vector<std::pair<std::pair<NodeID, NodeID>,
                          std::pair<distance_t, distance_t>>>& edge,
                          const std::vector<CutIndex>& ci,
                          NodeID tid)
{
    int tott = 0, tott2 = 0, tot3 = 0;

    // alias for our radix�\heap based PQ
    using RadixPQ = radix_heap::pair_radix_heap<distance_t, NodeID>;

    constexpr distance_t MAX_DIST = 200000000;

    for (size_t i = start; i < end; i++) {
        NodeID cut = cutt[i];

        RadixPQ q;   // ԭ���� q
        RadixPQ q2;  // ԭ���� q2
        RadixPQ q3;  // ԭ���� q3

        bool ff = false;
        std::vector<NodeID> affected_bd;
        ContractionLabel cu = con_index.labels[cut];
        uint16_t cut_level = cu.cut_index.cut_level();
        uint16_t cof = node_data[cut].cut_offset;
        distance_t mx = 0, mn = std::numeric_limits<distance_t>::infinity();
        NodeID time_now = static_cast<NodeID>(2 * (i + 1));

        for (size_t j = 0; j < edge.size(); j++) {
            auto e = edge[j];
            NodeID a = e.first.first, b = e.first.second;
            distance_t bfw = e.second.first, w = e.second.second;

            if (!con_index.isAncestor(edge_fg[j], cut)
                && !con_index.isAncestor(cut, edge_fg[j]))
            {
                continue;
            }

            distance_t disa = con_index.get_distance(cut, a);
            distance_t disb = con_index.get_distance(cut, b);
            mx = std::max(mx, std::min(disa, disb));
            if (disa > disb) ff = true;

            if (con_index.isAncestor(cut, a)
                && con_index.isAncestor(cut, b))
            {
                // ԭ���� q3.push(SearchNode(ff ? disb + bfw : disa + bfw, (ff ? a : b)));
                q3.push(ff ? (disb + bfw) : (disa + bfw),
                        ff ? a : b);
            }

            for (NodeID bd : node_data[cut].border) {
                if (node_data[bd].isbdtime[tid] == time_now)
                    continue;
                FlatCutIndex &x = con_index.labels[bd].cut_index;
                uint16_t x_offset = get_offset(x.dist_index(), cut_level);
                const distance_t* x_ptr = x.distances() + x_offset + cof;
                distance_t r = (ff ? disb : disa)
                + bfw
                + con_index.get_distance(bd, (ff ? a : b));
                if (*x_ptr + con_index.labels[bd].distance_offset == r) {
                    node_data[bd].isbdtime[tid] = time_now;
                    // ԭ���� q3.push(SearchNode(r, bd));
                    q3.push(r, bd);
                    affected_bd.push_back(bd);
                    mn = std::min(mn, *x_ptr + con_index.labels[bd].distance_offset);
                }
            }
        }

        if (q3.empty()) continue;

        std::vector<NodeID> wd;
        wd.reserve(affected_bd.size());
        //FlatCutIndex cnd;

        while (!q3.empty()) {
            // SearchNode next = q3.top(); q3.pop();
            distance_t dist_u = q3.top_key();
            NodeID     u = q3.top_value();
            q3.pop();

            auto &cnd = con_index.labels[u].cut_index;
            size_t off_u = get_offset(cnd.dist_index(), cut_level) + cof;
            if (cnd.distances()[off_u] != dist_u) {
                wd.push_back(u);
                continue;
            }

            node_data[u].viss[tid] |= 2;
            for (Neighbor n : node_data[u].neighbors) {
                auto &cnc = con_index.labels[u].cut_index;
                if (node_data[n.node].vistime[tid] != time_now) {
                    node_data[n.node].viss[tid] = 1;
                    node_data[n.node].vistime[tid] = time_now;
                }
                else {
                    node_data[n.node].viss[tid] |= 1;
                }
                if (cnc.cut_level() < cut_level) continue;

                distance_t disw =
                (n.distance >= MAX_DIST)
                ? edge[n.distance - MAX_DIST].second.first
                : n.distance;
                distance_t new_d = dist_u + disw;

                if (node_data[n.node].time[tid] != time_now
                    || node_data[n.node].distances[tid] > new_d) {
                    node_data[n.node].distances[tid] = new_d;
                node_data[n.node].time[tid] = time_now;
                // ԭ���� q3.push(SearchNode(new_d, n.node));
                q3.push(new_d, n.node);
                    }
            }
        }

        for (NodeID bd : wd) {
            if (node_data[bd].vistime[tid] == time_now
                && node_data[bd].viss[tid] == 1)
            {
                q2.push(node_data[bd].distances[tid], bd);
            }
        }

        node_data[cut].distances[tid] = 0;
        node_data[cut].time[tid] = time_now;
        // ԭ���� q.push(SearchNode(mn, cut));
        q.push(mn, cut);

        std::vector<SearchNode> s;
        // using DUAL to get the border dist.
        // con_index(fd,file_cg);
        // get_border_dist(fd,cut,affected_bd);
        /*
         *       // or A*
         *               while (!q.empty()) {
         *           SearchNode next = q.top();
         *           q.pop();
         *           if (node_data[next.node].vistime[tid] == time_now)
         *               continue;
         *           node_data[next.node].vistime[tid] = time_now;
         *           if (node_data[next.node].dist[tid] > mx) {
         *               for (NodeID bd : affected_bd) {
         *                   distance_t g = node_data[next.node].dist[tid] + con_index.get_distance(next.node, bd);
         *                   if (g < node_data[bd].dist[tid])
         *                       node_data[bd].dist[tid] = min(node_data[bd].dist[tid], g);
    }
    continue;
    }
    if (node_data[next.node].isbdtime[tid] == time_now) {
        s.push_back(SearchNode(node_data[next.node].dist[tid], next.node));
        auto it = std::find(affected_bd.begin(), affected_bd.end(), next.node);
        if (it != affected_bd.end())
            affected_bd.erase(it);
        if (affected_bd.empty())
            break;
        vector<NodeID> change;
        while (!q.empty()) {
            change.push_back(q.top().node);
            q.pop();
    }
    for (NodeID x : change) {
        distance_t best = infinity;
        for (NodeID bd : affected_bd)
            best = min(best, node_data[x].dist[tid] + con_index.get_distance(x, bd));
        if (best < infinity)
            q.push(SearchNode(best, x));
    }
    }
    for (Neighbor n : node_data[next.node].neighbors) {
        distance_t new_dist = node_data[next.node].dist[tid] + n.distance;
        if (node_data[n.node].disttime[tid] != time_now || new_dist < node_data[n.node].dist[tid]) {
            node_data[n.node].dist[tid] = new_dist;
            node_data[n.node].disttime[tid] = time_now;
            distance_t best = infinity;
            for (NodeID bd : affected_bd)
                best = min(best, new_dist + con_index.get_distance(n.node, bd));
        if (best < infinity)
            q.push(SearchNode(best, n.node));
    }
    }
    }
    */

        // ����׶Σ��� affected_bd ֱ��תΪ q2
        for (NodeID x : affected_bd) {
            s.push_back(SearchNode(node_data[x].dist[tid], x));
        }
        for (auto sn : s) {
            node_data[sn.node].distances[tid] = sn.distance;
            node_data[sn.node].disttime[tid] = time_now;
            // ԭ���� q2.push(sn);
            q2.push(sn.distance, sn.node);
        }

        // �����׶Σ����� q2����� while (!q2.empty())
        int local_tott = 0;
        while (!q2.empty()) {
            // ԭ����
            // NodeID node = q2.top().node;
            // distance_t ww = q2.top().distance;
            distance_t ww = q2.top_key();
            NodeID     node = q2.top_value();
            q2.pop();

            if (node_data[node].visstime[tid] == time_now)
                continue;
            node_data[node].visstime[tid] = time_now;

            for (Neighbor e : node_data[node].neighbors) {
                NodeID v = e.node;
                distance_t disw;
                if (node_data[v].time[tid] != time_now
                    || node_data[v].viss[tid] != 3)
                    continue;
                if (e.distance >= MAX_DIST)
                    disw = edge[e.distance - MAX_DIST].second.second;
                else
                    disw = e.distance;
                if (node_data[v].disttime[tid] != time_now
                    || node_data[v].distances[tid] > ww + disw) {
                    node_data[v].distances[tid] = ww + disw;
                node_data[v].disttime[tid] = time_now;
                // ԭ���� q2.push(SearchNode(node_data[v].distances[tid], v));
                q2.push(node_data[v].distances[tid], v);
                    }
            }
        }
        // cout << local_tott << " " << tot+affected_bd.size() << endl;
    }
    // cout << tott << " " << tott2 << " " << tot3 << endl;
}


void Graph::change_remain_one2(size_t start, size_t end, std::vector<NodeID>& cutt, ContractionIndex& con_index, distance_t w, distance_t bfw, NodeID a, NodeID b, const std::vector<CutIndex>& ci, vector<pair<pair<NodeID, NodeID>, distance_t>>& change, const std::unordered_map<NodeID, distance_t>& get_dis, distance_t md, NodeID tid, NodeID etime) {
    //int tot = 0;
    for (size_t t = start; t < end; t++) {
        NodeID node = cutt[t], u = a, v = b, cut = cutt[t];
        NodeID time = 2 * (t + 1) + etime * 147;

        distance_t disa = con_index.get_distance(node, u);
        distance_t disb = con_index.get_distance(node, v);

        bool fg = (disa < disb);

        queue<NodeID> q;

        if (con_index.isAncestor(node, u) && con_index.isAncestor(node, v)) {
            node_data[(fg ? v : u)].visstime[tid] = time;
            q.push(fg ? v : u);
        }

        ContractionLabel cu = con_index.labels[node];
        uint16_t cut_level = cu.cut_index.cut_level();
        uint16_t cof = node_data[node].cut_offset;

        size_t szbd = node_data[node].border.size();
        for (size_t sz = node_data[node].distances[32]; sz < szbd; sz++) {
            NodeID bd = node_data[node].border[sz];
            distance_t e = infinity;
            FlatCutIndex &cbd = con_index.labels[bd].cut_index;
            distance_t p = *(cbd.distances() + get_offset(cbd.dist_index(), cut_level) + cof);//con_index.get_distance(node, bd, cof + 1);
            //uint16_t x_offset = get_offset(x.dist_index(), cu.cut_level());
            //if (x.dist_index()[cu.cut_level()] <= x_offset + node_data[cut].cut_offset)continue;
            //distance_t* x_ptr = x.distances() + x_offset + node_data[cut].cut_offset;
            if (node_data[bd].time[32] == etime) {
                e = node_data[bd].distance;
            }
            else if ((fg ? disa : disb) + w + md >= p) {
                continue;
            }
            if (e == infinity) {
                e = con_index.get_distance(bd, (fg ? v : u));
            }
            if (p > (fg ? disa : disb) + w + e) {
                q.push(bd);
                node_data[bd].visstime[tid] = time;
            }
        }

        // �������Ϊ�գ�����
        if (q.empty()) continue;
        // Dijkstra ����
        while (!q.empty()) {
            NodeID next = q.front();
            q.pop();
            if (con_index.is_contracted(next)) continue;
            //tot++;
            // �����ھӽڵ�
            for (const Neighbor& n : node_data[next].neighbors) {
                if (con_index.labels[n.node].cut_index.cut_level() < cut_level) continue;
                if (node_data[n.node].visstime[tid] == time) continue;

                node_data[n.node].visstime[tid] = time;
                FlatCutIndex &cnn = con_index.labels[n.node].cut_index;
                distance_t p = *(cnn.distances() + get_offset(cnn.dist_index(), cut_level) + cof) //con_index.get_distance(cut, n.node, cof + 1)
                , e = infinity;

                if (node_data[n.node].time[32] == etime) {
                    e = node_data[n.node].distance;
                }
                else if ((fg ? disa : disb) + w + md >= p) {
                    continue;
                }
                if (e == infinity) {
                    e = con_index.get_distance(n.node, (fg ? v : u));
                }
                if (p > (fg ? disa : disb) + w + e) {
                    q.push(n.node);
                }
            }
        }
    }
    //cout << tot << endl;
}





void Graph::change_remain_one(
    size_t start, size_t end,
    std::vector<NodeID>& cutt,
    ContractionIndex& con_index,
    distance_t w, distance_t bfw,
    NodeID a, NodeID b,
    const std::vector<CutIndex>& ci,
    std::vector<std::pair<std::pair<NodeID, NodeID>, distance_t>>& change,
    const std::unordered_map<NodeID, distance_t>& get_dis,
    distance_t md, NodeID tid, NodeID etime)
{
    // 1. �ֲ�������
    int tottt = 0;
    const distance_t INF = std::numeric_limits<distance_t>::infinity();

    for (size_t i = start; i < end; i++) {
        NodeID cut = cutt[i];
        int time_now = 2 * (i + 1) + etime * 137;

        // ��ǰ�õ� cut ��Ӧ�� label��level��offset
        const auto& label_cut = con_index.labels[cut];
        uint16_t cut_level = label_cut.cut_index.cut_level();
        uint16_t cof = node_data[cut].cut_offset;

        // ���� far_point �����¾���
        distance_t disa = con_index.get_distance(cut, a);
        distance_t disb = con_index.get_distance(cut, b);
        bool ff = (disa > disb);
        NodeID far_point = ff ? a : b;
        distance_t newfardist = (ff ? disb : disa) + bfw;

        // �ж� ss ����
        NodeID ss = 0;
        if (con_index.isAncestor(cut, far_point) &&
            newfardist == (ff ? disa : disb)) {
            ss = far_point;
            }

            // Radix ��
            using RadixPQ = radix_heap::pair_radix_heap<distance_t, NodeID>;
        RadixPQ q, q2;

        // �ռ��߽��
        std::vector<NodeID> affected_bd;
        distance_t mn = INF;

        // --- �����߽��׶� ---
        {
            distance_t we = md;
            const auto& bd_vec = node_data[cut].border;
            size_t start_idx = node_data[cut].distances[32];
            for (size_t idx = start_idx; idx < bd_vec.size(); idx++) {
                NodeID bd = bd_vec[idx];

                // ȡ ee
                distance_t ee = (node_data[bd].time[32] == etime)
                ? node_data[bd].distance
                : INF;

                // cut_index ���Ӧ�ľ���ָ��
                const FlatCutIndex& x = con_index.labels[bd].cut_index;
                uint16_t x_offset = get_offset(x.dist_index(), cut_level);
                const distance_t* x_ptr = x.distances() + x_offset + cof;

                // ��������Ҫ���µ�
                if (ee == INF && (newfardist + we) > (*x_ptr + con_index.labels[bd].distance_offset))
                    continue;

                // ���� r
                distance_t r = newfardist +
                (ee != INF ? ee
                : con_index.get_distance(bd, far_point));

                if ((*x_ptr + con_index.labels[bd].distance_offset) == r) {
                    // ��ǲ�����
                    node_data[bd].isbdtime[tid] = time_now;
                    node_data[bd].dist[tid] = r + w - bfw;
                    node_data[bd].disttime[tid] = time_now;

                    // �ռ���������С mn
                    affected_bd.push_back(bd);
                    mn = std::min(mn, r + w - bfw);
                }
            }
            // ������� affected_bd Ҳ�� ss������
            if (affected_bd.empty() && ss == 0) continue;
        }

        // --- �� ss �� affected_bd ��ѹ����� q3 ---
        std::queue<NodeID> q3;
        if (ss) {
            q3.push(ss);
            node_data[ss].visstime[tid] = time_now;
        }
        for (NodeID bd : affected_bd) {
            q3.push(bd);
            node_data[bd].visstime[tid] = time_now;
        }

        // --- ��ɢ�׶Σ�ֻ�ö��� q3 ��ǿɴ�(boundary��ɢ) ---
        while (!q3.empty()) {
            NodeID u = q3.front(); q3.pop();
            for (const auto& nei : node_data[u].neighbors) {
                NodeID v = nei.node;
                if (node_data[v].visstime[tid] == time_now
                    || !con_index.isAncestor(cut, v))
                    continue;
                node_data[v].visstime[tid] = time_now;

                distance_t ee = (node_data[v].time[32] == etime)
                ? node_data[v].distance
                : con_index.get_distance(v, far_point);

                const FlatCutIndex& x = con_index.labels[v].cut_index;
                uint16_t x_offset = get_offset(x.dist_index(), cut_level);
                distance_t dist_cut_v = *(x.distances() + x_offset + cof);

                if (dist_cut_v == newfardist + ee) {
                    node_data[v].viss[tid] = true;
                    q3.push(v);
                }
                else {
                    // ���� boundary���Ϳ��� next-phase��ֱ��
                    // �ŵ� q2 �У�distance �� dist_cut_v
                    q2.push(dist_cut_v, v);
                }
            }
        }
        // --- ��һ�׶����� q���� mn �� cut ��ʼ�� ---
        if(affected_bd.size())
        {
            node_data[cut].dist[tid] = 0;
            node_data[cut].disttime[tid] = time_now;
            q.push(mn, cut);
            distance_t disu = std::min(disa, disb);
            //tottt++;
            while (!q.empty()) {
                auto key = q.top_key();
                auto node = q.top_value();
                q.pop();

                if (node_data[node].vistime[tid] == time_now) continue;
                node_data[node].vistime[tid] = time_now;

                // �����ǰ dist �Ѿ����� disu����ֻ���� boundary
                if (node_data[node].dist[tid] > disu) {
                    for (NodeID bd : affected_bd) {
                        distance_t cand = node_data[node].dist[tid]
                        + con_index.get_distance(node, bd);
                        if (cand < node_data[bd].dist[tid]) {
                            node_data[bd].dist[tid] = cand;
                        }
                    }
                    continue;
                }

                // ����������� boundary����Ҳ���� q2
                if (node_data[node].isbdtime[tid] == time_now) {
                    q2.push(node_data[node].dist[tid], node);
                    // �� affected_bd ���Ƴ� node
                    affected_bd.erase(
                        std::remove(affected_bd.begin(), affected_bd.end(), node),
                                      affected_bd.end());
                    if (affected_bd.empty()) break;
                }

                // ������ɢ���ھӣ���������µĸ���·��������� q
                for (const auto& nei : node_data[node].neighbors) {
                    NodeID v = nei.node;
                    distance_t nd = node_data[node].dist[tid] + nei.distance;
                    if (node_data[v].disttime[tid] != time_now
                        || nd < node_data[v].dist[tid]) {
                        node_data[v].dist[tid] = nd;
                    node_data[v].disttime[tid] = time_now;

                    // ͬ������ boundary ���۲����� q
                    distance_t best = INF;
                    for (NodeID bd : affected_bd) {
                        distance_t cand = nd + con_index.get_distance(v, bd);
                        best = std::min(best, cand);
                    }
                    if (best < INF) {
                        q.push(best, v);
                    }
                        }
                }
            }
        }

        // --- �ڶ��׶Σ����� q2����������Ҫ�ĵı߼��� change����������ɢ ---
        while (!q2.empty()) {
            auto key = q2.top_key();
            auto node = q2.top_value();
            q2.pop();

            if (!con_index.isAncestor(cut, node)) continue;
            //if(ss)

            const FlatCutIndex& x = con_index.labels[node].cut_index;
            uint16_t x_offset = get_offset(x.dist_index(), cut_level);
            distance_t dist_cut_node = *(x.distances() + x_offset + cof);

            if (dist_cut_node + w - bfw == key) continue;

            // ������ change
            change.emplace_back(std::make_pair(cut, node), key);

            // ��ɢ�����ھӣ���������
            for (const auto& nei : node_data[node].neighbors) {
                NodeID v = nei.node;
                if (node_data[v].visstime[tid] != time_now
                    || !node_data[v].viss[tid])
                    continue;

                distance_t nd = key + nei.distance;
                if (node_data[v].time[tid] != time_now
                    || nd < node_data[v].distances[tid]) {
                    node_data[v].distances[tid] = nd;
                node_data[v].time[tid] = time_now;
                q2.push(nd, v);
                    }
            }
        }
    }

}


void Graph::update_up_par(size_t start, size_t end, const std::vector<NodeID>& nodes, ContractionIndex& con_index, const vector<pair<pair<NodeID, NodeID>, pair<distance_t,distance_t>>>& edge,const std::vector<CutIndex>& ci,vector<NodeID>& nd_cg)
{

    vector<NodeID> affected_cut;
    int kk = 0,kp=0;
    vector<distance_t>tmp;
    tmp.reserve(200);
    for (size_t i = start; i < end; ++i) {
        NodeID cut = nodes[i];
        uint16_t cof = node_data[cut].cut_offset;
        uint16_t cut_level = con_index.labels[cut].cut_index.cut_level();
        if (node_data[nodes[i]].neighbors.size() == 1)continue;
        distance_t max_cut_bd_distance=0;
        tmp.clear();
        //unordered_map<NodeID, bool>isbd;
        int tf = 0;
        for (NodeID bd:node_data[cut].border) {
            auto &cbd = con_index.labels[bd].cut_index;
            tmp.pb(*(cbd.distances() + get_offset(cbd.dist_index(), cut_level) + cof));//con_index.get_distance(cut, bd,cof+1));
            max_cut_bd_distance = max(max_cut_bd_distance, tmp[tmp.size()-1]);//max_cut_bd_before change.
            //isbd[bd] = true;
        }
        vector<pair<NodeID,bool>>change_sg;
        bool change = false;
        //int tot = 0;
        for (size_t j = 0;j<edge.size();j++) {
            auto &e = edge[j];
            //NodeID nu = edge_fg[j];
            //if(min(Q[node_data[nu].tt].st,Q[node_data[cut].tt].st)>max(Q[node_data[nu].tt].ed , Q[node_data[cut].tt].ed))continue;
            const FlatCutIndex& ci_u = con_index.labels[edge_fg[j]].cut_index;
            const FlatCutIndex& ci_v = con_index.labels[cut].cut_index;
            uint64_t bv_u = *(ci_u.partition_bitvector());
            uint64_t bv_v = *(ci_v.partition_bitvector());
            uint16_t cla = ci_u.cut_level(), cld = ci_v.cut_level();
            // shifting by 64 does not work, so need to check for cla == 0
            if (cla != 0 && ((bv_u ^ bv_v) >> 6 << (64 - cla) != 0)&&cld!=0&& ((bv_u ^ bv_v) >> 6 << (64 - cld) != 0))continue;

            /*
             *           if (!con_index.isAncestor(edge_fg[j], cut) && !con_index.isAncestor(cut, edge_fg[j]))
             *           {
             *               continue;
        }*/
            distance_t dist_cut_a, dist_cut_b,bef_w=e.second.first,w=e.second.second;
            NodeID a = e.first.first;
            NodeID b = e.first.second;
            bool in_g = false;
            if (con_index.isAncestor(cut, a) && con_index.isAncestor(cut, b)) {
                in_g = true;
                dist_cut_a = con_index.get_distance(cut, a, cof + 1);
            }
            else
                dist_cut_a = con_index.get_distance(cut, a);
            if (!in_g&&dist_cut_a > max_cut_bd_distance) {
                continue;
            }//��߲���

            if (!in_g)
                dist_cut_b = con_index.get_distance(cut, b);
            else dist_cut_b = con_index.get_distance(cut, b, cof+1);
            if (dist_cut_a!=dist_cut_b+bef_w&&dist_cut_b!=dist_cut_a+bef_w) {
                continue;
            }
            bool fg1=false;
            if (dist_cut_a > dist_cut_b)fg1 = 1;

            if (!in_g) {
                size_t now = 0;
                if (dist_cut_b > max_cut_bd_distance)continue;
                bool fg2 = false;
                for (NodeID bd : node_data[cut].border) {

                    if (con_index.get_distance(bd,(fg1?a:b))+(fg1?dist_cut_a:dist_cut_b) == tmp[now]) {
                        fg2 = true;
                        break;
                        //change_sg.pb(j);
                    }
                    ++now;
                }
                if (!fg2)continue;
            }
            change_sg.pb({ j,fg1 });
            change = true;
        }

        if (!change)continue;

        if (con_index.is_contracted(cut))continue;
        if (node_data[cut].cut_offset >= 10000)continue;
        int totk = 0;
        int mk =0;
        /*
         *       queue<pair<distance_t, NodeID>>q;
         *       unordered_map<NodeID, bool>vis, vis2;
         *       for (auto eid : change_sg) {
         *           auto e = edge[eid.first];
         *           NodeID to = eid.second ? e.first.first : e.first.second;
         *           q.push({ con_index.get_distance(cut,to),to });
         *           vis[to] = 1;
    }

    priority_queue<SearchNode>q2; unordered_map<NodeID, distance_t>dist;

    while (!q.empty()) {
        NodeID x = q.front().second;
        distance_t ww = q.front().first;
        q.pop();
        totk++;
        if (totk >= mk)break;distance_t distnow = con_index.get_distance(cut, x);
        //if (!con_index.isAncestor(cut,x))continue;
        for (auto e : node_data[x].neighbors) {
            distance_t disw;
            if (e.distance >= (distance_t)200000000) {
                disw = edge[e.distance - (distance_t)200000000].second.first;
    }
    else disw = e.distance;
    distance_t todist = con_index.get_distance(cut, e.node);
        if (!vis[e.node] && distnow+ disw == todist) {
            vis[e.node] = true;
            q.push({ todist, e.node });
    }
    }
    }*/
        //cout << vis.size() << endl;

        //cout << vis.size() << endl;
        /*
         *       if (totk < mk) {
         *
         *           for (auto x : vis) {
         *               if (x.second == true)
         *                   dist[x.first] = infinity;
         *               else {
         *                   distance_t todist = con_index.get_distance(cut, x.first);
         *                   q2.push(SearchNode(todist, x.first));//����Ӱ��ı߽��
         *                   dist[x.first] = todist;
    }
    }
    //cout << q2.size() << endl;
    //cout << q2.size() << endl;

    while (!q2.empty()) {
        NodeID node = q2.top().node;
        distance_t ww = q2.top().distance;
        q2.pop();
        if (vis2[node])continue;
        vis2[node] = 1;
        //update (cut->u)

        if (con_index.isAncestor(cut, node)) {

            //upd.pb({ {cut,node},ww });
            if (con_index.is_contracted(node))continue;
            if (node_data[cut].cut_offset >= 10000)continue;
            FlatCutIndex x = con_index.labels[node].cut_index;
        FlatCutIndex cu = con_index.labels[cut].cut_index;
        uint16_t x_offset = get_offset(x.dist_index(), cu.cut_level());
        if (x.dist_index()[cu.cut_level()] <= x_offset + node_data[cut].cut_offset)continue;
        distance_t* x_ptr = x.distances() + x_offset + node_data[cut].cut_offset;
        *x_ptr = ww;

    }
    //q2.top().distance;
    //cout << cut << " " << node << " " << dist[node] << " "<<vis[node]<<endl;

    for (auto e : node_data[node].neighbors) {
        NodeID v = e.node;
        if (vis[v] != 1)continue;
        distance_t disw;
        if (e.distance >= (distance_t)200000000) {
            disw = edge[(e.distance - 200000000)].second.second;
    }
    else disw = e.distance;
    if (disw + ww < dist[v]) {
        //cout << disw << " " << q2.top().distance << " "<<dist[v]<<endl;
        dist[v] = disw + ww;
        //cout << dist[v] << " " << vis[v] << "?" << endl;
        q2.push(SearchNode(disw + ww, v));

    }
    }
    }
    }
    else
        */
        {
            nd_cg.pb(cut);
            continue;
            /*
             *           priority_queue<SearchNode> q, q2;
             *           distance_t vs = 0;
             *           bool ff = 0;
             *           unordered_map<NodeID, distance_t>dist;
             *           unordered_map<NodeID, bool>vis, isbd;
             *           vector<NodeID>affected_bd;
             *           ContractionLabel cu = con_index.labels[cut];
             *           //cout<< cu.cut_index.cut_level() << endl;
             *           uint16_t cut_level = cu.cut_index.cut_level();
             *           uint16_t cof = node_data[cut].cut_offset;
             *           //cout << cut_level << endl;
             *           //cout << "sd" << endl;
             *
             *           distance_t mx = 0, mn = infinity;
             *           {
             *               for (auto e : edge) {
             *                   NodeID a = e.first.first, b = e.first.second;
             *                   distance_t bfw = e.second.first, w = e.second.second;
             *                   distance_t disa = con_index.get_distance(cut, a), disb = con_index.get_distance(cut, a);
             *                   if (disa > disb)ff = true;
             *                   for (NodeID bd : node_data[cut].border)
             *                   {
             *                       //cout << bd << endl;
             *                       FlatCutIndex x = con_index.labels[bd].cut_index;
             *                       uint16_t x_offset = get_offset(x.dist_index(), cut_level);
             *                       const distance_t* x_ptr = x.distances() + x_offset + cof;
             *                       distance_t r = (ff ? disb : disa) + bfw + con_index.get_distance(bd, (ff ? a : b));
             *                       if (*x_ptr + con_index.labels[bd].distance_offset == r)
             *                       {
             *                           isbd[bd] = true;
             *
        }
        }
        }
        for (NodeID bd : node_data[cut].border) {
            if (isbd[bd]) {
                FlatCutIndex x = con_index.labels[bd].cut_index;
                uint16_t x_offset = get_offset(x.dist_index(), cut_level);
                const distance_t* x_ptr = x.distances() + x_offset + cof;
                mn = min(mn, *x_ptr + con_index.labels[bd].distance_offset);
                affected_bd.push_back(bd);
        }
        else {
            FlatCutIndex x = con_index.labels[bd].cut_index;
            uint16_t x_offset = get_offset(x.dist_index(), cut_level);
            const distance_t* x_ptr = x.distances() + x_offset + cof;
            q2.push({ *x_ptr + con_index.labels[bd].distance_offset,bd });
            dist[bd] = *x_ptr + con_index.labels[bd].distance_offset;
        }
        }
        }
        int w = 0, fg = 0;
        int tot = 0;//return;
        dist[cut] = 0;
        int k = affected_bd.size();
        q.push(SearchNode(mn, cut));
        if (!affected_bd.size())continue;
        //cout << affected_bd.size() << endl;
        //continue;
        while (!q.empty())
        {
        SearchNode next = q.top();
        q.pop();
        if (vis[next.node])continue;
        vis[next.node] = true;
            tot++;
            //if (!con_index.isAncestor(cut, next.node))continue;
            //if (next.distance > mx)break;
            if (isbd[next.node]) {
                auto it = std::find(affected_bd.begin(), affected_bd.end(), next.node);
                isbd[next.node] = false;
                vector<NodeID>change;
                affected_bd.erase(it);
                q2.push({ dist[next.node],next.node });
                while (!q.empty())
                {
                NodeID nd = q.top().node;
                q.pop();
                change.push_back(nd);
        }
        for (NodeID x : change) {
            distance_t f = infinity;
            for (auto bd : affected_bd) {
                f = min(f, dist[x] + con_index.get_distance(x, bd));
        }
        q.push({ f,x });
        }
        }
        if (!affected_bd.size())break;
        /*
         *               if (1) {
         *                   bool pp = false;
         *                   for (auto bd : affected_bd) {
         *                       if (bd == next.node) {
         *                           pp = true;
         *                           auto it = std::find(affected_bd.begin(), affected_bd.end(), bd);
         *                           affected_bd.erase(it);
         *                           q2.push({ next.distance,bd });
         *                           break;
        }
        FlatCutIndex x = con_index.labels[bd].cut_index;
        uint16_t x_offset = get_offset(x.dist_index(), cut_level);
        const distance_t* x_ptr = x.distances() + x_offset + cof;
        if (next.distance + con_index.get_distance(next.node, bd) <= *x_ptr + w - bfw)pp = true;
        }
        if (!pp)continue;
        }
        */

        //˵������Ӱ��
        //if (!con_index.isAncestor(cut, next.node) && next.distance == con_index.get_distance(cut,next.node))continue;
        /*
         *               for (Neighbor n : node_data[next.node].neighbors)
         *               {
         *                   // update distance and enque
         *                   distance_t disw;
         *                   if (n.distance >= (distance_t)200000000) {
         *                       disw = edge[(n.distance - 200000000)].second.second;
        }
        else disw = n.distance;
        distance_t new_dist = dist[next.node] + disw, f = infinity;
            if (!dist.count(n.node) || new_dist < dist[n.node])
            {
            dist[n.node] = new_dist;
            for (auto bd : affected_bd) {
                f = min(f, dist[n.node] + con_index.get_distance(n.node, bd));
        }
        q.push(SearchNode(f, n.node));
        }
        }
        }
        kp++;
        kk+=tot;
        continue;

        {
        int tott = 0;
        unordered_map<NodeID, bool>vis, vis2;
        unordered_map<NodeID, distance_t>dist;
        q2.push(SearchNode(0, cut));
        while (!q2.empty()) {
            NodeID node = q2.top().node;
            distance_t ww = q2.top().distance;
            q2.pop();
            tott++;
            if (vis2[node])continue;
            vis2[node] = true;
            /*
             *                   if (!con_index.isAncestor(cut, node)) {
             *                       continue;
        }
        */
            /*
             *                   if (tott >= 200000)
             *                   {
             *                       cout << "o"<<endl;
             *                       cout << border_num << " " << bdsz << endl;
             *                       continue;
        }
        */
            /*
             *                   if (!con_index.isAncestor(cut, node)) {
             *                       continue;
        }
        */
            /*
             *                   if (con_index.isAncestor(cut, node)) {
             *                       //upd.pb({ {cut,node},ww });
             *
             *                       if (con_index.is_contracted(node))continue;
             *                       if (node_data[cut].cut_offset >= 10000)continue;
             *                       FlatCutIndex x = con_index.labels[node].cut_index;
             *                       FlatCutIndex cu = con_index.labels[cut].cut_index;
             *                       uint16_t x_offset = get_offset(x.dist_index(), cu.cut_level());
             *                       if (x.dist_index()[cu.cut_level()] <= x_offset + node_data[cut].cut_offset)continue;
             *                       distance_t* x_ptr = x.distances() + x_offset + node_data[cut].cut_offset;
             *x_ptr = ww;
             *
             *
        }


        //cout << cut << " " << node << " " << dist[node] << " "<<vis[node]<<endl;
        //update (cut->u)
        for (auto e : node_data[node].neighbors) {
            NodeID v = e.node;
            distance_t disw;
            FlatCutIndex vv = con_index.labels[v].cut_index;
            if (vv.cut_level() > cut_level)continue;
            if (e.distance >= (distance_t)200000000) {
                disw = edge[(e.distance - 200000000)].second.second;
        }
        else disw = e.distance;
        if (!dist.count(v) || disw + ww < dist[v]) {
            //cout << disw << " " << q2.top().distance << " "<<dist[v]<<endl;
            dist[v] = disw + ww;
            //cout << dist[v] << " " << vis[v] << "?" << endl;
            q2.push(SearchNode(disw + ww, v));

        }
        }
        }*/
            //cout << tott << endl;
        }
    }
}
//cout << q2.size() << endl;

// cout << 1;*/

            //cout << kp << " " << kk << endl;

            struct pair_hash {
                std::size_t operator()(const std::pair<NodeID, NodeID>& p) const {
                    return std::hash<NodeID>()(p.first) ^ (std::hash<NodeID>()(p.second) << 1);
                }
            };
            void Graph::process_nodes(size_t start, size_t end, const std::vector<NodeID>& nodes, ContractionIndex& con_index,
                                      NodeID u, NodeID v, distance_t w, distance_t before_w,
                                      const std::vector<CutIndex>& ci,const std::unordered_map<NodeID,distance_t> &get_dis,distance_t md, vector<NodeID>& upd, vector<pair<pair<NodeID, NodeID>, distance_t>> &change,NodeID etime) {
                int TT = 0; //util::start_timer();
                int kk = 0;
                {
                    for (size_t i = start; i < end; ++i) {
                        NodeID node = nodes[i];
                        if (node_data[node].cut_offset >= 10000)continue;
                        if (con_index.is_contracted(node))continue;
                        bool p = false;
                        if (con_index.isAncestor(node, u) && con_index.isAncestor(node, v)) {
                            p = true;
                        }//continue;
                        //continue;

                        if (!p) {

                            if (node_data[node].time[32]!=etime)
                            {
                                if (!node_data[node].mxbordis)continue;
                                //cout << "sd" << endl;
                                if (md + w >= node_data[node].mxbordis)
                                    continue;
                                if (con_index.get_distance(node, u) >= node_data[node].mxbordis)continue;

                            }
                            else if (node_data[node].distance + w >= node_data[node].mxbordis)
                            {
                                continue;
                            }
                        }

                        //continue;
                        //continue;
                        //cout << w << " " << before_w << endl;
                        //cout << w << "  " << before_w << endl;
                        if (w < before_w) {
                            //upd.push_back(node); continue;
                            ContractionLabel ccu = con_index.labels[u];
                            distance_t disa = con_index.get_distance(node, u), disb = con_index.get_distance(node, v);
                            bool fg = false; if (disa < disb)fg = true;
                            NodeID cut = node;
                            vector<NodeID>affected_poi;
                            //cout << "sd" << endl;
                            //return;
                            //    util::start_timer();
                            unordered_map<NodeID, distance_t>dist2;
                            //TT++;
                            //continue;
                            if (p)
                            {
                                affected_poi.push_back((fg ? v : u));
                                continue;
                            }

                            //return;

                            ContractionLabel cu=con_index.labels[node];
                            uint16_t cut_level = cu.cut_index.cut_level();
                            uint16_t cof = node_data[node].cut_offset;
                            //continue;
                            //TT += node_data[node].border.size();
                            //continue;
                            int bdsz = node_data[node].border.size();
                            for (size_t sz=0;sz<bdsz;sz++)
                            {
                                NodeID bd = node_data[node].border[sz];
                                distance_t e = infinity;
                                FlatCutIndex x = con_index.labels[bd].cut_index;
                                uint16_t x_offset = get_offset(x.dist_index(), cut_level);
                                const distance_t* x_ptr = x.distances() + x_offset + cof;
                                //continue;
                                if (node_data[bd].time[32]==etime) { e=node_data[bd].distance; }
                                else if ((fg ? disa : disb) + w + md >= *x_ptr + con_index.labels[bd].distance_offset)
                                {
                                    //cout << "sd" << endl;
                                    continue;
                                }
                                if (e==infinity) {
                                    e = con_index.get_distance(bd, (fg ? v : u));
                                    /*
                                     *                   node_data[bd].time[32] = etime;
                                     *                   node_data[bd].distance = e;
                                     */
                                }

                                if (*x_ptr + con_index.labels[bd].distance_offset > (fg ? disa : disb) + w + e)
                                {
                                    node_data[node].distances[32] = sz;
                                    affected_poi.push_back(bd); break;
                                }
                            }
                            //continue;
                            //cout << util::stop_timer() << endl;
                            //cout << affected_poi.size() << endl;
                            if (!affected_poi.size())continue;
                            upd.push_back(node); continue;
                            NodeID a, b;

                            //run_dij_upd_down(con_index, ci, affected_poi, u, v, w, node,change);
                            //else
                            //run_dij_upd_down(con_index, ci, affected_poi, v, u, w, node,change);
                            //cout << upd.size() << endl;
                            //cout << util::stop_timer() << endl;

                            priority_queue<SearchNode> q;
                            for (NodeID x : affected_poi) {
                                q.push(SearchNode(dist2[x], x));
                            }
                            distance_t mb = 0;

                            //return;
                            // init queue
                            // dijkstra
                            //if (cut_level ) { cout << "?"; exit(0); }
                            //return;
                            int tot = 0;
                            //cout << cof << endl;

                            while (!q.empty())
                            {
                                SearchNode next = q.top();
                                q.pop();
                                if (con_index.is_contracted(next.node))continue;
                                if (!con_index.isAncestor(cut, next.node))continue;
                                if (con_index.get_distance(cut, next.node,cof+1) <= node_data[next.node].distance)continue;
                                tot++;
                                change.pb({ {cut,next.node},dist2[next.node] });
                                {
                                    /*
                                     *                       //upd.pb({ {cut,next.node},next.distance });
                                     *                       FlatCutIndex x = con_index.labels[next.node].cut_index;
                                     *                       uint16_t x_offset = get_offset(x.dist_index(), cut_level);
                                     *                       //if (cof + x_offset >= x.dist_index()[cut_level])
                                     *                       // if (cut_level && x.dist_index()[cut_level] < x.dist_index()[cut_level - 1])continue;
                                     *                       //cout << x.dist_index()[cut_level]<<" "<< x_offset << " " << cof << " " << x.dist_index()[cut_level] - x_offset << " !" << cut << " " << next.node << " " << x.distances() + x_offset + cof << " " << a << "  " << b << " " << cut_level << endl;
                                     *                       if (cof + x_offset < x.dist_index()[cut_level] && cof < 10000) {
                                     *                           //x.distances();
                                     *                           distance_t* x_ptr = x.distances() + x_offset + cof;
                                     *x_ptr = next.distance;
                                }
                                */
                                }
                                //upd.pb({ {cut,next.node},node_data[next.node].distance});
                                for (Neighbor n : node_data[next.node].neighbors)
                                {
                                    // filter neighbors nodes not belonging to subgraph or having higher landmark level
                                    if (con_index.labels[n.node].cut_index.cut_level() > cut_level||(con_index.labels[n.node].cut_index.cut_level() == cut_level && node_data[n.node].cut_offset<cof))
                                        continue;
                                    // update distance and enque
                                    distance_t new_dist = next.distance + n.distance;
                                    if (!dist2.count(n.node) || new_dist < dist2[n.node])
                                    {
                                        dist2[n.node] = new_dist;
                                        q.push(SearchNode(new_dist, n.node));
                                    }
                                }
                            }
                        }
                        else {//cout << "sd" << endl;

                            int to = 0;
                            //for (NodeID bd : node_data[node].border)if (con_index.get_distance(node, bd) == (fg ? disa : disb) + before_w + con_index.get_distance(bd, (fg ? v : u)))to=1;
                            //if (!to)continue;
                            NodeID a, b;
                            distance_t disa = con_index.get_distance(node, u), disb = con_index.get_distance(node, v);
                            NodeID cut = node; bool fg = false; if (disa < disb)fg = true;
                            if (fg)a = u, b = v; else a = v, b = u;
                            distance_t bfw = before_w;
                            queue<SearchNode> q;
                            unordered_map<NodeID, bool>vis;
                            unordered_map<NodeID, uint8_t>type;
                            int fz = 0;int tott = 1;//return;
                            //cout << fz << endl;
                            /*
                             *               q.push(SearchNode(con_index.get_distance(cut, a) + bfw, b));
                             *
                             *               type[b] = 3;
                             *
                             *               // init queue
                             *               // dijkstra
                             *               //cout << b << endl;
                             *
                             *               while (!q.empty())
                             *               {
                             *                   SearchNode next = q.front();
                             *                   q.pop();
                             *                   tott++;
                             *                   if (tott > fz)break;
                             *                   type[next.node] |= 1;
                             *
                             *                   for (Neighbor n : node_data[next.node].neighbors)
                             *                   {
                             *                       type[n.node] |= 2;
                             *                       // update distance and enque
                             *                       distance_t new_dist = next.distance + n.distance;
                             *                       //cout << con_index.get_distance(cut, next.node) + n.distance << " " << con_index.get_distance(cut, n.node) << endl;
                             *                       if (!vis[n.node] && con_index.get_distance(cut, next.node) + n.distance == con_index.get_distance(cut, n.node))
                             *                       {
                             *                           //cout << n.node << endl;
                             *                           vis[n.node] = 1;
                             *                           q.push(SearchNode(new_dist, n.node));
                        }
                        }
                        }
                        */
                            //cout << tott << "?" << endl;
                            //exit(0);
                            //return;
                            //cout << tott <<endl;
                            if (tott > fz) {
                                //continue;


                                NodeID cut = node;
                                priority_queue<SearchNode> q, q2;
                                distance_t vs = 0;
                                bool ff = 0;
                                unordered_map<NodeID, distance_t>dist;
                                unordered_map<NodeID, bool>vis, isbd;
                                vector<NodeID>affected_bd;
                                //distance_t disa = con_index.get_distance(cut, a), disb = con_index.get_distance(cut, b);
                                if (disa > disb)ff = true;
                                ContractionLabel cu = con_index.labels[cut];
                                //cout<< cu.cut_index.cut_level() << endl;
                                uint16_t cut_level = cu.cut_index.cut_level();
                                uint16_t cof = node_data[cut].cut_offset;
                                //cout << cut_level << endl;
                                //cout << "sd" << endl;
                                distance_t fardist = (ff ? disb : disa) + bfw;
                                distance_t mx = 0, mn = infinity;
                                bool fe = 0;
                                {
                                    //continue;
                                    distance_t we,ee;
                                    size_t alsz = node_data[node].border.size();
                                    for (size_t sz=0;sz<alsz;sz++)
                                    {
                                        NodeID bd = node_data[node].border[sz];
                                        ee = infinity;
                                        if (node_data[bd].time[32]!=etime)
                                        {
                                            we = md;
                                        }
                                        else { we = node_data[bd].distance; ee = we; }
                                        //cout << bd << endl;
                                        FlatCutIndex x = con_index.labels[bd].cut_index;
                                        uint16_t x_offset = get_offset(x.dist_index(), cut_level);
                                        const distance_t* x_ptr = x.distances() + x_offset + cof;//continue;
                                        if ((ff ? disb : disa) + bfw +we> *x_ptr + con_index.labels[bd].distance_offset)
                                        {
                                            //cout << "o" << endl;
                                            continue;
                                        }
                                        distance_t r = fardist;
                                        if (ee==infinity) r += con_index.get_distance(bd, (ff ? a : b));
                                        else r += ee;
                                        if (*x_ptr + con_index.labels[bd].distance_offset == r)
                                        {
                                            node_data[node].distances[32] = sz;
                                            //affected_bd.push_back(bd);
                                            fe = true;
                                            break;
                                            /*
                                             *                               isbd[bd] = true;
                                             *                               mn = min(mn, *x_ptr + con_index.labels[bd].distance_offset);
                                             *
                                             *                               dist[bd] = *x_ptr + con_index.labels[bd].distance_offset + w - bfw;
                                             */
                                        }

                                    }
                                    if (!fe)
                                    {
                                        //cout << "!" << endl;
                                        continue;
                                    }
                                    {
                                        upd.pb(node); continue;
                                    }
                                    //cout << "?" << endl;
                                    for (NodeID bd : node_data[cut].border)
                                    {
                                        FlatCutIndex x = con_index.labels[bd].cut_index;
                                        uint16_t x_offset = get_offset(x.dist_index(), cut_level);
                                        const distance_t* x_ptr = x.distances() + x_offset + cof;
                                        if (!isbd[bd])q2.push({ *x_ptr + con_index.labels[bd].distance_offset,bd });
                                    }
                                }
                                q.push(SearchNode(mn, cut));
                                //return;
                                //cout << affected_bd.size() << " "<< node_data[cut].border.size()<<" "<<mx<<" "<<node_data[cut].mxbordis<<endl;
                                {
                                    //cout << "?" << endl;
                                    //if (affected_bd.size() <= 50)return;
                                    //else cout <<  cut_level << endl;
                                    //else cout << affected_bd.size() << endl;
                                    //return;
                                    distance_t disu = min(disa, disb);
                                    int w = 0, fg = 0;
                                    int tot = 0; //continue;
                                    while (!q.empty())
                                    {
                                        SearchNode next = q.top();
                                        q.pop();
                                        if (vis[next.node])continue;
                                        vis[next.node] = true;
                                        //tot++;
                                        if (dist[next.node] > mx)break;
                                        if (dist[next.node] > disu) {
                                            distance_t f = infinity, g = 0;
                                            for (auto bd : affected_bd) {
                                                g = dist[next.node] + con_index.get_distance(next.node, bd);
                                                if (g < dist[bd])
                                                {
                                                    dist[bd] = min(dist[bd], g);
                                                }
                                            }
                                            continue;
                                        }
                                        if (isbd[next.node]) {
                                            auto it = std::find(affected_bd.begin(), affected_bd.end(), next.node);
                                            isbd[next.node] = false;
                                            q2.push({ dist[next.node],next.node });
                                            vector<NodeID>change;
                                            affected_bd.erase(it);
                                            while (!q.empty())
                                            {
                                                NodeID nd = q.top().node;
                                                q.pop();
                                                change.push_back(nd);
                                            }
                                            for (NodeID x : change) {
                                                distance_t f = infinity, g = 0;
                                                for (auto bd : affected_bd) {
                                                    g = dist[x] + con_index.get_distance(x, bd);
                                                    if (g < dist[bd])
                                                    {
                                                        f = min(f, g);
                                                    }
                                                }
                                                if (f != infinity)q.push({ f,x });
                                            }
                                        }
                                        for (Neighbor n : node_data[next.node].neighbors)
                                        {
                                            // update distance and enque
                                            distance_t new_dist = dist[next.node] + n.distance, f = infinity, g = 0;
                                            if (!dist.count(n.node) || new_dist < dist[n.node])
                                            {
                                                dist[n.node] = new_dist;
                                                for (auto bd : affected_bd) {
                                                    g = dist[n.node] + con_index.get_distance(n.node, bd);
                                                    if (g < dist[bd])
                                                    {
                                                        f = min(f, g);
                                                    }
                                                }
                                                if (f != infinity)
                                                    q.push(SearchNode(f, n.node));
                                            }
                                        }
                                    }
                                    //cout << tot << endl;
                                    for (NodeID x : affected_bd) {
                                        q2.push({ dist[x],x });

                                    }
                                    //cout << "sd" <<  endl;
                                    //cout << tot << endl;
                                    //cout << tot <<" "<< affected_bd.size() << " " << node_data[cut].border.size() << " " << mx << " " << node_data[cut].mxbordis <<  endl; //return;
                                    vis.clear();
                                    dist.clear();
                                    q2.push({ 0,cut });
                                    while (!q2.empty())
                                    {
                                        SearchNode next = q2.top();
                                        q2.pop();
                                        if (vis[next.node])continue;
                                        vis[next.node] = true;
                                        if (!con_index.isAncestor(cut, next.node))continue;
                                        change.pb({ {cut,next.node},dist[next.node] });//cout << "!" << endl;
                                        //˵������Ӱ��
                                        //if (!con_index.isAncestor(cut, next.node) && next.distance == con_index.get_distance(cut,next.node))continue;
                                        for (Neighbor n : node_data[next.node].neighbors)
                                        {
                                            // update distance and enque
                                            distance_t new_dist = next.distance + n.distance;
                                            ContractionLabel nn = con_index.labels[n.node];
                                            //cout<< cu.cut_index.cut_level() << endl;
                                            if (nn.cut_index.cut_level() > cut_level)continue;
                                            if (!dist.count(n.node) || new_dist < dist[n.node])
                                            {
                                                dist[n.node] = new_dist;
                                                q2.push(SearchNode(new_dist, n.node));
                                            }
                                        }
                                    }
                                    //cout << "!" << tott << endl;

                                    //if (!g)cout << tot << endl;

                                }

                                //
                            }
                            else {
                                priority_queue<SearchNode> q;
                                //if (type.size())
                                //  cout << type.size() << endl;
                                unordered_map<NodeID, distance_t>dist;
                                unordered_map<NodeID, bool>vis;
                                //cout << fg_Node.size() << endl;
                                for (auto x : type) {
                                    if (x.second == 2) {
                                        NodeID border = x.first;
                                        dist[border] = con_index.get_distance(cut, border);
                                        q.push(SearchNode(dist[border], border));
                                    }
                                }
                                //return;
                                q.push({ 0,cut });
                                while (!q.empty())
                                {
                                    SearchNode next = q.top();
                                    q.pop();
                                    if (type[next.node] == 3 && con_index.isAncestor(cut, next.node)) {
                                        change.pb({ { cut,next.node },dist[next.node] });
                                    }
                                    if (next.distance >= node_data[next.node].mxbordis + w - bfw && !con_index.isAncestor(cut, next.node))continue;
                                    for (Neighbor n : node_data[next.node].neighbors)
                                    {
                                        if (type[n.node] != 3)continue;
                                        // update distance and enque
                                        distance_t new_dist = next.distance + n.distance;
                                        if (!dist.count(n.node) || new_dist < dist[n.node])
                                        {
                                            dist[n.node] = new_dist;
                                            q.push(SearchNode(new_dist, n.node));
                                        }
                                    }
                                }
                            }


                        }

                        /*
                         *           if (w < before_w) {
                         *               if (dist_node_u + w < dist_node_v) {
                         *                   local_affected_set1.push_back(node);
                    }
                    else if (dist_node_v + w < dist_node_u) {
                        local_affected_set2.push_back(node);
                    }
                    }
                    else {
                        if (dist_node_u + before_w == dist_node_v) {
                            local_affected_set1.push_back(node);
                    }
                    else if (dist_node_v + before_w == dist_node_u) {
                        local_affected_set2.push_back(node);
                    }
                    }*/
                    }
                }//cout << kk << endl;
                //cout << util::stop_timer() << endl;
                //cout << TT << endl;
                                      }

                                      void Graph::get_cut(int id, int level,
                                                          const ContractionIndex& con_index, vector<NodeID> upd, vector<NodeID> edge,
                                                          const vector<pair<pair<NodeID, NodeID>, pair<distance_t, distance_t>>>& E, vector<NodeID>& ans) {

                                          // Use stack to simulate the recursion
                                          struct StackFrame {
                                              int id;
                                              int level;
                                              vector<NodeID> upd;
                                              vector<NodeID> edge;
                                              vector<NodeID> edgL, edgR;
                                              size_t edge_index;
                                          };

                                          // Stack to store states of the "recursive" call
                                          stack<StackFrame> stack;
                                          stack.push({ id, level, upd, edge, {}, {}, 0 });

                                          while (!stack.empty()) {
                                              // Pop the current state from the stack
                                              StackFrame frame = stack.top();
                                              stack.pop();

                                              int current_id = frame.id;
                                              int current_level = frame.level;
                                              vector<NodeID>& current_upd = frame.upd;
                                              vector<NodeID>& current_edge = frame.edge;
                                              vector<NodeID>& current_edgL = frame.edgL;
                                              vector<NodeID>& current_edgR = frame.edgR;
                                              size_t edge_index = frame.edge_index;

                                              // Process all edges in the current "level"
                                              while (edge_index < current_edge.size()) {
                                                  NodeID edg = current_edge[edge_index];
                                                  edge_index++;

                                                  NodeID u = E[edg].first.first;
                                                  NodeID v = E[edg].first.second;
                                                  distance_t w = E[edg].second.first;
                                                  uint64_t bv = con_index.labels[u].cut_index.partition();
                                                  int cut_level = min(con_index.labels[v].cut_index.cut_level(),
                                                                      con_index.labels[u].cut_index.cut_level());
                                                  int fg = 0;

                                                  // Level check
                                                  if (current_level == cut_level) {
                                                      ans[edg] = current_id;
                                                      continue;
                                                  }

                                                  // Check for cut conditions
                                                  for (size_t j = 0; j < Q[current_id].cut.size(); j++) {
                                                      NodeID bd = Q[current_id].cut[j];
                                                      int cof2 = node_data[bd].cut_offset;
                                                      for (size_t i = 0; i < current_upd.size(); i++) {
                                                          NodeID x = current_upd[i];
                                                          if (x == bd) continue;
                                                          int cof1 = node_data[x].cut_offset;
                                                          if (con_index.get_distance(x, bd, cof1 + 1) ==
                                                              con_index.get_distance(x, u, cof1 + 1) + w + con_index.get_distance(bd, v, cof2 + 1)) {
                                                              fg = 1; break;
                                                              }
                                                              if (con_index.get_distance(x, bd, cof1 + 1) ==
                                                                  con_index.get_distance(x, v, cof1 + 1) + w + con_index.get_distance(bd, u, cof2 + 1)) {
                                                                  fg = 1; break;
                                                                  }
                                                      }
                                                  }

                                                  // If a valid edge is found
                                                  if (fg) {
                                                      ans[edg] = current_id;
                                                      continue;
                                                  }

                                                  // Push children to stack for further processing
                                                  if (1 & (bv >> current_level)) {
                                                      if (Q[current_id].rchild) {
                                                          current_edgR.push_back(edg);
                                                      }
                                                      else {
                                                          ans[edg] = current_id;
                                                      }
                                                  }
                                                  else {
                                                      if (Q[current_id].lchild) {
                                                          current_edgL.push_back(edg);
                                                      }
                                                      else {
                                                          ans[edg] = current_id;
                                                      }
                                                  }
                                              }

                                              // If there are any edges in the right child, push the state for the right child
                                              if (Q[current_id].rchild && !current_edgR.empty()) {
                                                  stack.push({ Q[current_id].rchild, current_level + 1, current_upd, current_edgR, {}, {}, 0 });
                                              }

                                              // If there are any edges in the left child, push the state for the left child
                                              if (Q[current_id].lchild && !current_edgL.empty()) {
                                                  stack.push({ Q[current_id].lchild, current_level + 1, current_upd, current_edgL, {}, {}, 0 });
                                              }
                                          }
                                                          }


                                                          void Graph::get_cutdis(size_t start,size_t end,vector<NodeID> &C,vector<NodeID> &CALL, const ContractionIndex& con_index,NodeID u,NodeID v,distance_t w,NodeID &ans) {
                                                              int ft = 0; size_t p = 0;

                                                              const FlatCutIndex& uu = con_index.labels[u].cut_index;
                                                              const FlatCutIndex& vv = con_index.labels[v].cut_index;
                                                              //uint16_t x_offset = get_offset(x.dist_index(), cu.cut_level());
                                                              //if (x.dist_index()[cu.cut_level()] <= x_offset + node_data[cut].cut_offset)continue;
                                                              //distance_t* x_ptr = x.distances() + x_offset + node_data[cut].cut_offset;
                                                              for (size_t sz = 0; sz < C.size(); sz++) {
                                                                  NodeID Cid = C[sz];
                                                                  bool fg = 0;
                                                                  //cout << Cid << endl;
                                                                  p += Q[Cid].cut.size();
                                                                  for (size_t i = start; i < min(p, end); i++) {
                                                                      NodeID cuta = CALL[i];
                                                                      int cof1 = node_data[cuta].cut_offset;
                                                                      const FlatCutIndex& ca = con_index.labels[cuta].cut_index;
                                                                      uint16_t cuta_level = ca.cut_level(), cutb_level;
                                                                      distance_t disu = *(uu.distances() + get_offset(uu.dist_index(), cuta_level) + cof1);
                                                                      distance_t disv = *(vv.distances() + get_offset(vv.dist_index(), cuta_level) + cof1);
                                                                      for (NodeID cutb : Q[Cid].cut) {
                                                                          if (!con_index.isAncestor(cuta, cutb))continue;
                                                                          const FlatCutIndex& cb = con_index.labels[cutb].cut_index;
                                                                          cutb_level = cb.cut_level();
                                                                          int cof2 = node_data[cutb].cut_offset;
                                                                          if ((*(cb.distances() + get_offset(cb.dist_index(), cuta_level) + cof1)) == disu + w + *(vv.distances() + get_offset(vv.dist_index(), cutb_level) + cof2)) {
                                                                              ans = sz;
                                                                              return;
                                                                          }
                                                                          /*if (con_index.get_distance(cuta, cutb, cof1 + 1) >
                                                                           *                   con_index.get_distance(cuta, u, cof1 + 1) + w + con_index.get_distance(cutb, v, cof2 + 1))
                                                                           *               {
                                                                           *                   ans = sz;
                                                                           *                   return;
                                                                      }
                                                                      */
                                                                          if ((*(cb.distances() + get_offset(cb.dist_index(), cuta_level) + cof1)) == disv + w + *(uu.distances() + get_offset(uu.dist_index(), cutb_level) + cof2)) {
                                                                              ans = sz;
                                                                              return;
                                                                          }
                                                                          /*
                                                                           *               if (con_index.get_distance(cuta, cutb, cof1 + 1) >
                                                                           *                   con_index.get_distance(cuta, v, cof1 + 1) + w + con_index.get_distance(cutb, u, cof2 + 1))
                                                                           *               {
                                                                           *                   ans = sz;
                                                                           *                   return;
                                                                      }
                                                                      */
                                                                      }
                                                                  }
                                                              }
                                                              ans = C.size() - 1;
                                                              return;
                                                          }

                                                          void Graph::get_cutdis2(size_t start, size_t end, vector<NodeID>& C, vector<NodeID>& CALL, const ContractionIndex& con_index, NodeID u, NodeID v, distance_t w, NodeID& ans) {
                                                              int ft = 0; size_t p = 0;

                                                              const FlatCutIndex& uu = con_index.labels[u].cut_index;
                                                              const FlatCutIndex& vv = con_index.labels[v].cut_index;
                                                              //uint16_t x_offset = get_offset(x.dist_index(), cu.cut_level());
                                                              //if (x.dist_index()[cu.cut_level()] <= x_offset + node_data[cut].cut_offset)continue;
                                                              //distance_t* x_ptr = x.distances() + x_offset + node_data[cut].cut_offset;
                                                              for (size_t sz = 0; sz < C.size(); sz++) {
                                                                  NodeID Cid = C[sz];
                                                                  bool fg = 0;
                                                                  //cout << Cid << endl;
                                                                  p += Q[Cid].cut.size();
                                                                  for (size_t i = start; i < min(p, end); i++) {
                                                                      NodeID cuta = CALL[i];
                                                                      int cof1 = node_data[cuta].cut_offset;
                                                                      const FlatCutIndex &ca = con_index.labels[cuta].cut_index;
                                                                      uint16_t cuta_level = ca.cut_level(),cutb_level;
                                                                      distance_t disu = *(uu.distances() + get_offset(uu.dist_index(), cuta_level) + cof1);
                                                                      distance_t disv = *(vv.distances() + get_offset(vv.dist_index(), cuta_level) + cof1);
                                                                      for (NodeID cutb : Q[Cid].cut) {
                                                                          if (!con_index.isAncestor(cuta, cutb))continue;
                                                                          const FlatCutIndex &cb = con_index.labels[cutb].cut_index;
                                                                          cutb_level = cb.cut_level();
                                                                          int cof2 = node_data[cutb].cut_offset;
                                                                          if ((*(cb.distances() + get_offset(cb.dist_index(), cuta_level) + cof1)) > disu + w + *(vv.distances() + get_offset(vv.dist_index(), cutb_level) + cof2)) {
                                                                              ans = sz;
                                                                              return;
                                                                          }
                                                                          /*if (con_index.get_distance(cuta, cutb, cof1 + 1) >
                                                                           *                   con_index.get_distance(cuta, u, cof1 + 1) + w + con_index.get_distance(cutb, v, cof2 + 1))
                                                                           *               {
                                                                           *                   ans = sz;
                                                                           *                   return;
                                                                      }
                                                                      */
                                                                          if ((*(cb.distances() + get_offset(cb.dist_index(), cuta_level) + cof1)) > disv + w + *(uu.distances() + get_offset(uu.dist_index(), cutb_level) + cof2)) {
                                                                              ans = sz;
                                                                              return;
                                                                          }
                                                                          /*
                                                                           *               if (con_index.get_distance(cuta, cutb, cof1 + 1) >
                                                                           *                   con_index.get_distance(cuta, v, cof1 + 1) + w + con_index.get_distance(cutb, u, cof2 + 1))
                                                                           *               {
                                                                           *                   ans = sz;
                                                                           *                   return;
                                                                      }
                                                                      */
                                                                      }
                                                                  }
                                                              }
                                                              ans = C.size() - 1;
                                                              return;
                                                          }

                                                          unordered_map<NodeID, int>m;
                                                          void Graph::update_edge(const vector<CutIndex>& ci, ContractionIndex &con_index, NodeID u, NodeID v, distance_t w, NodeID etime) {
                                                              distance_t before_w, tmp = infinity;
                                                              if (con_index.isAncestor(v, u)) swap(u, v);
                                                              for (const auto& n : node_data[u].neighbors) {
                                                                  if (n.node == v) {
                                                                      tmp = std::min(tmp, n.distance);
                                                                  }
                                                              }
                                                              if (tmp == infinity) return;
                                                              before_w = tmp;
                                                              if (w == before_w) return;

                                                              add_edge(u, v, w, true, true);

                                                              int d = 0, tt = 0;
                                                              unordered_map<NodeID, distance_t> get_dis;
                                                              vector<NodeID> upd[35];
                                                              vector<pair<pair<NodeID, NodeID>, distance_t>> change[35];

                                                              // submit find_affected_point as fire-and-forget (same semantics as detached thread)
                                                              distance_t md = 0;
                                                              auto f_md = global_pool().submit(&Graph::find_affected_point, this,
                                                                                               std::cref(con_index), std::ref(md),
                                                                                               u, v, w, before_w, etime);

                                                              NodeID id;
                                                              vector<NodeID> C, CALL;
                                                              NodeID nowid = Q[0].lchild;
                                                              int level = -1;
                                                              while (1) {
                                                                  if (!con_index.isAncestor(nowid, u)) break;
                                                                  level++;
                                                                  C.pb(nowid);
                                                                  for (NodeID xx : Q[nowid].cut) {
                                                                      CALL.pb(xx);
                                                                  }
                                                                  uint64_t bv = con_index.labels[u].cut_index.partition();
                                                                  if (1 & (bv >> level)) {
                                                                      if (Q[nowid].rchild) nowid = Q[nowid].rchild;
                                                                      else break;
                                                                  } else {
                                                                      if (Q[nowid].lchild) nowid = Q[nowid].lchild;
                                                                      else break;
                                                                  }
                                                              }

                                                              size_t num_threads = 31;
                                                              NodeID to = 59;

                                                              if (w > before_w) {
                                                                  vector<NodeID> nd; nd.resize(num_threads + 1);
                                                                  size_t chunk_size = (CALL.size() + num_threads - 1) / num_threads;

                                                                  // parallel get_cutdis
                                                                  {
                                                                      auto& pool = global_pool();
                                                                      std::vector<std::future<void>> futs;
                                                                      futs.reserve(num_threads);
                                                                      for (size_t t = 0; t < num_threads; ++t) {
                                                                          size_t start = t * chunk_size;
                                                                          size_t end = std::min(start + chunk_size, CALL.size());
                                                                          futs.emplace_back(
                                                                              pool.submit(&Graph::get_cutdis, this,
                                                                                          start, end, std::ref(C), std::ref(CALL), std::ref(con_index),
                                                                                          u, v, before_w, std::ref(nd[t]))
                                                                          );
                                                                      }
                                                                      for (auto& f : futs) f.get();
                                                                  }

                                                                  for (size_t t = 0; t < num_threads; t++)
                                                                      to = std::min(to, nd[t]);

                                                                  id = C[to];
                                                                  int tot = 0;
                                                                  FlatCutIndex& uu = con_index.labels[u].cut_index;
                                                                  FlatCutIndex& vv = con_index.labels[v].cut_index;
                                                                  for (int i = 0; i < to; i++) {
                                                                      for (auto x : Q[C[i]].cut) {
                                                                          uint16_t cof = node_data[x].cut_offset;
                                                                          uint16_t cut_level = con_index.labels[x].cut_index.cut_level();
                                                                          if (before_w == std::abs((int)*(uu.distances() + get_offset(uu.dist_index(), cut_level) + cof) - (int)*(vv.distances() + get_offset(vv.dist_index(), cut_level) + cof)))
                                                                              upd[(++tot) % num_threads].push_back(x);
                                                                      }
                                                                  }
                                                              } else {
                                                                  size_t num_threads2 = 31;
                                                                  vector<NodeID> nd; nd.resize(num_threads2 + 1);
                                                                  size_t chunk_size = (CALL.size() + num_threads2 - 1) / num_threads2;

                                                                  // parallel get_cutdis2
                                                                  {
                                                                      auto& pool = global_pool();
                                                                      std::vector<std::future<void>> futs;
                                                                      futs.reserve(num_threads2);
                                                                      for (size_t t = 0; t < num_threads2; ++t) {
                                                                          size_t start = t * chunk_size;
                                                                          size_t end = std::min(start + chunk_size, CALL.size());
                                                                          futs.emplace_back(
                                                                              pool.submit(&Graph::get_cutdis2, this,
                                                                                          start, end, std::ref(C), std::ref(CALL), std::ref(con_index),
                                                                                          u, v, w, std::ref(nd[t]))
                                                                          );
                                                                      }
                                                                      for (auto& f : futs) f.get();
                                                                  }

                                                                  for (size_t t = 0; t < num_threads2; t++)
                                                                      to = std::min(to, nd[t]);

                                                                  id = C[to];
                                                                  int tott = 0;
                                                                  FlatCutIndex& uu = con_index.labels[u].cut_index;
                                                                  FlatCutIndex& vv = con_index.labels[v].cut_index;
                                                                  for (int i = 0; i < to; i++) {
                                                                      for (auto x : Q[C[i]].cut) {
                                                                          uint16_t cof = node_data[x].cut_offset;
                                                                          uint16_t cut_level = con_index.labels[x].cut_index.cut_level();
                                                                          if (w < std::abs((int)*(uu.distances() + get_offset(uu.dist_index(), cut_level) + cof) - (int)*(vv.distances() + get_offset(vv.dist_index(), cut_level) + cof)))
                                                                              upd[(++tott) % num_threads2].push_back(x);
                                                                      }
                                                                  }
                                                              }

                                                              // process_nodes over the id segment (parallel)
                                                              if (Q[id].ed - Q[id].st + 1 >= 100) {
                                                                  size_t chunk_size = (Q[id].ed - Q[id].st + 1 + num_threads - 1) / num_threads;
                                                                  auto& pool = global_pool();
                                                                  std::vector<std::future<void>> futs;
                                                                  futs.reserve(num_threads);
                                                                  for (size_t t = 0; t < num_threads; ++t) {
                                                                      size_t start = t * chunk_size + Q[id].st;
                                                                      size_t end = std::min(start + chunk_size, Q[id].ed + 1);
                                                                      futs.emplace_back(
                                                                          pool.submit(&Graph::process_nodes, this, start, end, std::ref(ct), std::ref(con_index),
                                                                                      u, v, w, before_w, std::ref(ci), std::ref(get_dis), md, std::ref(upd[t]),
                                                                                      std::ref(change[t]), etime)
                                                                      );
                                                                  }
                                                                  for (auto& f : futs) f.get();
                                                              } else {
                                                                  size_t local_threads = 4;
                                                                  size_t chunk_size = (Q[id].ed - Q[id].st + 1 + local_threads - 1) / local_threads;
                                                                  auto& pool = global_pool();
                                                                  std::vector<std::future<void>> futs;
                                                                  futs.reserve(local_threads);
                                                                  for (size_t t = 0; t < local_threads; ++t) {
                                                                      size_t start = t * chunk_size + Q[id].st;
                                                                      size_t end = std::min(start + chunk_size, Q[id].ed + 1);
                                                                      futs.emplace_back(
                                                                          pool.submit(&Graph::process_nodes, this, start, end, std::ref(ct), std::ref(con_index),
                                                                                      u, v, w, before_w, std::ref(ci), std::ref(get_dis), md, std::ref(upd[t]),
                                                                                      std::ref(change[t]), etime)
                                                                      );
                                                                  }
                                                                  for (auto& f : futs) f.get();
                                                              }
                                                              // change_remain_one / change_remain_one2 (parallel)
                                                                                                                          f_md.get();
                                                              if (w > before_w) {
                                                                  int fsz = 0;
                                                                  for (size_t t = 0; t < num_threads; ++t) fsz += (int)upd[t].size();
                                                                  vector<NodeID> cut_change;
                                                                  cut_change.reserve(fsz + 3);
                                                                  for (size_t t = 0; t < num_threads; ++t)
                                                                      cut_change.insert(cut_change.end(), upd[t].begin(), upd[t].end());

                                                                  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                                                                  std::shuffle(cut_change.begin(), cut_change.end(), std::default_random_engine(seed));

                                                                  size_t chunk_size = (cut_change.size() + num_threads - 1) / num_threads;
                                                                  auto& pool = global_pool();
                                                                  std::vector<std::future<void>> futs;
                                                                  futs.reserve(num_threads);
                                                                  for (size_t t = 0; t < num_threads; ++t) {
                                                                      size_t start = t * chunk_size;
                                                                      size_t end = std::min(start + chunk_size, cut_change.size());
                                                                      futs.emplace_back(
                                                                          pool.submit(&Graph::change_remain_one, this, start, end, std::ref(cut_change), std::ref(con_index),
                                                                                      w, before_w, u, v, std::ref(ci), std::ref(change[t]), std::ref(get_dis), md, (int)t, etime)
                                                                      );
                                                                  }
                                                                  for (auto& f : futs) f.get();
                                                              } else {
                                                                  int fsz = 0;
                                                                  for (size_t t = 0; t < num_threads; ++t) fsz += (int)upd[t].size();
                                                                  vector<NodeID> cut_change;
                                                                  cut_change.reserve(fsz + 3);
                                                                  for (size_t t = 0; t < num_threads; ++t)
                                                                      cut_change.insert(cut_change.end(), upd[t].begin(), upd[t].end());

                                                                  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                                                                  std::shuffle(cut_change.begin(), cut_change.end(), std::default_random_engine(seed));

                                                                  size_t chunk_size = (cut_change.size() + num_threads - 1) / num_threads;
                                                                  auto& pool = global_pool();
                                                                  std::vector<std::future<void>> futs;
                                                                  futs.reserve(num_threads);
                                                                  for (size_t t = 0; t < num_threads; ++t) {
                                                                      size_t start = t * chunk_size;
                                                                      size_t end = std::min(start + chunk_size, cut_change.size());
                                                                      futs.emplace_back(
                                                                          pool.submit(&Graph::change_remain_one2, this, start, end, std::ref(cut_change), std::ref(con_index),
                                                                                      w, before_w, u, v, std::ref(ci), std::ref(change[t]), std::ref(get_dis), md, (int)t, etime)
                                                                      );
                                                                  }
                                                                  for (auto& f : futs) f.get();
                                                              }

                                                              // keep original semantics

                                                              return;
                                                              remove_edge(u, v);
                                                              add_edge(u, v, before_w, true, true);
                                                              // (The code below remains unreachable as in your original function.)
                                                              for (size_t i = 0; i < num_threads; i++) {
                                                                  for (auto p : change[i]) {
                                                                      NodeID cut = p.first.first;
                                                                      NodeID node = p.first.second;
                                                                      if (cut == node) continue;
                                                                      if (con_index.is_contracted(node)) continue;
                                                                      distance_t ww = p.second;
                                                                      if (node_data[cut].cut_offset >= 10000) continue;
                                                                      FlatCutIndex x = con_index.labels[node].cut_index;
                                                                      FlatCutIndex cu = con_index.labels[cut].cut_index;
                                                                      uint16_t x_offset = get_offset(x.dist_index(), cu.cut_level());
                                                                      if (x.dist_index()[cu.cut_level()] <= x_offset + node_data[cut].cut_offset) continue;
                                                                      distance_t* x_ptr = x.distances() + x_offset + node_data[cut].cut_offset;
                                                                      *x_ptr = ww;
                                                                  }
                                                              }
                                                          }

                                                          void Graph::Contract(ContractionIndex& con_index, NodeID s, distance_t d) {

                                                              priority_queue<SearchNode> q;

                                                              q.push(SearchNode(d, s));
                                                              while (!q.empty()) {
                                                                  SearchNode next = q.top();
                                                                  q.pop();
                                                                  con_index.labels[next.node].distance_offset = next.distance;
                                                                  for (Neighbor n : node_data[next.node].neighbors) {
                                                                      if (con_index.labels[n.node].parent == next.node)
                                                                          q.push(SearchNode(next.distance + n.distance, n.node));
                                                                  }
                                                              }
                                                          }


                                                          void Graph::update_many_edge_up(const vector<CutIndex>& ci, ContractionIndex& con_index, vector<pair<pair<NodeID,NodeID>,pair<distance_t,distance_t>>> &edge) {

                                                              vector<NodeID>upd, edgid,ans;
                                                              edge_fg.resize(edge.size()+5);
                                                              for (size_t i = 0; i < edge.size();i++) {
                                                                  auto e = edge[i];
                                                                  distance_t w = e.second.second;
                                                                  distance_t before_w, tmp = infinity;
                                                                  NodeID u = e.first.first, v = e.first.second;
                                                                  before_w = e.second.first;
                                                                  if (w == before_w)continue;
                                                                  remove_edge(u, v);
                                                                  add_edge(u, v, 200000000+i, true, true);
                                                              }
                                                              unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                                                              std::shuffle(nodes.begin(), nodes.end(), std::default_random_engine(seed));
                                                              {
                                                                  size_t num_threads = 32;
                                                                  std::vector<std::thread> threads;
                                                                  std::vector<std::vector<NodeID>> edgg(num_threads);
                                                                  size_t chunk_size = (edge.size() + num_threads - 1) / num_threads;
                                                                  for (size_t t = 0; t < num_threads; ++t) {
                                                                      size_t start = t * chunk_size;
                                                                      size_t end = std::min(start + chunk_size, edge.size());
                                                                      for (size_t j = start; j < end; j++)
                                                                      {
                                                                          edgg[t].pb(j);
                                                                          //if(j<=50)cout << j << " ";
                                                                      }
                                                                  }

                                                                  for (size_t t = 0; t < num_threads; ++t)
                                                                      threads.emplace_back(&Graph::get_cut, this, Q[0].lchild, 0,std::cref(con_index),
                                                                                           std::cref(upd), std::cref(edgg[t]), std::cref(edge), std::ref(edge_fg));
                                                                      //get_cut(Q[0].lchild, 0, con_index, upd, edgid, edge, ans);
                                                                      for (auto& t : threads) {
                                                                          t.join(); // Wait for all threads to complete
                                                                      }
                                                                      //threads.clear();

                                                              }
                                                              //for(int i=0;i<100;i++)cout << edge_fg[i] << endl;

                                                              //�Ȳ����ж϶���һ����㣬��Щ�߻�Ӱ���������ľ��롣|Nodes|*|e|*|querytimes|/����Ч�ʡ�Ȼ��ֱ�ӿ�ʼdijkstra
                                                              //util::start_timer();
                                                              //return;
                                                              size_t num_threads = 32;
                                                              std::vector<std::thread> threads;
                                                              std::vector<std::vector<NodeID>> ndcg(num_threads);
                                                              size_t chunk_size = (nodes.size() + num_threads - 1) / num_threads;
                                                              for (size_t t = 0; t < num_threads; ++t) {
                                                                  size_t start = t * chunk_size;
                                                                  size_t end = std::min(start + chunk_size, nodes.size());
                                                                  threads.emplace_back(&Graph::update_up_par,this, start, end, std::cref(nodes), std::ref(con_index),
                                                                                       std::cref(edge), std::cref(ci),std::ref(ndcg[t]));
                                                              }

                                                              for (auto& th : threads) {
                                                                  th.join();
                                                              }
                                                              //return;
                                                              //cout << util::stop_timer() << endl;
                                                              threads.clear();
                                                              vector<NodeID>change_cut;

                                                              for (size_t i = 0; i < num_threads; i++) {
                                                                  for (NodeID x : ndcg[i]) {
                                                                      change_cut.pb(x);
                                                                  }
                                                              }
                                                              seed = std::chrono::system_clock::now().time_since_epoch().count();
                                                              std::shuffle(change_cut.begin(), change_cut.end(), std::default_random_engine(seed));
                                                              //return;

                                                              cout << change_cut.size() <<"!!"<< endl;//return;
                                                              //cout << change_cut.size() << endl;
                                                              util::start_timer();
                                                              chunk_size = (change_cut.size() + num_threads - 1) / num_threads;
                                                              for (size_t t = 0; t < num_threads; ++t) {
                                                                  size_t start = t * chunk_size;
                                                                  size_t end = std::min(start + chunk_size, change_cut.size());
                                                                  threads.emplace_back(&Graph::change_remain, this, start, end, std::ref(change_cut), std::ref(con_index),
                                                                                       std::ref(edge), std::ref(ci),t);
                                                              }
                                                              for (auto& th : threads) {
                                                                  th.join();
                                                              }
                                                              cout << util::stop_timer() << endl;

                                                              //con_index.get_distance(node, node);
                                                              /*
                                                               *   for (size_t i = 0; i < num_threads; i++) {
                                                               *       for (auto p : upd[i])
                                                               *       {
                                                               *           NodeID cut = p.first.first;
                                                               *           NodeID node = p.first.second;
                                                               *           if (cut == node)continue;
                                                               *           //
                                                               *           if (con_index.is_contracted(node))continue;
                                                               *           //if (con_index.isAncestor(cut, node))continue;
                                                               *           distance_t ww = p.second;
                                                               *           if (node_data[cut].cut_offset >= 10000)continue;
                                                               *           FlatCutIndex x = con_index.labels[node].cut_index;
                                                               *           FlatCutIndex c = con_index.labels[cut].cut_index;
                                                               *           uint16_t x_offset = get_offset(x.dist_index(), c.cut_level());
                                                               *           distance_t* x_ptr = x.distances() + x_offset+ node_data[cut].cut_offset;
                                                               *           distance_t* x_end = x_ptr + x.dist_index()[ci[cut].cut_level] - x_offset;
                                                               *           if (x.dist_index()[ci[cut].cut_level]  <= x_offset + node_data[cut].cut_offset) {
                                                               *               continue;
                                                               *               cout << x.dist_index()[ci[cut].cut_level] << " " << x_offset + node_data[cut].cut_offset << endl;
                                                               *
                                                               *               cout << cut << "% " << node << " $" << ww <<" #"<< node_data[cut].cut_offset<<" ?"<<c.cut_level()<<" !"<<x_offset<<" $"<< x.dist_index()[ci[cut].cut_level]<<
                                                               *                   "  "<< x_offset + node_data[cut].cut_offset <<endl;
                                                               *               exit(0);
                                                          }


                                                          /*while (x_ptr != x_end)
                                                           *           {
                                                           *               x_ptr++;
                                                          }
                                                          //con_index.get_distance(node,node);
                                                          if (x_ptr < x.distances() || x_ptr >= x_end) {
                                                              continue;
                                                              std::cerr << "Pointer out of bounds: x_ptr is invalid!" << std::endl;
                                                              // ��ӡ������Ϣ��������λ����
                                                              std::cerr << "x_ptr: " << x_ptr << ", x.distances(): " << x.distances() << ", x_end: " << x_end << std::endl;
                                                          }
                                                          else {//
                                                              // ��ȫ����
                                                              *x_ptr = ww;
                                                          }
                                                          }

                                                          }

                                                          */


                                                          for (size_t i = 0; i < edge.size(); i++) {
                                                              auto e = edge[i];
                                                              distance_t w = e.second.second;
                                                              distance_t before_w, tmp = infinity;
                                                              NodeID u = e.first.first, v = e.first.second;
                                                              before_w = e.second.first;
                                                              if (w == before_w)continue;
                                                              remove_edge(u, v);
                                                              add_edge(u, v, w, true, true);
                                                          }
                                                          }




                                                          void get_edge(std::istream& in,
                                                                        std::vector<std::pair<std::pair<NodeID, NodeID>, std::pair<distance_t, distance_t>>>& edge,
                                                                        size_t fg,
                                                                        const ContractionIndex& con_index)
                                                          {
                                                              (void)in;
                                                              (void)con_index;

                                                              edge.clear();

                                                              std::ifstream fin("update.txt");
                                                              if (!fin) {
                                                                  std::cerr << "[get_edge] ERROR: cannot open update.txt\n";
                                                                  return;
                                                              }

                                                              const bool read_all = (fg == 0);
                                                              size_t count = 0;

                                                              std::string line;
                                                              NodeID u, v;
                                                              long long ow, nw;

                                                              while (std::getline(fin, line)) {
                                                                  if (line.empty() || line[0] == '#') continue;
                                                                  std::istringstream ss(line);
                                                                  if (!(ss >> u >> v >> ow >>nw)) continue;
                                                                  //cout << u << " " << v << " " << ow << endl;
                                                                  edge.push_back({ {u, v}, { static_cast<distance_t>(ow), static_cast<distance_t>(1.5*ow) } });
                                                                  ++count;
                                                                  if (!read_all && count >= fg) break;
                                                              }
                                                          }




                                                          int main(int argc, char *argv[])
                                                          {
                                                              \
                                                              if (argc < 2)
                                                              {
                                                                  cout << "syntax: " << argv[0] << " [balance] <filename> ... <filename>" << endl;
                                                                  return 0;
                                                              }
                                                              // check for balance parameter
                                                              double balance = atof(argv[1]);
                                                              int file_start = 2;
                                                              if (balance == 0.0)
                                                              {
                                                                  balance = 0.15;
                                                                  file_start = 1;
                                                              }

                                                              #ifdef NO_SHORTCUTS
                                                              cout << "shortcuts disabled" << endl;
                                                              #elif defined(ALL_SHORTCUTS)
                                                              cout << "redundant shortcuts enabled" << endl;
                                                              #else
                                                              cout << "shortcuts enabled" << endl;
                                                              #endif

                                                              #ifdef PRUNING
                                                              cout << "pruning enabled" << endl;
                                                              #else
                                                              cout << "pruning disabled" << endl;
                                                              #endif

                                                              #ifdef CONTRACT2D
                                                              cout << "path contraction enabled" << endl;
                                                              #else
                                                              cout << "path contraction disabled" << endl;
                                                              #endif

                                                              #ifdef MULTI_THREAD
                                                              cout << "multi-threading enabled" << endl;
                                                              cout << "threads supported by hardware: " << thread::hardware_concurrency() << endl;
                                                              #else
                                                              cout << "multi-threading disabled" << endl;
                                                              #endif

                                                              #ifdef NDEBUG
                                                              srand(time(nullptr));
                                                              #endif
                                                              for (int f = file_start; f < argc; f++)
                                                              {
                                                                  const char* filename = argv[f];
                                                                  vector<ResultData> results;
                                                                  bool use_buckets;
                                                                  for (size_t i = 0; i < repeats; i++)
                                                                  {
                                                                      cout << endl << "reading graph from " << filename << endl;
                                                                      fstream fs(filename);
                                                                      Graph g;
                                                                      read_graph(g, fs);
                                                                      fs.close();
                                                                      cout << "read " << g.node_count() << " vertices and " << g.edge_count() << " edges" << flush;
                                                                      distance_t diameter = g.diameter(true);
                                                                      cout << " (diameter=" << g.diameter(false) << "|" << diameter << ")" << endl;
                                                                      // check for redundant edges
                                                                      vector<Edge> redundant_edges;
                                                                      util::start_timer();
                                                                      g.get_redundant_edges(redundant_edges);
                                                                      #ifdef REMOVE_REDUNDANT
                                                                      for (Edge e : redundant_edges)
                                                                          g.remove_edge(e.a, e.b);
                                                                      cout << "removed " << redundant_edges.size() << " redundant edges in " << util::stop_timer() << "s" << endl;
                                                                      #else
                                                                      cout << "found " << redundant_edges.size() << " redundant edges in " << util::stop_timer() << "s" << endl;
                                                                      #endif
                                                                      #ifdef CONTRACT
                                                                      util::start_timer();
                                                                      size_t old_size = g.node_count();
                                                                      vector<Neighbor> closest;
                                                                      g.contract(closest);
                                                                      cout << "contracted to " << g.node_count() << " vertices (" << g.node_count() * 100 / max<size_t>(1, old_size) << "%) and "
                                                                      << g.edge_count() << " edges in " << util::stop_timer() << "s" << endl;
                                                                      #endif
                                                                      #ifdef CONTRACT2D
                                                                      size_t deg2nodes = 0;
                                                                      for (NodeID node : g.get_nodes())
                                                                          if (g.degree(node) == 2)
                                                                              deg2nodes++;
                                                                      cout << deg2nodes << " of these vertices (" << deg2nodes * 100 / max<size_t>(1, g.node_count()) << "%) have degree 2" << endl;
                                                                      #endif
                                                                      #ifdef NDEBUG
                                                                      g.randomize();
                                                                      #endif
                                                                      ResultData result = {};
                                                                      // construct index
                                                                      Graph::show_progress(true);
                                                                      vector<CutIndex> ci;
                                                                      util::start_timer();

                                                                      size_t shortcuts = g.create_cut_index(ci, balance);
                                                                      #ifdef PRUNING
                                                                      result.pruning_2hop = get_2hop_pruning(ci);
                                                                      result.pruning_3hop = get_3hop_pruning(ci);
                                                                      result.pruning_tail = get_tail_pruning(ci);
                                                                      #endif
                                                                      #ifdef CONTRACT
                                                                      ContractionIndex con_index(ci, closest);
                                                                      #else
                                                                      ContractionIndex con_index(ci);
                                                                      #endif
                                                                      result.index_time = util::stop_timer();
                                                                      result.index_size = con_index.size() / MB;
                                                                      result.label_count = con_index.label_count();
                                                                      result.max_label_count = con_index.max_label_count();
                                                                      result.index_height = con_index.height();
                                                                      result.avg_cut_size = con_index.avg_cut_size();
                                                                      result.max_cut_size = con_index.max_cut_size();
                                                                      cout << "created index of size " << result.index_size << " MB in " << result.index_time << "s using " << shortcuts << " shortcuts" << endl;
                                                                      cout << "partition tree contains " << con_index.non_empty_cuts() << " non-empty cuts (" << 100 * con_index.non_empty_cuts() / con_index.uncontracted_count() << "% of uncontracted vertices)" << endl;
                                                                      #ifdef PRUNING
                                                                      size_t unpruned_labels = max<size_t>(1, result.label_count + result.pruning_tail);
                                                                      cout << "3-HOP pruning could remove " << result.pruning_3hop << " labels (" << result.pruning_3hop * 100 / unpruned_labels << "%)" << endl;
                                                                      cout << "2-HOP pruning could remove " << result.pruning_2hop << " labels (" << result.pruning_2hop * 100 / unpruned_labels << "%)" << endl;
                                                                      cout << "tail pruning *has* removed " << result.pruning_tail << " labels (" << result.pruning_tail * 100 / unpruned_labels << "%)" << endl;
                                                                      #endif
                                                                      g.reset(); // needed for distance testing
                                                                      g.init(con_index);
                                                                      // show memory consumption
                                                                      rusage usage;
                                                                      if (getrusage(RUSAGE_SELF, &usage) != -1)
                                                                          cout << "maximum memory used: " << usage.ru_maxrss / 1024 << " MB" << endl;
                                                                      //cout << "sd2" << "??"<<endl;
                                                                      vector<pair<NodeID,NodeID>> queries;
                                                                      for (size_t i = 0; i < nr_queries; i++)
                                                                      {
                                                                          pair<NodeID, NodeID> p = g.random_pair();
                                                                          queries.push_back(p);
                                                                      }
                                                                      if (queries.size() > nr_query_tests)
                                                                          queries.resize(nr_query_tests);
                                                                      util::make_set(queries);

                                                                      // test query speed by distance, as for H2H / P2H
                                                                      use_buckets = diameter >= bucket_min * nr_buckets;
                                                                      if (use_buckets)
                                                                      {
                                                                          cout << "generating queries by distance: " << flush;
                                                                          vector<vector<pair<NodeID,NodeID>>> query_buckets(nr_buckets);
                                                                          util::start_timer();
                                                                          g.random_pairs(query_buckets, bucket_min, bucket_size, con_index);
                                                                          cout << " in " << util::stop_timer() << "s" << endl;
                                                                          for (size_t bucket = 0; bucket < query_buckets.size(); bucket++)
                                                                          {
                                                                              util::start_timer();
                                                                              for (pair<NodeID,NodeID> q : query_buckets[bucket])
                                                                                  con_index.get_distance(q.first, q.second);
                                                                              result.bucket_query_times.push_back(util::stop_timer());
                                                                              result.bucket_hoplinks.push_back(con_index.avg_hoplinks(query_buckets[bucket]));
                                                                              int tot = 0,p=0;
                                                                              for (pair<NodeID, NodeID>q : query_buckets[bucket]) {
                                                                                  FlatCutIndex a, b;
                                                                                  a = con_index.labels[q.first].cut_index, b = con_index.labels[q.second].cut_index;
                                                                                  size_t cut_level = PBV::lca_level(*a.partition_bitvector(), *b.partition_bitvector());
                                                                                  uint16_t a_offset = get_offset(a.dist_index(), cut_level);
                                                                                  uint16_t b_offset = get_offset(b.dist_index(), cut_level);
                                                                                  p += min(a.dist_index()[cut_level] - a_offset, b.dist_index()[cut_level] - b_offset);
                                                                                  tot++;
                                                                              }
                                                                              cout << "ran " << query_buckets[bucket].size() << " queries (bucket " << bucket << ") in " << result.bucket_query_times.back() << "s (hoplinks=" << result.bucket_hoplinks.back() << ")" <<" "<<p/tot<< endl;
                                                                          }
                                                                      }
                                                                      results.push_back(result);

                                                                      vector < pair<pair<NodeID, NodeID>, pair<distance_t, distance_t>>>edge,edge2;
                                                                      fstream ft(filename);
                                                                      get_edge(ft, edge,100,con_index);//1������up��-1������down

                                                                      cout << edge.size()<< endl;
                                                                      ft.close();
                                                                      srand(time(0));

                                                                      long long p = 0;
                                                                      int tot = 0;

                                                                      util::start_timer();
                                                                      int e=edge.size();
                                                                      for (int i = 1; i <= e; i++) {
                                                                          int x = i-1;
                                                                          auto e = edge[x];
                                                                          NodeID u = e.first.first, v = e.first.second;
                                                                          distance_t w = e.second.second;
                                                                          distance_t before_w = e.second.first;
                                                                          //cout << w << " "<<before_w<<endl;
                                                                          //exit(0);

                                                                          ContractionLabel a = con_index.labels[u], b = con_index.labels[v];
                                                                          if (con_index.is_contracted(u) || con_index.is_contracted(v)) {
                                                                              g.remove_edge(u, v);
                                                                              g.add_edge(u, v, w, true);
                                                                              if (a.distance_offset > b.distance_offset)
                                                                                  g.Contract(con_index, u, b.distance_offset + w);
                                                                              else if (a.distance_offset < b.distance_offset)
                                                                                  g.Contract(con_index, v, a.distance_offset +w);

                                                                          }
                                                                          else
                                                                          {
                                                                              //     continue;
                                                                              g.update_edge(ci, con_index, edge[x].first.first, edge[x].first.second, edge[x].second.second,i);
                                                                          }
                                                                      }
                                                                      cout << "test update time: "<<util::stop_timer() <<endl;
                                                                  }
                                                              }

                                                              return 0;
                                                          }
