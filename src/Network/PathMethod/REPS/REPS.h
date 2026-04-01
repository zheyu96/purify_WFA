#include "../PathMethodBase/PathMethod.h"

using namespace std;

class REPS : public PathMethod {
public:
    REPS();
    ~REPS();
    void build_paths(Graph graph, vector<SDpair> requests);
};
