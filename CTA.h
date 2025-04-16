#ifndef CTA_H
#define CTA_H

#include <vector>
#include <string>

using namespace std;

struct Agent {
    string name;
    float latitude;
    float longitude;
    float altitude;
};

struct Waypoint {
    float latitude;
    float longitude;
    float altitude;
};

struct Task {
    string name;
    vector<Waypoint> waypoints;
};

struct CTAResult {
    vector<vector<vector<double>>> rewardMatrices;
};

CTAResult CTA(vector<Agent>& agents, vector<Task>& tasks);
double get_distance_latlon(double current_lat, double current_lon, double target_lat, double target_lon);

#endif