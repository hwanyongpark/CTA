#include "CTA.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#define CONSTANTS_RADIUS_OF_EARTH 6371000.0

using namespace std;

double get_distance_latlon(double current_lat, double current_lon, double target_lat, double target_lon){
    const double current_lat_rad = current_lat / 180.0 * M_PI;
    const double target_lat_rad = target_lat / 180.0 * M_PI;

    const double d_lat = target_lat_rad - current_lat_rad;
    const double d_lon = (target_lon - current_lon) / 180.0 * M_PI;

    const double a = sin(d_lat / 2.0) * sin(d_lat / 2.0) + sin(d_lon / 2.0) * sin(d_lon / 2.0) * cos(current_lat_rad) * cos(target_lat_rad);
    const double c = atan2(sqrt(a), sqrt(1.0 - a));

    return static_cast<double>(CONSTANTS_RADIUS_OF_EARTH * 2.0 * c);
}

CTAResult CTA(vector<Agent>& agents, vector<Task>& tasks)
{
    CTAResult result;

    const double UAV_MAX_SPEED = 5.0;
    const double bat = 32767.0;

    int n = (int)agents.size();  
    int m = (int)tasks.size();   

    vector<vector<double>> rewardMatrix(n, vector<double>(m, 0.0));

    for(int i=0; i<n; i++){
        Agent& ag = agents[i];
        for(int j=0; j<m; j++){
            Task& tk = tasks[j];
            if(tk.waypoints.empty()) {
                rewardMatrix[i][j] = 0.0;
                continue;
            }
            // 거리 계산
            double distAtoFirst = get_distance_latlon(
                ag.latitude, ag.longitude,
                tk.waypoints[0].latitude, tk.waypoints[0].longitude
            );
            double distInsideTask = 0.0;
            for(size_t w=0; w+1<tk.waypoints.size(); w++){
                distInsideTask += get_distance_latlon(
                    tk.waypoints[w].latitude, tk.waypoints[w].longitude,
                    tk.waypoints[w+1].latitude, tk.waypoints[w+1].longitude
                );
            }
            double totalDist = distAtoFirst + distInsideTask;
            double t = totalDist / UAV_MAX_SPEED;

            double reward = 100.0 * exp(-0.001*t);
                    //           / (1.0 + exp(-10.0*(bat - 32767.0)));
            rewardMatrix[i][j] = reward;
        }
    }
    result.rewardMatrices.push_back(rewardMatrix);

    return result;
}

int main() {
    //Agent 정의
    vector<Agent> agents = {
        {"Agent1", 47.397751f, 8.545608f, 0.0f},
        {"Agent2", 47.397751f, 8.545661f, 0.0f},
        {"Agent3", 47.397751f, 8.545728f, 0.0f}
    };

    //Task 정의 (각 Task에 여러 waypoint 가능)
    vector<Task> tasks;

    {
        Task t;
        t.name = "Task0";
        t.waypoints.push_back({47.3978589831117f, 8.54573672902687f, 10.0f});
        t.waypoints.push_back({47.3979160168883f, 8.54573672902687f, 10.0f});
        t.waypoints.push_back({47.3979160168883f, 8.54566327097313f, 10.0f});
        t.waypoints.push_back({47.3978589831117f, 8.54566327097313f, 10.0f});
        t.waypoints.push_back({47.3978589831117f, 8.54573672902687f, 10.0f});
        tasks.push_back(t);
    }
    {
        Task t;
        t.name = "Task1";
        t.waypoints.push_back({47.3979339831118f, 8.54573672902687f, 10.0f});
        t.waypoints.push_back({47.3979910168883f, 8.54573672902687f, 10.0f});
        t.waypoints.push_back({47.3979339831118f, 8.54566327097313f, 10.0f});
        t.waypoints.push_back({47.3979910168883f, 8.54566327097313f, 10.0f});
        t.waypoints.push_back({47.3979339831118f, 8.54573672902687f, 10.0f});
        tasks.push_back(t);
    }
    {
        Task t;
        t.name = "Task2";
        t.waypoints.push_back({47.3978589831117f, 8.54583672902687, 10.0f});
        t.waypoints.push_back({47.3979160168883f, 8.54583672902687f, 10.0f});
        t.waypoints.push_back({47.3978589831117f, 8.54576327097313f, 10.0f});
        t.waypoints.push_back({47.3979160168883f, 8.54576327097313f, 10.0f});
        t.waypoints.push_back({47.3978589831117f, 8.54583672902687f, 10.0f});
        tasks.push_back(t);
    }
    {
        Task t;
        t.name = "Task3";
        t.waypoints.push_back({47.3979339831118f, 8.54583672902687f, 10.0f});
        t.waypoints.push_back({47.3979910168883f, 8.54583672902687f, 10.0f});
        t.waypoints.push_back({47.3979910168883f, 8.54576327097313f, 10.0f});
        t.waypoints.push_back({47.3979339831118f, 8.54576327097313f, 10.0f});
        t.waypoints.push_back({47.3979339831118f, 8.54583672902687f, 10.0f});
        tasks.push_back(t);
    }
    CTAResult result = CTA(agents, tasks);

    cout << "==== Reward Matrices ====\n";
    for(size_t k = 0; k < result.rewardMatrices.size(); ++k) {
        const auto& matrix = result.rewardMatrices[k];
        cout << " Reward Matrix:\n";
        if(matrix.empty()) {
            cout << "  (empty)\n";
            continue;
        }
        // 행 × 열
        for(size_t i = 0; i < matrix.size(); ++i) {
            for(size_t j = 0; j < matrix[i].size(); ++j) {
                cout << matrix[i][j] << " ";
            }
            cout << "\n";
        }
        cout << endl;
    }

    return 0;
}
