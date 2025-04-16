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

vector<int> hungarian_algorithm(const vector<vector<double>>& squareCost) {
    int size = squareCost.size();
    if (size == 0) return {};

    vector<vector<double>> matrix = squareCost;  
    vector<bool> rowCover(size, false), colCover(size, false);
    vector<vector<bool>> star(size, vector<bool>(size, false));
    vector<vector<bool>> prime(size, vector<bool>(size, false));
    vector<std::pair<int,int>> path(2*size);

    int step = 1;
    int pathRow = -1, pathCol = -1;

    auto step1 = [&]() {
        // 각 행에서 최소값을 빼 0 만듦
        for(int i=0; i<size; i++){
            double mn = std::numeric_limits<double>::infinity();
            for(int j=0; j<size; j++){
                if(matrix[i][j] < mn) mn = matrix[i][j];
            }
            for(int j=0; j<size; j++){
                matrix[i][j] -= mn;
            }
        }
        step = 2;
    };

    auto step2 = [&]() {
        // 0인 곳에 별표, 해당 행과 열을 커버
        for(int i=0; i<size; i++){
            for(int j=0; j<size; j++){
                if(std::fabs(matrix[i][j])<1e-12 && !rowCover[i] && !colCover[j]) {
                    star[i][j] = true;
                    rowCover[i] = true;
                    colCover[j] = true;
                }
            }
        }
        // 커버 풀기
        std::fill(rowCover.begin(), rowCover.end(), false);
        std::fill(colCover.begin(), colCover.end(), false);
        step = 3;
    };

    auto step3 = [&]() {
        // 별표가 표시된 열을 커버
        int colCount = 0;
        for(int j=0; j<size; j++){
            for(int i=0; i<size; i++){
                if(star[i][j]) {
                    colCover[j] = true;
                    break;
                }
            }
        }
        // 커버된 열 수 체크
        for(int j=0; j<size; j++){
            if(colCover[j]) colCount++;
        }
        if(colCount >= size) {
            step = 7; // 매칭 완료
        } else {
            step = 4;
        }
    };

    // step4~6에서 0을 찾아 프라임 붙이기, 교대 경로 찾기 등
    std::function<void()> step4 = [&]() {
        while(true){
            int row=-1, col=-1;
            // 덮이지 않은 행/열에서 0을 찾으면 prime 표시
            for(int i=0; i<size; i++){
                if(!rowCover[i]) {
                    for(int j=0; j<size; j++){
                        if(!colCover[j] && std::fabs(matrix[i][j])<1e-12) {
                            row = i; 
                            col = j; 
                            goto FOUND;
                        }
                    }
                }
            }
            FOUND:;
            if(row == -1){
                // 못 찾으면 step6
                step = 6;
                return;
            }
            prime[row][col] = true;
            // 같은 행에 star가 있는지?
            int starCol = -1;
            for(int j=0; j<size; j++){
                if(star[row][j]) {
                    starCol = j; 
                    break;
                }
            }
            if(starCol != -1){
                // 행 커버, star col 언커버
                rowCover[row] = true;
                colCover[starCol] = false;
            } else {
                // 교대 경로 만들 준비
                pathRow = row;
                pathCol = col;
                step = 5;
                return;
            }
        }
    };

    auto step5 = [&]() {
        // 교대 경로(path) 구성
        int count=0;
        path[count] = {pathRow, pathCol};
        bool done = false;
        while(!done){
            // 1) path[count].second 열에 star가 있는 행 r 찾기
            int r_=-1;
            for(int i=0; i<size; i++){
                if(star[i][ path[count].second ]) {
                    r_=i;
                    break;
                }
            }
            if(r_<0){ // 없으면 끝
                done = true;
                break;
            }
            count++;
            path[count] = {r_, path[count-1].second};
            // 2) 그 행에 prime이 있는 열 c 찾기
            int c_=-1;
            for(int j=0; j<size; j++){
                if(prime[r_][j]) {
                    c_=j;
                    break;
                }
            }
            count++;
            path[count] = {r_, c_};
        }
        // 교대 경로 뒤집기(star->off, prime->on)
        for(int i=0; i<=count; i++){
            auto [rr, cc] = path[i];
            if(star[rr][cc]) {
                star[rr][cc] = false;
            } else {
                star[rr][cc] = true;
                prime[rr][cc] = false;
            }
        }
        // prime 지우고 cover 해제
        for(int i=0; i<size; i++){
            std::fill(prime[i].begin(), prime[i].end(), false);
        }
        std::fill(rowCover.begin(), rowCover.end(), false);
        std::fill(colCover.begin(), colCover.end(), false);
        step=3;
    };

    auto step6 = [&]() {
        // 덮이지 않은 셀 중 최솟값 d 찾기
        double minv = std::numeric_limits<double>::infinity();
        for(int i=0; i<size; i++){
            if(!rowCover[i]) {
                for(int j=0; j<size; j++){
                    if(!colCover[j] && matrix[i][j]<minv) {
                        minv = matrix[i][j];
                    }
                }
            }
        }
        // rowCover=true인 행은 +d, colCover=false인 열은 -d
        for(int i=0; i<size; i++){
            for(int j=0; j<size; j++){
                if(rowCover[i]) matrix[i][j]+=minv;
                if(!colCover[j]) matrix[i][j]-=minv;
            }
        }
        step=4;
    };

    // 메인 루프
    bool done=false;
    while(!done){
        switch(step){
        case 1: step1(); break;
        case 2: step2(); break;
        case 3: step3(); break;
        case 4: step4(); break;
        case 5: step5(); break;
        case 6: step6(); break;
        case 7: done=true; break;
        }
    }

    // 결과 추출
    vector<int> assign(size, -1);
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(star[i][j]) {
                assign[i]=j;
                break;
            }
        }
    }
    return assign;
}

vector<int> hungarian_algorithm_with_padding(const vector<vector<double>>& costMatrix)
{
    int n = costMatrix.size();
    if(n==0) return {};
    int m = costMatrix[0].size();
    int size = std::max(n,m);

    // 1) 0-padding된 정사각행렬 만들기
    vector<vector<double>> squareCost(size, vector<double>(size, 0.0));
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            squareCost[i][j] = costMatrix[i][j];
        }
    }

    // 2) 정사각 헝가리안 알고리즘 실행
    vector<int> fullAssignment = hungarian_algorithm(squareCost);

    // 3) 실제 n개 행(에이전트)에 해당하는 결과만 추려서 반환
    //    col >= m이면 더미 열로 간주 -> -1
    vector<int> assignment(n, -1);
    for(int i=0; i<n; i++){
        int c = fullAssignment[i];
        if(c >= 0 && c < m) {
            assignment[i] = c;
        } else {
            assignment[i] = -1;
        }
    }
    return assignment;
}

CTAResult CTA(vector<Agent>& agents, vector<Task>& tasks)
{
    CTAResult result;
    result.assignmentResult.resize(agents.size());

    const double UAV_MAX_SPEED = 5.0;
    const double bat = 32767.0;

    while(!tasks.empty()) {
        int n = (int)agents.size();  
        int m = (int)tasks.size();   

        vector<vector<double>> rewardMatrix(n, vector<double>(m, 0.0));
        vector<vector<double>> disMatrix(n, vector<double>(m, 0.0));

        for(int i=0; i<n; i++){
            Agent& ag = agents[i];
            for(int j=0; j<m; j++){
                Task& tk = tasks[j];
                if(tk.waypoints.empty()) {
                    rewardMatrix[i][j] = 0.0;
                    continue;
                }
                // 거리 계산
                double distAtoFirst = sqrt(pow(get_distance_latlon(
                    ag.latitude, ag.longitude,
                    tk.waypoints[0].latitude, tk.waypoints[0].longitude
                ),2)+pow(tk.waypoints[0].altitude,2)-pow(ag.altitude,2));
                double distInsideTask = 0.0;
                for(size_t w=0; w+1<tk.waypoints.size(); w++){
                    distInsideTask += sqrt(pow(get_distance_latlon(
                        tk.waypoints[w].latitude, tk.waypoints[w].longitude,
                        tk.waypoints[w+1].latitude, tk.waypoints[w+1].longitude
                    ),2)+pow(tk.waypoints[w+1].altitude,2)-pow(tk.waypoints[w].altitude,2));
                }
    
                double totalDist = distAtoFirst + distInsideTask;
                disMatrix[i][j] = totalDist;
                double t = totalDist / UAV_MAX_SPEED;

                double reward = 100.0 * exp(-0.001*t);
                     //           / (1.0 + exp(-10.0*(bat - 32767.0)));
                rewardMatrix[i][j] = reward;
            }
        }
        result.rewardMatrices.push_back(rewardMatrix);
        result.disMatrices.push_back(disMatrix);
        // 2) costMatrix(n x m)
        vector<vector<double>> costMatrix(n, vector<double>(m, 0.0));
        for(int i=0; i<n; i++){
            for(int j=0; j<m; j++){
                costMatrix[i][j] = -rewardMatrix[i][j]; 
            }
        }

        // 3) 헝가리안(0-padding 버전) 실행
        vector<int> assignment = hungarian_algorithm_with_padding(costMatrix);

        // 4) 결과 반영
        vector<int> tasksToRemove;
        tasksToRemove.reserve(n);
        for(int i=0; i<n; i++){
            int col = assignment[i];
            if(col>=0 && col<m) {
                result.assignmentResult[i].push_back(tasks[col].name);

                // 에이전트 위치 = 태스크 마지막 waypoint
                if(!tasks[col].waypoints.empty()){
                    auto& lastWp = tasks[col].waypoints.back();
                    agents[i].latitude  = lastWp.latitude;
                    agents[i].longitude = lastWp.longitude;
                    agents[i].altitude  = lastWp.altitude;
                }
                // 할당된 태스크는 목록에서 제거
                tasksToRemove.push_back(col);
            }
        }

        if(!tasksToRemove.empty()){
            sort(tasksToRemove.begin(), tasksToRemove.end(), std::greater<int>());
            for(int idx : tasksToRemove){
                if(idx>=0 && idx<(int)tasks.size()) {
                    tasks.erase(tasks.begin()+idx);
                }
            }
        } else {
            // 더 이상 할당할 태스크 없음 -> break
            break;
        }
    }

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

    // 4) 결과 출력

    // (a) 각 반복에서의 리워드 매트릭스
    cout << "==== Reward Matrices per Iteration ====\n";
    for(size_t k = 0; k < result.rewardMatrices.size(); ++k) {
        const auto& matrix = result.rewardMatrices[k];
        cout << "[Iteration " << k << "] Reward Matrix:\n";
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

    // (b) 최종 할당 결과
    cout << "==== Final Assignment Result ====\n";
    for(size_t i = 0; i < result.assignmentResult.size(); ++i) {
        cout << "[Agent " << (i+1) << " - " << agents[i].name << "] Assigned Tasks: ";
        for(const auto& tname : result.assignmentResult[i]) {
            cout << tname << " ";
        }
        cout << endl;
    }

    cout << "==== Distance Matrix ====\n";
    for(size_t k = 0; k < result.disMatrices.size(); ++k) {
        const auto& matrix = result.disMatrices[k];
        cout << " Distance Matrix:\n";
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
