#include "dbscan.h"

namespace ofec {
    int DBSCAN::run(){
        int clusterID = 1;
        //std::vector<Point>::iterator iter;
        for (size_t i = 0;i<m_points.size(); ++i) {
            if (m_points[i]->clusterID == UNCLASSIFIED) {
                if (expandCluster(*m_points[i], clusterID) != FAILURE) {
                    clusterID += 1;
                }
            }
        }
        return 0;
    }

    int DBSCAN::expandCluster(Point &point, int clusterID){
        std::vector<int> clusterSeeds = calculateCluster(point);
        if (clusterSeeds.size() < m_minPoints){
            point.clusterID = NOISE;
            return FAILURE;
        }
        else{
            int index = 0, indexCorePoint = 0;
            for (size_t j = 0;j< clusterSeeds.size();++j) {
                m_points[clusterSeeds[j]]->clusterID = clusterID;
                for (size_t i = 0; i < point.m_variables.size(); ++i) {//找到当前核心点，并去除
                    if (m_points[clusterSeeds[j]]->m_variables[i] != point.m_variables[i]) {
                        break;
                    }
                    if (i == point.m_variables.size() - 1) {
                        indexCorePoint = index;
                    }
                }
                ++index;
            }
            clusterSeeds.erase(clusterSeeds.begin() + indexCorePoint);

            for (std::vector<int>::size_type i = 0, n = clusterSeeds.size(); i < n; ++i){
                std::vector<int> clusterNeighors = calculateCluster(*m_points[clusterSeeds[i]]);

                if (clusterNeighors.size() >= m_minPoints){
                    for (size_t j = 0;j< clusterNeighors.size(); ++j) {
                        if (m_points[clusterNeighors[j]]->clusterID == UNCLASSIFIED || m_points[clusterNeighors[j]]->clusterID == NOISE){
                            if (m_points[clusterNeighors[j]]->clusterID == UNCLASSIFIED){
                                clusterSeeds.push_back(clusterNeighors[j]);
                                n = clusterSeeds.size();
                            }
                            m_points[clusterNeighors[j]]->clusterID = clusterID;
                        }
                    }
                }
            }
            return SUCCESS;
        }
    }

    std::vector<int> DBSCAN::calculateCluster(const Point &point){
        int index = 0;
        std::vector<int> clusterIndex;
        for (size_t i = 0;i<m_points.size();++i) {
            if (calculateDistance(point, *m_points[i]) <= m_epsilon){
                clusterIndex.push_back(index);
            }
            index++;
        }
        return clusterIndex;
    }

    inline Real DBSCAN::calculateDistance(const Point &pointCore, const Point &pointTarget){
        Real sum = 0.;
        for (size_t i = 0; i < pointCore.m_variables.size(); ++i) {
            sum += (std::pow(pointCore.m_variables[i] - pointTarget.m_variables[i],2));
        }
        return std::sqrt(sum);
    }
}