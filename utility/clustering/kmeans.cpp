#include "kmeans.h"
#include "../../core/problem/solution.h"
#include "../functional.h"
#include "../../core/problem/continuous/continuous.h"

namespace ofec {
	void KMeansClu::setData(const std::vector<const SolutionBase*> sols, Problem *pro) {
		if (dynamic_cast<Continuous *>(pro) == nullptr)
			throw MyExcept("K-means clustering is only feasible in continuous search space");
		m_samples.resize(sols.size());
		for (size_t i = 0; i < sols.size(); ++i)
			m_samples[i] = dynamic_cast<const Solution<>*>(sols[i])->variable().vect();
	}

	void KMeansClu::clustering(size_t K) {

	}

	void KMeansClu::clustering(const std::list<const SolutionBase*> &init_central_sols) {
		std::vector<std::vector<Real>> centers;
		for (auto sol : init_central_sols)
			centers.push_back(dynamic_cast<const Solution<>*>(sol)->variable().vect());
		updateClusters(centers);
		updateCenters(centers);
		while (true) {
			Real old_sum_inner_distance = m_sum_inner_distance;
			updateClusters(centers);
			if (abs(old_sum_inner_distance - m_sum_inner_distance) / old_sum_inner_distance < 0.01)
				break;
			updateCenters(centers);
		}
		updateClusterCenters(centers);
	}

	void KMeansClu::updateClusters(const std::vector<std::vector<Real>> &centers) {
		m_clusters.clear();
		m_clusters.resize(centers.size());
		m_sum_inner_distance = 0;
		for (size_t i = 0; i < m_samples.size(); ++i) {
			size_t id_nearest_cluster = 0;
			Real min_distance = euclideanDistance(m_samples[i].begin(), m_samples[i].end(), centers[0].begin());
			for (size_t k = 1; k < centers.size(); ++k) {
				Real distance = euclideanDistance(m_samples[i].begin(), m_samples[i].end(), centers[k].begin());
				if (distance < min_distance) {
					min_distance = distance;
					id_nearest_cluster = k;
				}
			}
			m_clusters[id_nearest_cluster].push_back(i);
			m_sum_inner_distance += min_distance;
		}
	}

	void KMeansClu::updateCenters(std::vector<std::vector<Real>> &centers) {
		for (size_t k = 0; k < centers.size(); ++k) {
			for (size_t j = 0; j < centers[k].size(); ++j)
				centers[k][j] = 0;
			for (size_t i : m_clusters[k]) {
				for (size_t j = 0; j < centers[k].size(); ++j)
					centers[k][j] += m_samples[i][j];
			}
			for (size_t j = 0; j < centers[k].size(); ++j)
				centers[k][j] /= m_clusters[k].size();
		}
	}

	void KMeansClu::updateClusterCenters(const std::vector<std::vector<Real>> &centers) {
		m_cluster_centers.resize(m_clusters.size());
		for (size_t k = 0; k < m_clusters.size(); ++k) {
			m_cluster_centers[k] = m_clusters[k][0];
			Real min_distance = euclideanDistance(centers[k].begin(), centers[k].end(), m_samples[m_clusters[k][0]].begin());
			for (size_t i = 1; i < m_clusters[k].size(); ++i) {
				Real distance = euclideanDistance(centers[k].begin(), centers[k].end(), m_samples[m_clusters[k][i]].begin());
				if (distance < min_distance) {
					min_distance = distance;
					m_cluster_centers[k] = m_clusters[k][i];
				}
			}
		}
	}



	KmeansPoint::KmeansPoint(size_t point_id, std::vector<Real>& p_values, std::string p_name) {
		id_point = point_id;
		total_values = p_values.size();

		for (size_t i = 0; i < total_values; i++)
			values.push_back(p_values[i]);

		name = p_name;
		id_cluster = -1;
	}

	KmeansCluster::KmeansCluster(size_t cluster_id, KmeansPoint point) {
		id_cluster = cluster_id;

		size_t total_values = point.getTotalValues();

		for (size_t i = 0; i < total_values; i++)
			central_values.push_back(point.getValue(i));

		points.push_back(point);
	}

	void KmeansCluster::addPoint(KmeansPoint point) {
		points.push_back(point);
	}

	bool KmeansCluster::removePoint(size_t id_point) {
		size_t total_points = points.size();

		for (size_t i = 0; i < total_points; i++){
			if (points[i].getID() == id_point){
				points.erase(points.begin() + i);
				return true;
			}
		}
		return false;
	}

	Real KmeansCluster::getCentralValue(size_t index) {
		return central_values[index];
	}

	void KmeansCluster::setCentralValue(size_t index, Real value) {
		central_values[index] = value;
	}

	KmeansPoint KmeansCluster::getPoint(size_t index) {
		return points[index];
	}

	size_t KmeansCluster::getTotalPoints() {
		return points.size();
	}

	size_t KmeansCluster::getID() {
		return id_cluster;
	}


	KMeans::KMeans(size_t k, size_t num_points, size_t dims, size_t iterations):K(k),total_points(num_points),total_values(dims),max_iterations(iterations) {
		
	}

	size_t KMeans::getIDNearestCenter(KmeansPoint point) {
		Real sum = 0.0, min_dist;
		size_t id_cluster_center = 0;

		for (size_t i = 0; i < total_values; i++) {
			sum += std::pow(clusters[0].getCentralValue(i) - point.getValue(i), 2.0);
		}

		min_dist = std::sqrt(sum);

		for (size_t i = 1; i < K; i++) {
			Real dist;
			sum = 0.0;

			for (size_t j = 0; j < total_values; j++) {
				sum += std::pow(clusters[i].getCentralValue(j) - point.getValue(j), 2.0);
			}

			dist = std::sqrt(sum);

			if (dist < min_dist) {
				min_dist = dist;
				id_cluster_center = i;
			}
		}

		return id_cluster_center;
	}

	void KMeans::run(std::vector<KmeansPoint>& points, Random* rnd) {
		if (K > total_points)
			return;

		std::vector<size_t> prohibited_indexes;

		// choose K distinct values for the centers of the clusters
		for (size_t i = 0; i < K; i++) {
			while (true) {
				size_t index_point = (size_t)std::floor(total_points * rnd->uniform.next());

				if (std::find(prohibited_indexes.begin(), prohibited_indexes.end(), index_point) == prohibited_indexes.end()) {
					prohibited_indexes.push_back(index_point);
					points[index_point].setCluster(i);
					KmeansCluster cluster(i, points[index_point]);
					clusters.push_back(cluster);
					break;
				}
			}
		}

		size_t iter = 1;

		while (true) {
			bool done = true;

			// associates each point to the nearest center
			for (size_t i = 0; i < total_points; i++) {
				size_t id_old_cluster = points[i].getCluster();
				size_t id_nearest_center = getIDNearestCenter(points[i]);

				if (id_old_cluster != id_nearest_center) {
					if (id_old_cluster != -1)
						clusters[id_old_cluster].removePoint(points[i].getID());

					points[i].setCluster(id_nearest_center);
					clusters[id_nearest_center].addPoint(points[i]);
					done = false;
				}
			}

			// recalculating the center of each cluster
			for (size_t i = 0; i < K; i++) {
				for (size_t j = 0; j < total_values; j++) {
					size_t total_points_cluster = clusters[i].getTotalPoints();
					Real sum = 0.0;

					if (total_points_cluster > 0) {
						for (size_t p = 0; p < total_points_cluster; p++)
							sum += clusters[i].getPoint(p).getValue(j);
						clusters[i].setCentralValue(j, sum / total_points_cluster);
					}
				}
			}

			if (done == true || iter >= max_iterations) {
				std::cout << "Break in iteration " << iter << std::endl;
				break;
			}

			iter++;
		}
	}

	void KMeans::run(std::vector<std::vector<Real>>& all_data, Random* rnd) {
		if (K > total_points)
			return;

		std::vector<size_t> prohibited_indexes;
		std::vector<KmeansPoint> points;
		for (size_t i = 0; i < all_data.size(); ++i) {
			auto value = all_data[i];
			KmeansPoint p(i, value);
			points.emplace_back(p);
		}

		// choose K distinct values for the centers of the clusters
		for (size_t i = 0; i < K; i++) {
			while (true) {
				size_t index_point = (size_t)std::floor(total_points * rnd->uniform.next());

				if (std::find(prohibited_indexes.begin(), prohibited_indexes.end(), index_point) == prohibited_indexes.end()) {
					prohibited_indexes.push_back(index_point);
					points[index_point].setCluster(i);
					KmeansCluster cluster(i, points[index_point]);
					clusters.push_back(cluster);
					break;
				}
			}
		}

		size_t iter = 1;

		while (true) {
			bool done = true;

			// associates each point to the nearest center
			for (size_t i = 0; i < total_points; i++) {
				size_t id_old_cluster = points[i].getCluster();
				size_t id_nearest_center = getIDNearestCenter(points[i]);

				if (id_old_cluster != id_nearest_center) {
					if (id_old_cluster != -1)
						clusters[id_old_cluster].removePoint(points[i].getID());

					points[i].setCluster(id_nearest_center);
					clusters[id_nearest_center].addPoint(points[i]);
					done = false;
				}
			}

			// recalculating the center of each cluster
			for (size_t i = 0; i < K; i++) {
				for (size_t j = 0; j < total_values; j++) {
					size_t total_points_cluster = clusters[i].getTotalPoints();
					Real sum = 0.0;

					if (total_points_cluster > 0) {
						for (size_t p = 0; p < total_points_cluster; p++)
							sum += clusters[i].getPoint(p).getValue(j);
						clusters[i].setCentralValue(j, sum / total_points_cluster);
					}
				}
			}

			if (done == true || iter >= max_iterations) {
				std::cout << "Break in iteration " << iter << std::endl;
				break;
			}

			iter++;
		}
	}

	void KMeans::showClusters() {
		// shows elements of clusters
		for (size_t i = 0; i < K; i++) {
			size_t total_points_cluster = getClusters()[i].getTotalPoints();

			std::cout << "Cluster " << getClusters()[i].getID() + 1 << std::endl;
			for (size_t j = 0; j < total_points_cluster; j++) {
				std::cout << "Point " << getClusters()[i].getPoint(j).getID() + 1 << ": ";
				for (size_t p = 0; p < total_values; p++)
					std::cout << clusters[i].getPoint(j).getValue(p) << " ";

				std::string point_name = getClusters()[i].getPoint(j).getName();

				if (point_name != "")
					std::cout << "- " << point_name;

				std::cout << std::endl;
			}

			std::cout << "Cluster center values: ";

			for (size_t j = 0; j < total_values; j++)
				std::cout << getClusters()[i].getCentralValue(j) << " ";

			std::cout << std::endl;
		}
	}
}