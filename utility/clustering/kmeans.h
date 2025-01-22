#ifndef OFEC_KMEANS_CLUSTERING_H
#define OFEC_KMEANS_CLUSTERING_H

#include "../../core/problem/encoding.h"
#include "../../utility/random/newran.h"
#include "../../core/problem/solution.h"
#include "../../core/problem/problem.h"

namespace ofec {
	class KMeansClu {
	private:
		std::vector<std::vector<Real>> m_samples;
		std::vector<std::vector<size_t>> m_clusters;
		std::vector<size_t> m_cluster_centers;
		Real m_sum_inner_distance;

		void updateClusters(const std::vector<std::vector<Real>> &centers);
		void updateCenters(std::vector<std::vector<Real>> &centers);
		void updateClusterCenters(const std::vector<std::vector<Real>> &centers);

	public:
		void setData(const std::vector<const SolutionBase*> sols, Problem *pro);
		void clustering(size_t K);
		void clustering(const std::list<const SolutionBase*> &init_central_sols);
		const std::vector<std::vector<size_t>>& clusters() const { return m_clusters; }
		const std::vector<size_t>& clusterCenters() const { return m_cluster_centers; }
	};




	//another realization
	class KmeansPoint{
	private:
		size_t id_point, id_cluster;
		std::vector<Real> values;
		size_t total_values;// the number dimensions
		std::string name;

	public:
		KmeansPoint(size_t point_id, std::vector<Real>& p_values, std::string p_name = "");

		size_t getID() {
			return id_point;
		}

		void setCluster(size_t cluster_id) {
			id_cluster = cluster_id;
		}

		size_t getCluster() {
			return id_cluster;
		}

		Real getValue(size_t index) {
			return values[index];
		}

		size_t getTotalValues() {
			return total_values;
		}

		void addValue(Real value) {
			values.push_back(value);
		}

		std::string getName() {
			return name;
		}
	};

	class KmeansCluster{
	private:
		size_t id_cluster;
		std::vector<Real> central_values;
		std::vector<KmeansPoint> points;

	public:
		KmeansCluster(size_t cluster_id, KmeansPoint point);

		void addPoint(KmeansPoint point);

		bool removePoint(size_t id_point);

		Real getCentralValue(size_t index);

		void setCentralValue(size_t index, Real value);

		KmeansPoint getPoint(size_t index);

		size_t getTotalPoints();

		size_t getID();
	};

	class KMeans{
	private:
		size_t K; // the number of clusters
		size_t total_values;// the number of dimension
		size_t total_points;// the number of points
		size_t max_iterations;// the number of iterations
		std::vector<KmeansCluster> clusters;

		// return ID of nearest center (uses euclidean distance)
		size_t getIDNearestCenter(KmeansPoint point);

	public:
		KMeans(size_t k, size_t num_points, size_t total_dims, size_t iterations);

		void run(std::vector<KmeansPoint>& points, Random* rnd);

		void run(std::vector<std::vector<Real>>& all_data, Random* rnd);

		std::vector<KmeansCluster> getClusters() {
			return clusters;
		}

		void showClusters();
	};
}

#endif // !OFEC_KMEANS_CLUSTERING_H
