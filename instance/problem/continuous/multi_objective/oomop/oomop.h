
/*********************************************************************
* This MOP is constructed by object-oriented construction method
**********************************************************************/

// created by tanqingshan on June 20th, 2022

#ifndef OFEC_OOMOP_H
#define OFEC_OOMOP_H

#include "components/mpb_class.h"
#include "components/ps_shapes.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"
#include "../../../../../core/global.h"

namespace ofec {
	enum class PSType { line, circle, power, sin, sigmoid };
	enum class PFType { linear, convex, extreConvex, concave, extreConcave, mix, multiple_mix, discontinuity };

	class OOMOP : public Continuous, public MetricsMOP {
	protected:
		std::vector<std::shared_ptr<Mpb_class>> m_mpb_pri;
		std::vector<std::shared_ptr<PS_shape>> m_ps_shapes;
		std::vector<std::vector<std::shared_ptr<Peak>>> m_ps_peaks;
		Real m_max_value_of_last_obj = 0;
		std::vector<std::pair<Real, Real>> m_range_optima;

		//parameters
		std::vector<int> m_num_pri_var;
		size_t m_num_pub_var;
		std::vector<int> m_pri_peak_num;
		std::vector<Real> m_pri_peak_height_rate;//for each obj
		std::vector<Real> m_pri_peak_slope_range;//common range
		std::vector<Real> m_pri_peak_width_range;//common range
		std::vector<Real> m_pri_norm;//for each obj
		std::vector<Real> m_pub_norm;
		size_t m_num_ps;
		size_t m_num_global_ps;
		std::vector<Real> m_global_ps_span;//for each global ps
		std::vector<Real> m_ps_decay_rate;//for each class of global ps
		std::vector<Real> m_ps_radiate_slope;//for each ps

		std::vector<int> m_ps_type;//for each global ps
		std::vector<int> m_pf_type;//for each global ps

		std::vector<std::pair<Real, size_t>> m_ps_record;//record ps decay rate and attached global ps

	public:
		OOMOP() {}
		void initialize_() override;
		Real getMaxValueLastObj() { return m_max_value_of_last_obj; }
		void setMaxValueLastObj(Real v) { m_max_value_of_last_obj = v; }
		
		std::vector<std::pair<size_t, bool>> attachPS(std::vector<Real>& sol);

		std::vector<int> getPriVarNum() { return m_num_pri_var; }
		size_t getPubVarNum() { return m_num_pub_var; }
		std::vector<std::shared_ptr<Mpb_class>> getMpbPri() { return m_mpb_pri; }
		std::vector<std::shared_ptr<PS_shape>> getPShape() { return m_ps_shapes; }
		std::vector<std::vector<std::shared_ptr<Peak>>> getPSPeak() { return m_ps_peaks; }

		std::vector<std::pair<Real, Real>> getRangeOptima() { return m_range_optima; }
		size_t getGlobalPSNum() { return m_num_global_ps; }

		bool checkBoundary(std::vector<std::pair<Real,Real>>& dim_bounds,std::vector<Real>& point);

		void sample2HighDim(std::vector<std::vector<Real>>& low_dim_samples, Real offset, size_t sample_num, size_t dim_index);

		virtual void generatePF(){}
		virtual void calPubValue(std::vector<Real>& x, std::vector<Real>& pub_objs){}
		virtual std::vector<std::vector<Real>> samplePubOptima(size_t dim_num,size_t num_g_ps){
			std::vector<std::vector<Real>> temp;
			std::vector<Real> temp_sol;
			temp_sol.push_back(0.);
			temp.push_back(temp_sol);
			return temp;
		}
	};
}

#endif // !OFEC_OOMOP_H