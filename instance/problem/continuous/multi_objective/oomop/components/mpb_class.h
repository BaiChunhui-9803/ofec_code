
/*************************************************************************
*  this class define a MPB, which includes peaks with different peaks
*************************************************************************/


// created by tanqingshan on June 19th, 2022.

#ifndef OFEC_MPB_CLASS_H
#define OFEC_MPB_CLASS_H

#include"peak_shapes.h"

namespace ofec {
	class Mpb_class {
	private:
		std::vector<std::shared_ptr<Peak>> m_peaks;
	public:
		Mpb_class(const std::vector<std::shared_ptr<Peak>>& peaks):m_peaks(peaks) {}
		//std::vector<std::shared_ptr<Peak>> getPeaks() { return m_peaks; }
		std::shared_ptr<Peak> getPeak(size_t idx) { return m_peaks[idx]; }
		std::vector<size_t> getMaxPeakIndex() {
			Real temp = -1. * 10e14;
			std::vector<size_t> idx;
			for (size_t i = 0; i < m_peaks.size(); ++i) {
				if (m_peaks[i]->getHeight() > temp) {
					temp = m_peaks[i]->getHeight();
					//idx = i;
				}
			}
			for (size_t i = 0; i < m_peaks.size(); ++i) {
				if (m_peaks[i]->getHeight() == temp) {
					idx.push_back(i);
				}
			}
			return idx;
		}

		std::vector<size_t> getMinPeakIndex() {
			Real temp = 1. * 10e14;
			std::vector<size_t> idx;
			for (size_t i = 0; i < m_peaks.size(); ++i) {
				if (m_peaks[i]->getHeight() < temp) {
					temp = m_peaks[i]->getHeight();
				}
			}
			for (size_t i = 0; i < m_peaks.size(); ++i) {
				if (m_peaks[i]->getHeight() == temp) {
					idx.push_back(i);
				}
			}
			return idx;
		}

		Real calMpbValue(const std::vector<Real>& sol) {
			std::vector<Real> peak_values;
			for (size_t i = 0; i < m_peaks.size(); ++i) {
				peak_values.push_back(m_peaks[i]->calPeakValue(sol));
			}
			return maxElement(peak_values);
		}
	};
}

#endif // !OFEC_MPB_CLASS_H




