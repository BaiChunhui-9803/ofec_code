
/*
* this file define different PS shapes which can be initialed in OOMOP
*/

/*
* created by tanqingshan on June 19th, 2022.
*/

#ifndef OFEC_PS_SHAPES_H
#define OFEC_PS_SHAPES_H

#include<vector>
#include "../../../../../../utility/functional.h"


namespace ofec {
	class PS_shape {
	protected:
		std::vector<std::vector<Real>> m_points;//the start and end point
		std::vector<Real> m_offset;
		Real m_ratio;   //ps decay rate
		size_t m_shape_flag;
	public:
		PS_shape(const std::vector<std::vector<Real>>& points ,const std::vector<Real>& offset,Real ratio,size_t flag):m_points(points),m_offset(offset),m_ratio(ratio),m_shape_flag(flag){}
		
		virtual std::vector<Real> calValue(const Real& offset) {
			std::vector<Real> temp = { 0.1 };
			return temp;
		}
		Real calOneDimValue(const Real& offset, size_t idx);
		bool if_in_range(std::vector<Real>& sol);
		std::pair<std::vector<Real>, Real> dist_from_PS(const std::vector<Real>& sol, const std::vector<Real>& p_norms);

		std::vector<std::vector<Real>> getPoints() { return m_points; }
		void setpoints(const std::vector<std::vector<Real>>& points) {
			m_points = points;
		}
		std::vector<Real> getOffset() { return m_offset; }
		void setOffset(const std::vector<Real>& offset) {
			m_offset = offset;
		}
		Real getRatio() { return m_ratio; }
		void setRatio(const Real& ratio) {
			m_ratio = ratio;
		}
		size_t getShapeFlag() { return m_shape_flag; }
		void setShapeFlag(const size_t flag) {
			m_shape_flag = flag;
		}
	};

	class Line_shape:public PS_shape {
	private:
		std::vector<std::pair<Real, Real>> m_ps_boundary;
	public:
		Line_shape(const std::vector<std::vector<Real>>& points, const std::vector<Real>& offset, Real ratio,size_t flag=1) :\
			PS_shape(points, offset,ratio,flag){

		}
	
		//calculate the value of other dimensions, other ps shape need to design again
		std::vector<Real> calValue(const Real& offset) override;
	};

	class Trigo_shape :public PS_shape {
	private:
		Real m_amp;     //the amplitude of the trigo
		Real m_fre;     //the period 
		std::vector<std::pair<Real, Real>> m_ps_boundary;
	public:
		Trigo_shape(const std::vector<std::vector<Real>>& points, Real amp, Real fre, const std::vector<Real>& offset, Real ratio, size_t flag = 2) :PS_shape(points, offset, ratio,flag), \
			m_amp(amp), m_fre(fre){}

		//calculate the value of other dimensions, other ps shape need to design again
		std::vector<Real> calValue(const Real& offset);
		Real getAmp() { return m_amp; }
		void setAmp(const Real& amp) {
			m_amp = amp;
		}
		Real getFre() { return m_fre; }
		void setFre(const Real& fre) {
			m_fre = fre;
		}
	};

	class Power_shape :public PS_shape {
	private:
		Real m_power;     //the power of the function
		Real m_scale;     //the coefficient of the function 
		std::vector<std::pair<Real, Real>> m_ps_boundary;
	public:
		Power_shape(const std::vector<std::vector<Real>>& points, Real power, Real scale,const std::vector<Real>& offset, Real ratio,size_t flag=3) :\
			PS_shape(points, offset,ratio,flag), m_power(power), m_scale(scale){}

		//calculate the value of other dimensions, other ps shape need to design again
		std::vector<Real> calValue(const Real& offset);
		Real getPower() { return m_power; }
		void setPower(const Real& power) {
			m_power = power;
		}
		Real getScale() { return m_scale; }
		void setScale(const Real& scale) {
			m_scale = scale;
		}
	};


	/*class Circle_shape:public PS_shape {
	private:
		std::vector<Real> m_center;
		Real m_r;
		Real m_thick;
		size_t m_index;
		size_t m_shape_flag=8;
	public:
		Circle_shape(const std::vector<Real>& center, Real r, Real thick, size_t idx) :m_center(center), m_r(r), m_thick(thick), m_index(idx) {}
		bool if_in_range(std::vector<Real>& sol) override;
		std::vector<Real> getCenter() { return m_center; }
		Real getR() { return m_r; }
		Real getThick() { return m_thick; }
		size_t getIndex() { return m_index; }
		void setCenter(const std::vector<Real>& center) {
			m_center = center;
		}
		void setR(const Real& r) {
			m_r = r;
		}
		void setThick(const Real& thick) {
			m_thick = thick;
		}
		void setIndex(const size_t& idx) {
			m_index = idx;
		}
	};*/
}

#endif // !OFEC_PS_SHAPES_H