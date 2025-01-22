#include "moving_peak.h"
#include "../../../../../../core/global.h"
#include <algorithm>


namespace ofec {
	static Real basis_peak[5][7] = {
	  {8.0,  64.0,  67.0,  55.0,   4.0, 0.1, 50.0},
	  {50.0,  13.0,  76.0,  15.0,   7.0, 0.1, 50.0},
	  {9.0,  19.0,  27.0,  67.0, 24.0, 0.1, 50.0},
	  {66.0,  87.0,  65.0,  19.0,  43.0, 0.1, 50.0},
	  {76.0,  32.0,  43.0,  54.0,  65.0, 0.1, 50.0}
	};

	static Real twin_peak[7] = /* difference to first peak */
	{
	  1.0,  1.0,  1.0,  1.0,   1.0, 0.0, 0.0,
	};

	void MovingPeak::setSeverity() {
		for (int i = 0; i < m_peak.size(); i++) {
			m_height_severity[i] = m_random->uniform.nextNonStd(1.0,10.);//1.+9.*i/m_peak.size();//severity of height changes, larger numbers  mean larger severity. in the contex of ROOT, peaks have different values
			m_width_severity[i] = m_random->uniform.nextNonStd(0.1, 1.0);//1.+9.*i/m_peak.size();//severity of height changes, larger numbers  mean larger severity. in the contex of ROOT, peaks have different values
		}

	}

	//MovingPeak::MovingPeak(const ParameterMap &v) :Problem(v.at("problem name")), \
	//	UncertaintyContinuous(v.at("problem name"), v.at("number of variables"), v.at("numPeak")) {
	//

	//}
	//MovingPeak::MovingPeak(const std::string &name, int rDimNumber,  int rNumPeaks, Real  rChangingRatio,  int fre, Real vlength, bool rFlagDimChange, \
	//	 bool rFlagNumPeakChange,  int peakNumChangeMode,  bool flagNoise,  bool flagTimelinkage) :\
	//	Problem(name), UncertaintyContinuous(name, rDimNumber, rNumPeaks, 1) {

	//	setFrequency(fre);
	//	setFlagVarMemChange(rFlagDimChange);
	//	setNumPeakChangeMode(peakNumChangeMode);
	//	setFlagNumPeaksChange(rFlagNumPeakChange);
	//	setFlagNoisyFromVariable(flagNoise);
	//	setFlagTimeLinkage(flagTimelinkage);
	//	set_vlength(vlength);
	//	setNumChange(rChangingRatio);
	//	m_optima->setVariableGiven(true);
	//	m_num_peak_tracked = 0;
	//	
	//	initialize();
	//}
	void MovingPeak::initializeParameters() {
		int i = 0;
		m_variable_accuracy = 0.1;
		m_objective_accuracy=0.2;
		setDomain(0, 100);
		setInitialDomain(0, 100);
		setOptMode(OptimizeMode::kMaximize,0);
		addTag(ProblemTag::kMMOP);
		/***************************************
		//		m_F		Evaluation Function
		//		1		constant_basis_func()
		//		2		five_peak_basis_func()
		//		3		peak_function1()
		//		4		peak_function1()
		//		5		peak_function_hilly()
		//		6		peak_function_twin()
		**************************************
		in>>temp>>m_vlength; // distance by which the peaks are moved, severity
		lambda determines whether there is a direction of the movement, or whether
		they are totally random. For lambda = 1.0 each move has the same direction,
		while for lambda = 0.0, each move has a random direction
		//in>>temp>>lambda;
		//in>>temp>>m_useBasisFunction;  if set to 1, a static landscape (basis_function) is included in the fitness evaluation
		}*/

		m_F = 4;
		m_lambda = 0;
		//m_vlength = 1.0;
		m_min_height =30.0;
		m_max_height = 70.0;
		m_standardHeight = 50.0;
		m_min_width = 1.;
		m_max_width = 12.;
		m_standardWidth = 0.0;

		m_shift.resize(m_number_variables);
		m_pre_movement.resize(m_num_peaks, std::vector<Real>(m_number_variables));

		setSeverity();

		updateTimeLinkage();
    		for (i = 0; i < m_peak.size(); i++) {
			for (int j = 0; j < m_number_variables; j++) {
				m_peak[i][j] = 100.0 * m_random->uniform.next();
				m_pre_movement[i][j] = m_random->uniform.next() - 0.5;
			}
		}

		if (m_standardHeight <= 0.0) {
			for (i = 0; i < m_peak.size(); i++) m_height[i] = (m_max_height - m_min_height) *
				m_random->uniform.next()+ m_min_height;
		}
		else {
			for (i = 0; i < m_peak.size(); i++) m_height[i] = m_standardHeight;
		}

		if (m_standardWidth <= 0.0) {
			for (i = 0; i < m_peak.size(); i++)
				m_width[i] = (m_max_width - m_min_width)* m_random->uniform.next()+ m_min_width;
		}
		else {
			for (i = 0; i < m_peak.size(); i++)
				m_width[i] = m_standardWidth;
		}

		calculateGlobalOptima();

		for (i = 0; i < m_peak.size(); i++) m_num_tracking[i] = 0;
		for (i = 0; i < m_peak.size(); i++)
			std::copy(m_peak[i].begin(), m_peak[i].end(),m_pre_peak[i].begin());

		std::copy(m_height.begin(), m_height.end(), m_pre_height.begin());
		std::copy(m_width.begin(), m_width.end(), m_pre_width.begin());
	}

	void MovingPeak::initialize_()	{
		UncertaintyContinuous::initialize_();
		auto v(*m_param);

		
		if (v.has("shiftLength"))
			m_vlength = v.get<Real>("shiftLength");
		else 
			m_vlength = m_random->normal.next();
		setType(CT_Random);
		m_num_peak_tracked = 0;
		initializeParameters();
	}

	/* current_peak_calc determines the peak of the current best Solution */
	void MovingPeak::calculateCurrentPeak(const Real *gen) {
		int i;
		Real maximum = -100000.0, dummy;

		m_current_peak = 0;
		maximum = selectFunction(gen, 0);
		for (i = 1; i < m_peak.size(); i++) {
			dummy = selectFunction(gen, i);
			if (dummy > maximum) {
				maximum = dummy;
				m_current_peak = i;
			}
		}
	}

	void MovingPeak::evaluateObjective(Real* x_, std::vector<Real>& obj_) {
		Real *x = new Real[m_number_variables];
		for (size_t i = 0; i < m_number_variables; i++)
			x[i] = x_[i];
		if (this->m_flag_noisy_from_variable)	addNoise(x);

		Real maximum = -std::numeric_limits<Real>::max(), dummy;

		for (int i = 0; i < m_peak.size(); i++) {
			dummy = selectFunction(x, i);
			if (dummy > maximum)      maximum = dummy;
		}
		obj_[0] = maximum;

		delete[] x;
		x = 0;
	}


	/* dummy evaluation function allows to evaluate without being counted */
	Real MovingPeak::dummyEval(const Real *gen) {
		int i;
		Real maximum = -100000.0, dummy;

		for (i = 0; i < m_peak.size(); i++) {
			dummy = selectFunction(gen, i);
			if (dummy > maximum)    maximum = dummy;
		}
		return(maximum);
	}

	/* whenever this function is called, the peaks are changed */
	void MovingPeak::changeRandom() {
		int i = 0, j = 0;
		Real sum, sum2, offset;

		
		for (i = 0; i < m_peak.size(); i++) 	m_pre_peak[i] = m_peak[i];
		m_pre_height = m_height;
		m_pre_width = m_width;

		for (i = 0; i < m_peak.size(); i++) {
			if (m_flag_change[i] == false) continue;
			/* shift peak locations */
			sum = 0.0;
			for (j = 0; j < m_number_variables; j++) {
				m_shift[j] = m_random->uniform.next() - 0.5;
				sum += m_shift[j] * m_shift[j];
			}
			if (sum > 0.0)		sum = m_vlength / sqrt(sum);
			else  sum = 0.0;                        /* only in case of rounding errors */

			sum2 = 0.0;
			for (j = 0; j < m_number_variables; j++) {
				m_shift[j] = sum*(1.0 - m_lambda)*m_shift[j] + m_lambda*m_pre_movement[i][j];
				sum2 += m_shift[j] * m_shift[j];
			}
			if (sum2 > 0.0)sum2 = m_vlength / sqrt(sum2);
			else     sum2 = 0.0;                      /* only in case of rounding errors */

			for (j = 0; j < m_number_variables; j++) {
				m_shift[j] *= sum2;
				m_pre_movement[i][j] = m_shift[j];
				if (m_domain[j].limit.first > (m_peak[i][j] + m_pre_movement[i][j])) {
					m_peak[i][j] = 2.0*m_domain[j].limit.first - m_peak[i][j] - m_pre_movement[i][j];
					m_pre_movement[i][j] *= -1.0;
				}
				else if (m_domain[j].limit.second < (m_peak[i][j] + m_pre_movement[i][j])) {
					m_peak[i][j] = 2.0*m_domain[j].limit.second - m_peak[i][j] - m_pre_movement[i][j];
					m_pre_movement[i][j] *= -1.0;
				}
				else
					m_peak[i][j] += m_pre_movement[i][j];
			}

			/* change peak width */
			offset = m_random->normal.next() * m_width_severity[i];
			if ((m_width[i] + offset) < m_min_width)		m_width[i] = 2.0*m_min_width - m_width[i] - offset;
			else if ((m_width[i] + offset) > m_max_width)	m_width[i] = 2.0*m_max_width - m_width[i] - offset;
			else	m_width[i] += offset;

			if (m_counter > 1 && m_ratio_changing_peak < 1.0&&m_flag_global_optima[i]) continue;
			/* change peak height */

			offset = m_height_severity[i] * m_random->normal.next();

			if ((m_height[i] + offset) < m_min_height)	m_height[i] = 2.0*m_min_height - m_height[i] - offset;
			else if ((m_height[i] + offset) > m_max_height)	m_height[i] = 2.0*m_max_height - m_height[i] - offset;
			else	m_height[i] += offset;
		}
		calculateGlobalOptima();
		updateNumChange();
	}


	/* Basis Functions */

	/* This gives a constant value back to the eval-function that chooses the max of them */
	Real MovingPeak::constantBasisFunc(const Real *gen) {
		return 0.0;
	}

	Real MovingPeak::fivePeakBasisFunc(const Real *gen) {
		Real maximum = -100000.0, dummy = 0;
		for (int i = 0; i < 5; i++) {
			dummy = (gen[0] - basis_peak[i][0])*(gen[0] - basis_peak[i][0]);
			for (int j = 1; j < m_number_variables; j++)  dummy += (gen[j] - basis_peak[i][j])*(gen[j] - basis_peak[i][j]);
			dummy = basis_peak[i][m_number_variables + 1] - (basis_peak[i][m_number_variables] * dummy);
			if (dummy > maximum)       maximum = dummy;
		}
		return maximum;
	}

	/* Peak Functions */

	/* sharp peaks */
	Real MovingPeak::peakFunction1(const Real *gen, int peak_number) {

		Real dummy = (gen[0] - m_peak[peak_number][0])*(gen[0] - m_peak[peak_number][0]);
		for (int j = 1; j < m_number_variables; j++)
			dummy += (gen[j] - m_peak[peak_number][j])*(gen[j] - m_peak[peak_number][j]);

		return m_height[peak_number] / (1 + m_width[peak_number] * dummy);
	}

	Real MovingPeak::peakFunctionHilly(const Real *gen, int peak_number) {
		int j = 0;
		Real dummy = (gen[0] - m_peak[peak_number][0])*(gen[0] - m_peak[peak_number][0]);
		for (j = 1; j < m_number_variables; j++)
			dummy += (gen[j] - m_peak[peak_number][j])*(gen[j] - m_peak[peak_number][j]);

		return m_height[peak_number] - m_width[peak_number] * dummy - 0.01*sin(20.0*dummy);
	}

	Real MovingPeak::peakFunctionTwin(const Real  *gen, int peak_number) /* two twin peaks moving together */
	{
		int j;
		Real maximum = -100000.0, dummy = pow(gen[0] - m_peak[peak_number][0], 2);
		for (j = 1; j < m_number_variables; j++)
			dummy += pow(gen[j] - m_peak[peak_number][j], 2);

		dummy = m_height[peak_number] - m_width[peak_number] * dummy;

		maximum = dummy;
		dummy = pow(gen[0] - (m_peak[peak_number][0] + twin_peak[0]), 2);
		for (j = 1; j < m_number_variables; j++)
			dummy += pow(gen[j] - (m_peak[peak_number][j] + twin_peak[0]), 2);

		dummy = m_height[peak_number] + twin_peak[m_number_variables + 1] - ((m_width[peak_number] + twin_peak[m_number_variables])*dummy);
		if (dummy > maximum)
			maximum = dummy;

		return maximum;
	}


	Real MovingPeak::peakFunctionCone(const Real *gen, const int &peak_number) {

		Real val, dummy = 0;
		for (int j = 0; j < m_number_variables; j++) {
			val = gen[j] - m_peak[peak_number][j];
			dummy += val*val;
		}
		if (dummy != 0)  dummy = m_height[peak_number] - m_width[peak_number] * sqrt(dummy);
		else dummy = m_height[peak_number];
		return dummy;
	}
	Real MovingPeak::selectFunction(const Real  *gen, const int &peak_number) {
		Real dummy = 0;
		switch (m_F) {
		case 1: {
			dummy = constantBasisFunc(gen);
			break;
		}
		case 2: {
			dummy = fivePeakBasisFunc(gen);
			break;
		}
		case 3: {
			dummy = peakFunction1(gen, peak_number);
			break;
		}
		case 4: {
			dummy = peakFunctionCone(gen, peak_number);
			break;
		}
		case 5: {
			dummy = peakFunctionHilly(gen, peak_number);
			break;
		}
		case 6: {
			dummy = peakFunctionTwin(gen, peak_number);
			break;
		}
		}
		return dummy;
	}
	/* The following procedures may be used to change the Step size over time */


	void MovingPeak::changeStepsizeRandom() /* assigns vlength a value from a normal distribution */
	{
		m_vlength = m_random->normal.next();
	}

	void MovingPeak::changeStepsizeLinear() /* sinusoidal change of the stepsize, */
	{
		static	thread_local std::unique_ptr< int> counter;
		if (!counter.get()) counter.reset(new int(1));
		thread_local std::unique_ptr<Real> frequency;
		if (!frequency.get()) frequency.reset(new Real(3.14159 / 20.0));

		m_vlength = 1 + sin((Real)(*counter)*(*frequency));
		(*counter)++;
	}

	int MovingPeak::getRightPeak()  /* returns 1 if current best Solution is on highest m_peak, 0 otherwise */
	{
		bool flag = false;

		for (int i = 0; i < m_peak.size(); i++) {
			if (m_flag_global_optima[i] == true && m_current_peak == i) {
				flag = true;
				break;
			}
		}

		return flag;
	}
	void MovingPeak::setVlength(Real s) {
		m_vlength = s;
		m_params["vlength"] = s;
	}

	void MovingPeak::changeNumPeak() {

		m_num_peaks = m_temp_num_peak;
		MovingPeak mpb(*this);

		mpb.initialize();

		/*
		* 	//MovingPeak::MovingPeak(const std::string &name, int rDimNumber,  int rNumPeaks, Real  rChangingRatio,  int fre, Real vlength, bool rFlagDimChange, \
	//	 bool rFlagNumPeakChange,  int peakNumChangeMode,  bool flagNoise,  bool flagTimelinkage) :\
	//	Problem(name), UncertaintyContinuous(name, rDimNumber, rNumPeaks, 1)
		*/
		//MovingPeak mpb(m_name, m_number_variables, m_temp_num_peak, m_ratio_changing_peak, m_flag_variable_memory_change
		//	, m_flag_num_peak_change, m_mode, m_flag_noisy_from_variable, m_flag_time_linkage);
		mpb.copy(*this);
		mpb.calculateGlobalOptima();

		*this = std::move(mpb);

	}

	void MovingPeak::copy(const Problem &rP) {
		UncertaintyContinuous::copy(rP);

		auto& mpb = dynamic_cast<const MovingPeak &>(rP);
		
		int peaks = m_peak.size() < mpb.getNumPeak() ? m_peak.size() : mpb.getNumPeak();

		m_F = mpb.m_F;
		m_vlength = mpb.m_vlength;
		m_lambda = mpb.m_lambda;
		m_standardHeight = mpb.m_standardHeight;
		m_standardWidth = mpb.m_standardWidth;
		m_shift = mpb.m_shift;
				
		for (int i = 0; i < peaks; i++) {
			m_pre_movement[i] = mpb.m_pre_movement[i];
		}

	}

	Real MovingPeak::getVlength() {
		return m_vlength;
	}

	void MovingPeak::getDynamicObjective(const SolutionBase& s, std::vector<Real>& objs, int t)
	{
	}

	void MovingPeak::getNoisyObjective(const SolutionBase& s, std::vector<Real>& objs, Random *rnd)
	{
	}
}
