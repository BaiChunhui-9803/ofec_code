#include "elementary_function.h"


void ofec::fun_gen::elementaryFunctionParameters::resetFun(std::unique_ptr<elementary_function>& fun)const
{
	switch (m_type)
	{

	case elementary_fun_type::constant_fun:
		fun.reset(new constant_fun(m_from_x_range, m_to_y_range));
		break;
	case elementary_fun_type::random_fun:
		fun.reset(new random_fun(m_sample_num, m_from_x_range, m_to_y_range));
		break;
	case elementary_fun_type::cos_fun:
		fun.reset(new cos_fun(m_fun_num,m_numPI, m_from_x_range, m_to_y_range));
		break;
	default:
		break;
	}
}



void ofec::fun_gen::elementary_function_initialize(std::unique_ptr<elementary_function>& fun, elementaryFunctionParameters par, RandBase& randomBase)
{
	par.resetFun(fun);
	fun->initialize(randomBase);
}


