#pragma once
#include "ioh/problem/utils.hpp"
#include "pbo_problem.hpp"

namespace ioh
{
    namespace problem
    {
        namespace pbo
        {
            //! OneMaxDummy2 problem id 5
            class OneMaxDummy2 final: public PBOProblem<OneMaxDummy2>
            {
                std::vector<int> info_;
            protected:
                //! Evaluation method
                double evaluate(const std::vector<int> &x) override
                {
                    auto result = 0.0;
                    for (size_t i = 0; i != info_.size(); ++i)
                        result += x[info_[i]];

                    return result;
                }

            public:
                /**
                 * \brief Construct a new OneMax_Dummy2 object. Definition refers to
                 *https://doi.org/10.1016/j.asoc.2019.106027
                 *
                 * \param instance The instance number of a problem, which controls the transformation
                 * performed on the original problem.
                 * \param n_variables The dimensionality of the problem to created, 4 by default.
                 **/
                OneMaxDummy2(const int instance, const int n_variables) :
                    PBOProblem(5, instance, n_variables, "OneMaxDummy2"),
                    info_(utils::dummy(n_variables, 0.9, 10000))
                {
                    optimum_.x = std::vector<int>(n_variables,1);
                    optimum_.y = evaluate(optimum_.x);
                    optimum_.x = reset_transform_variables(optimum_.x);
                    optimum_.y = transform_objectives(optimum_.y);
                }
            };
        } // namespace pbo
    } // namespace problem
} // namespace ioh
