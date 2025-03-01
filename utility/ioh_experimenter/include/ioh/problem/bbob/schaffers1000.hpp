#pragma once

#include "schaffers10.hpp"

namespace ioh::problem::bbob
{
    //! Shaffers 1000 problem id 18
    template<typename P=BBOB>
    class Schaffers1000 final: public Schaffers<P>, BBOProblem<Schaffers1000>
    {
    public:
        /**
         * @brief Construct a new Schaffers 1 0 0 0 object
         * 
         * @param instance instance id
         * @param n_variables the dimension of the problem 
         */
        Schaffers1000(const int instance, const int n_variables) :
            Schaffers<P>(18, instance, n_variables,  "Schaffers1000", 1000.0)
        {
        }
    };
    template class Schaffers1000<BBOB>;
}
