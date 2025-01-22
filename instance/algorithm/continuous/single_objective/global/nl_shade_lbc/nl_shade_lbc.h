/********* Begin Register Information **********
{
    "name": "NL_SHADE_LBC",
    "identifier": "NL_SHADE_LBC",
    "tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------/
*
*********************************************************************************/

#ifndef OFEC_NL_SHADE_LBC_H
#define OFEC_NL_SHADE_LBC_H

#include "../../../../../../core/algorithm/algorithm.h"


namespace ofec {


    class NL_SHADE_LBC : virtual public Algorithm {
        OFEC_CONCRETE_INSTANCE(NL_SHADE_LBC)
    protected:

        void addInputParameters() {};
        void initialize_(Environment* env) override;
        void run_(Environment* env) override;

    protected:


    };


}

#endif