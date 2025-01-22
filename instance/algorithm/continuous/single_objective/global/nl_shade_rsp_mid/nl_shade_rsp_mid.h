/********* Begin Register Information **********
{
	"name": "NL_SHADE_RSP_MID",
	"identifier": "NL_SHADE_RSP_MID",
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

#ifndef OFEC_NL_SHADE_RSP_MID_H
#define OFEC_NL_SHADE_RSP_MID_H



#include "../../../../../../core/algorithm/algorithm.h"


namespace ofec {


    class NL_SHADE_RSP_MID : virtual public Algorithm {
        OFEC_CONCRETE_INSTANCE(NL_SHADE_RSP_MID)

    protected:


        void addInputParameters() {};
        void initialize_(Environment* env) override;
        void run_(Environment* env) override;

    protected:


        const int NMBR_OF_RUNS = 30;
        const int NUM_OF_DUMPS = 16;
        int POP_SIZE = 0;
        double fopt = 0;
        int oldPopSize = POP_SIZE;

    };


}

#endif