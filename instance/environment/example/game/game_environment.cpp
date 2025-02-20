/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
*
*********************************************************************************/

// Created: 19 Feb 2025
// Modified:  19 Feb 2025 by BAICHUNHUI

#include "game_environment.h"

namespace ofec {

    void GameEnvironment::addInputParameters() {
        m_input_parameters.add("changeFre", new RangedInt(m_frequency, 500, 10000, 1000));
    }


}