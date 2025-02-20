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

#ifndef GAME_ENVIRONMENT_H
#define GAME_ENVIRONMENT_H

#include "../../../../core/environment/environment.h"

namespace ofec {
#define CAST_GAMEENVIRONMENT(env) dynamic_cast<GameEnvironment*>(env)

    class GameEnvironment : virtual public Environment {
        OFEC_CONCRETE_INSTANCE(GameEnvironment)
    protected:
        int m_frequency = 0;

        void addInputParameters();
    };
}

#endif //GAME_ENVIRONMENT_H
