/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Xiaofang Wu
* Email: changhe.lw@google.com Or wuxiaofang@cug.edu.cn
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 11 Mar 2024


#ifndef OFEC_IGD_H
#define OFEC_IGD_H

#include "../functional.h"
#include "../../core/definition.h"
#include <vector>
#include <numeric>

/* return IGD (Inverted Generational Distance) of popularion to the Pareto front
        R. Venkata Rao and Vivek Patel,
        "A multi-objective improved teaching-learning based optimization algorithm for unconstrained and constrained optimization problems,"
        International Journal of Industrial Engineering Computations 5 (2014) 1:C22
        */

namespace ofec {

    template<typename TPop1, typename TPop2>
    Real IGD(const TPop1& ref_pop, const TPop2& pop){
        Real distance = 0;
        for (size_t i = 0; i < ref_pop.size();++i) {
            Real min_d = std::numeric_limits<Real>::max();
            for (int j = 0; j < pop.size(); ++j) {
                Real d = pop[j].objectiveDistance(ref_pop[i].objective());
                if (d < min_d)  
                    min_d = d;
            }
            distance += min_d;
        }
        return distance / ref_pop.size();
    }

    template<typename T>
    Real IGD(const std::vector<std::vector<T>>& ref_points,const std::vector<std::vector<T>>& points){
        Real distance = 0;
        for (size_t i = 0; i < ref_points.size(); ++i) {
            Real min_d = std::numeric_limits<Real>::max();
            for (int j = 0; j < points.size(); ++j) {
                Real d = euclideanDistance(points[j].begin(), points[j].end(), ref_points[i].begin());
                if (d < min_d)  min_d = d;
            }
            distance += min_d;
        }
        return distance / ref_points.size();
    }

    //Reference:
    /*H.Ishibuchi, H.Masuda, Y.Tanigaki, and Y.Nojima.Modified distance
    calculation in generational distance and inverted generational distance,
    Proceedings of the International Conference on Evolutionary
    Multi - Criterion Optimization, 2015, 110 - 125.*/
    template<typename TPop1, typename TPop2>
    Real IGDp(const TPop1& ref_pop, const TPop2& pop, const std::vector<OptimizeMode>& mode) {
        Real distance = 0;
        for (size_t i = 0; i < ref_pop.size(); ++i) {
            Real min_d = std::numeric_limits<Real>::max();
            auto ref_obj = ref_pop[i].objective();
            for (int j = 0; j < pop.size(); ++j) {
                auto pop_obj = pop[j].objective();
                Real d = 0.;
                for (size_t k = 0; k < pop_obj.size(); ++k) {
                    Real temp;
                    if (mode[k] == OptimizeMode::kMinimize) {
                        if (pop_obj[k] > ref_obj[k]) {
                            temp = pop_obj[k] - ref_obj[k];
                        }
                        else {
                            temp = 0.;
                        }
                    }
                    else if (mode[k] == OptimizeMode::kMaximize) {
                        if (pop_obj[k] > ref_obj[k]) {
                            temp = 0.;
                        }
                        else {
                            temp = ref_obj[k] - pop_obj[k];
                        }
                    }
                    d += (temp * temp);
                }
                d = std::sqrt(d);
                if (d < min_d)
                    min_d = d;
            }
            distance += min_d;
        }
        return distance / ref_pop.size();
    }

    template<typename T1,typename T2>
    Real IGDp(const std::vector<std::vector<T1>>& ref_points, const std::vector<std::vector<T2>>& points, const std::vector<OptimizeMode>& mode){
        Real distance = 0;
        for (size_t i = 0; i < ref_points.size(); ++i) {
            Real min_d = std::numeric_limits<Real>::max();
            auto ref_obj = ref_points[i];
            for (int j = 0; j < points.size(); ++j) {
                auto pop_obj = points[j];
                Real d = 0.;
                for (size_t k = 0; k < pop_obj.size(); ++k) {
                    Real temp;
                    if (mode[k] == OptimizeMode::kMinimize) {
                        if (pop_obj[k] > ref_obj[k]) {
                            temp = pop_obj[k] - ref_obj[k];
                        }
                        else {
                            temp = 0.;
                        }
                    }
                    else if (mode[k] == OptimizeMode::kMaximize) {
                        if (pop_obj[k] > ref_obj[k]) {
                            temp = 0.;
                        }
                        else {
                            temp = ref_obj[k] - pop_obj[k];
                        }
                    }
                    d += (temp * temp);
                }
                d = std::sqrt(d);
                if (d < min_d)
                    min_d = d;
            }
            distance += min_d;
        }
        return distance / ref_points.size();
    }

}

#endif
