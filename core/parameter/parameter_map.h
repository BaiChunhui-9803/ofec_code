/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Yiya Diao & Junchen Wang & Changhe Li
* Email: diaoyiyacug@gmail.com & wangjunchen.chris@gmail.com & changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*********************************************************************************/

#ifndef OFEC_PARAMETER_MAP_H
#define OFEC_PARAMETER_MAP_H

#include <utility>
#include <map>
#include "parameter_variant.h"
#include "../exception.h"
#include <ostream>

namespace ofec {
	/**
 	* @class ParameterMap
 	* @brief_byBCH A class for managing a collection of parameters with various types.
 	*
 	* This class provides a flexible way to store, retrieve, and manage parameters
 	* of different types using a unified interface. It is designed to support
 	* dynamic parameter handling in evolutionary computation frameworks.
 	*/
	class ParameterMap {
	private:
		std::map<std::string, ParameterVariant> m_map;
		using IteratorType = std::map<std::string, ParameterVariant>::iterator;
		using CIteratorType = std::map<std::string, ParameterVariant>::const_iterator;
		using RIterType = std::map<std::string, ParameterVariant>::reverse_iterator;
		using RCIterType = std::map<std::string, ParameterVariant>::const_reverse_iterator;

	public:

		virtual ~ParameterMap() = default;
		friend bool operator==(const ParameterMap& p1, const ParameterMap& p2);

		/**
 		* @brief_byBCH Overloaded operator[] to access a parameter by key.
 		* @param key The key of the parameter to access.
 		* @return A reference to the ParameterVariant associated with the key.
 		*/
		ParameterVariant& operator[](const std::string &key);
		const ParameterVariant& at(const std::string &key) const;

		/**
 		* @brief_byBCH Check if a parameter with the given key exists.
 		* @param key The key to check.
 		* @return True if the key exists, false otherwise.
 		*/
		bool has(const std::string &key) const;
		size_t erase(const std::string &key);

		/**
 		* @brief_byBCH Clear all parameters from the map.
 		*/
		void clear();

		/**
 		* @brief_byBCH Get an iterator to the beginning of the map.
 		* @return A const iterator to the beginning of the map.
 		*/
		CIteratorType begin() const { return m_map.begin(); }

 		/**
 		* @brief_byBCH Get an iterator to the end of the map.
 		* @return A const iterator to the end of the map.
 		*/
		CIteratorType end() const { return m_map.end(); }
		/**
 		* @brief_byBCH Get a const iterator to the beginning of the map.
 		* @return A const iterator to the beginning of the map.
 		*/
		CIteratorType cbegin() const { return m_map.cbegin(); }

		/**
 		* @brief_byBCH Get a const iterator to the end of the map.
 		* @return A const iterator to the end of the map.
 		*/
		CIteratorType cend() const { return m_map.cend(); }

		/**
 		* @brief_byBCH Get a parameter value by key, with type checking.
 		* @tparam T The expected type of the parameter.
 		* @param key The key of the parameter to retrieve.
 		* @return A const reference to the parameter value of type T.
 		* @throws Exception if the key does not exist or type mismatch occurs.
 		*/
		template <typename T>
		const T& get(const std::string &key) const {
			if (m_map.count(key) > 0) {
				return std::get<T>(m_map.at(key));
			}
			else {
				throw Exception("parameter \"" + key + "\" is not given.");
			}
		}

		/**
 		* @brief_byBCH Get a parameter value by key, with type checking and default value.
 		* @tparam T The expected type of the parameter.
 		* @param key The key of the parameter to retrieve.
 		* @param default_val The default value to return if the key does not exist.
 		* @return A const reference to the parameter value of type T, or the default value.
 		*/
		template <typename T>
		const T& get(const std::string& key, const T &default_val) const {
			if (m_map.count(key) > 0) {
				return std::get<T>(m_map.at(key));
			}
			else {
				return default_val;
			}
		}
	};

	/**
	 * @brief_byBCH Convert a vector to a series of parameters in a ParameterMap.
	 * @tparam T The type of elements in the vector.
	 * @param v The ParameterMap to store the parameters.
	 * @param name The base name for the parameters.
	 * @param vt The vector to convert.
	 */
	template< typename T>
	inline void vecToParam(ParameterMap &v, const std::string &name, const std::vector<T> &vt) {
		v[name + char(2) + "size"] = int(vt.size());
		for (unsigned i(0); i < vt.size(); ++i) {
			v[name + char(2) + std::to_string(i)] = vt[i];
		}
	}

	/**
 	* @brief_byBCH Convert a series of parameters in a ParameterMap to a vector.
 	* @tparam T The type of elements in the resulting vector.
 	* @param v The ParameterMap containing the parameters.
 	* @param name The base name for the parameters.
 	* @return A vector containing the parameter values.
 	*/
	template< typename T>
	inline std::vector<T> paramToVec(const ParameterMap &v, const std::string &name) {
		std::vector<T> vt(v.get<int>(name + char(2) + "size"));
		for (unsigned i(0); i < vt.size(); ++i) {
			vt[i] = v.get<T>(name + char(2) + std::to_string(i));
		}
		return vt;
	}

	/**
 	* @brief_byBCH Convert a vector of pairs to a series of parameters in a ParameterMap.
 	* @tparam T The type of the first element in the pairs.
 	* @tparam K The type of the second element in the pairs.
 	* @param v The ParameterMap to store the parameters.
 	* @param name The base name for the parameters.
 	* @param vptt The vector of pairs to convert.
 	*/
	template< typename T, typename K>
	void vecPairToParam(ParameterMap &v, const std::string &name, const std::vector<std::pair<T, K>> &vptt) {
		v[name + char(2) + "size"] = int(vptt.size());
		for (unsigned i(0); i < vptt.size(); ++i) {
			v[name + char(2) + std::to_string(i) + char(2) + "1"] = vptt[i].first;
			v[name + char(2) + std::to_string(i) + char(2) + "2"] = vptt[i].second;
		}
	}

	/**
 	* @brief_byBCH Convert a series of parameters in a ParameterMap to a vector of pairs.
 	* @tparam T The type of the first element in the resulting pairs.
 	* @tparam K The type of the second element in the resulting pairs.
 	* @param v The ParameterMap containing the parameters.
 	* @param name The base name for the parameters.
 	* @return A vector of pairs containing the parameter values.
 	*/
	template< typename T, typename K>
	std::vector<std::pair<T, K>>  paramToVecPair(const ParameterMap &v, const std::string &name) {
		std::vector<std::pair<T, K>> vptt(v.get<int>(name + char(2) + "size"));
		for (unsigned i(0); i < vptt.size(); ++i) {
			vptt[i].first = v.get<T>(name + char(2) + std::to_string(i) + char(2) + "1");
			vptt[i].second = v.get<K>(name + char(2) + std::to_string(i) + char(2) + "2");
		}
		return vptt;
	}
}
#endif /* OFEC_PARAMETER_MAP_H */
