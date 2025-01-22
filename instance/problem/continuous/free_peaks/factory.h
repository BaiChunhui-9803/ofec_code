#ifndef OFEC_FREE_PEAKS_FACTORY_H
#define OFEC_FREE_PEAKS_FACTORY_H

#include <map>
#include <functional>
#include "../../../../core/problem/problem.h"

namespace ofec {
	namespace free_peaks {
		template<typename Base>
		struct FactoryFP {
			template<typename T>
			struct Register_ {
				Register_(const std::string& key) {
					map_.emplace(
						key,
						[](Problem *pro, const std::string& subspace_name, const ParameterMap& param, Random *rnd) {
						return new T(pro, subspace_name, param, rnd);
					}
					);
				}
			};

			static unsigned size() {
				return map_.size();
			}

			static Base* produce(const std::string& key, Problem *pro, const std::string& subspace_name, const ParameterMap& param, Random *rnd) {
				auto it = map_.find(key);
				if (it == map_.end())
					THROW("the key is not exist!");
				return it->second(pro, subspace_name, param, rnd);
			}

			static Base* produceRandom(unsigned id, Problem *pro, const std::string& subspace_name, const ParameterMap& param, Random *rnd) {
				auto it(map_.begin());
				advance(it, id);
				return it->second(pro, subspace_name, param, rnd);
			}

			static Base* produceRandom(double randnum, Problem *pro, const std::string& subspace_name, const ParameterMap& param, Random *rnd) {
				unsigned id(map_.size() * randnum);
				auto it(map_.begin());
				advance(it, id);
				return it->second(pro, subspace_name, param, rnd);
			}

			static const std::map<std::string, std::function<Base*(Problem*, const std::string&, const ParameterMap&, Random *rnd)>>& get() {
				return map_;
			}
			FactoryFP() = default;
		private:
			FactoryFP& operator=(const FactoryFP&) = delete;
			FactoryFP& operator=(FactoryFP&&) = delete;
			FactoryFP(const FactoryFP&) = delete;
			FactoryFP(FactoryFP&&) = delete;
			static std::map<std::string, std::function<Base*(Problem*, const std::string&, const ParameterMap&, Random*)>>  map_;
		};

		template<typename Base>
		std::map<std::string, std::function<Base*(Problem*, const std::string&, const ParameterMap&, Random*)>> FactoryFP<Base>::map_;

#define REGISTER_FP(Base, Derived, key) FactoryFP<Base>::Register_<Derived> reg_##Derived(key);
	}
}

#endif // !OFEC_FREE_PEAKS_FACTORY_H
