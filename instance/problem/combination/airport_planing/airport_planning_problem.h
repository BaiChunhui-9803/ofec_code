/********* Begin Register Information **********
{
	"name": "ComOP_APP",
	"identifier": "ComOP_APP",
	"problem tags": [ "APP", "ComOP", "DOP" ]
}
*********** End Register Information **********/

#ifndef AIRPORT_PLANNING_H
#define AIRPORT_PLANNING_H

#include "../../../../core/problem/problem.h"
#include "../../../../core/problem/optima.h"
#include "../../../../core/problem/domain.h"
#include "../../../../core/global.h"
#include "../../../../core/problem/solution.h"
#include <ostream>
#include <string>
#include <regex>
#include <iomanip>      
#include <ctime>  
namespace ofec {

#define CAST_APP(pro) dynamic_cast<APP::airport_planning_problem*>(pr))

	namespace APP {

		struct Flight {
			std::string m_FltNum;
			int m_FltId=0;
			std::string m_strDptrDateTime;
			time_t m_DptrDateTime;
			Real m_DptrDateTimeSpan = 0;
			std::string m_DptrStn;
			int m_DptrStnId=0;
			std::string m_strArrvDateTime;
			time_t m_ArrvDateTime;
			Real m_ArrvDateTimeSpan = 0;
			std::string m_ArrvStn;
			int m_ArrvStnId=0;
			std::string m_Comp;
			std::pair<int, int> m_CompNum;
			int m_CompTotNum = 0;
			std::array<int,3> m_totalPos;
			std::array<int, 3> m_startPos;
		};
		struct Member {
			std::string m_EmpNo;
			int m_EmpId = 0;
			bool m_Captain = false;
			bool m_FirstOfficer = false;
			bool m_Deadhead = true;
			std::string m_Base;
			int m_BaseId = 0;
			Real m_DutyCostPerHour = 0;
			Real m_ParingCostPerHour = 0;
			std::array<bool, 3> m_position;

		};


		
		//class APPsol:public VariableBase {
		//public:
		//	std::vector<std::vector<int>> m_flight_membersids;
		//	APPsol() = default;
		//};
		//class APPSolType :public APPsol {
		//public:
		//	std::vector<std::vector<int>> m_flight_membersid;
		//	std::vector<std::vector<int>> m_member_flightseq;
		//};
		//struct MemberState {
		//	int  m_curPositionId = 0;
		//	int m_networkStateId = 0;
		//	Real m_arriveTime = 0;
		//	std::vector<std::pair<int,int>> m_flights;
		//	bool m_active = true;
		//	std::vector<int> m_loop_SP;
		//	std::vector<std::vector<std::pair<int,int>>> m_loopFlights;
		//	std::vector<std::pair<int,int>> m_loop;
		//	int m_id = 0;
		////	std::vector<int> m_loop_position;
		//};


		struct MemberTempState {
			int  m_curPositionId = 0;
			int m_networkStateId = 0;
			Real m_arriveTime = 0;
			bool m_active = true;
			int m_mid;
		};
		struct MemberFlightPlan {
			MemberTempState m_state;
			std::vector<std::vector<std::pair<int, int>>> m_loopFlights;

		};

		struct MemberFlightLoop {
			MemberTempState m_state;
			std::vector<std::pair<int, int>> m_loop;
		};

		struct FlightState {
			std::array<int, 3> m_totalPos;

			int sumTotalPos() const{
				int sum(0);
				for (int idx(0); idx < 3; ++idx) {
					sum += m_totalPos[idx];
				}
				return sum;
			}
			bool feasible()const {
				return m_totalPos[0] + m_totalPos[1] == 0;
			}
		};
		struct AAP_struct {
			std::vector<std::vector<int>> m_flight_membersid;
			std::vector<bool> m_possible_flight;
			std::vector<std::vector<int>> m_member_flightseq;

			void resetInfo();
			void resize(int fNum,int mNum);
		};

		struct SolState {
			
			std::vector<MemberFlightPlan> m_memberState;
			std::vector<FlightState> m_flightState;
			std::vector<std::vector<int>> m_flight_members;
			int m_unfeasibleFlightNum = 0;
			int m_feasibleMemberNum = 0;
			MemberFlightLoop m_curMemberLoop;
		};

		struct node {
			int m_node_id = 0;
			std::vector<int> m_pres;
			std::vector<int> m_nexts;

			void clear() {
				m_node_id = 0;
				m_pres.clear();
				m_nexts.clear();
			}
		};


		class airport_planning_problem :public Problem {
		public:


		protected:

			// 航班连接网络-
			std::vector<node> m_network;
			int m_startNetId = 0;
			int m_endNetId = 0;
			int m_totalEdges = 0;
			void networkInfo();
			void initNetwork();
			bool connectionFeasible(int fromFid, int toFid)const;
			bool connectStart(int Fid)const {
				return m_baseFid.find(m_flights[Fid].m_DptrStnId) != m_baseFid.end();
			}
			bool connectEnd(int Fid)const {
				return m_baseFid.find(m_flights[Fid].m_ArrvStnId) != m_baseFid.end();

			}

			const std::string m_filepath = g_working_dir + "/instance/problem/combination/airport_planing/";
			const std::string& m_fileformat = "%m/%d/%Y-%H:%M";
			// flight line id to FltNum   
			int m_airport_nameNum = 0;
			std::map<std::string, int> m_airport_name2id;
			// real flights and one virtual flights
			std::vector<Flight> m_flights;
			int m_numFlight=0;
			std::vector<Member> m_members;
			int m_numMember = 0;
			std::string m_startTimeInfo;
	  		time_t m_startTime;
			std::vector<int> m_flights_order;
			std::set<int> m_baseFid;
			//std::vector<Flight> m_flightsInfo_order;
			// 
			// 
			// the captian or first officer
			std::pair<std::vector<int>, std::vector<int>> m_C_FO;
			std::pair<std::set<int>, std::set<int>> m_set_C_FO;

			//parameters 同一单位全部转化为：分钟
			const Real m_time_span = 60;
			const Real m_hour_span = 12;
			//航段之间最小连接时间 MinCT = 40分钟
			Real m_MinCT = 40;
		    //一次执勤飞行时长最多不超过 MaxBlk = 600分钟
			Real m_MaxBlk = 600;
		    //执勤时长最多不超过 MaxDP = 720分钟
			Real m_MaxDP = 720;
			//相邻执勤之间的休息时间不少于 MinRest = 660分钟
			Real m_MinRest = 660;
			//每趟航班最多乘机人数 MaxDH = 5
			int m_MaxDH = 5;
			//排班周期单个机组人员任务环总时长不超过 MaxTAFB = 14400分钟
			Real m_MaxTAFB = 14400;
			//连续执勤天数不超过MaxSuccOn = 4天
			Real m_MaxSuccOn = 4* m_hour_span;
			//相邻两个任务环之间至少有MinVacDay = 2天休息
			Real m_MinVacDay = 2 * m_hour_span;

			inline int insertAirportName(const std::string& name) {
				if (m_airport_name2id.find(name) != m_airport_name2id.end()) {
					return m_airport_name2id.at(name);
				}
				else return m_airport_name2id[name] = ++m_airport_nameNum;
			}


			// filter:
			//2.每个机组人员的下一航段的起飞机场必须和上一航段的到达机场一致（执勤时间内，每个飞行员的上一个航班的ArrvStn和下一个航班的DptrStn相同）；
			//每个机组人员的下一航段的起飞机场的时间大于上一航段的到达机场时间
			void transferSolutions(AAP_struct& cur, const VariableBase& s1)const;
			// 某些飞机不满足飞行条件
			void transferSolutions(SolutionBase& s,  SolState& s1)const;

			
			// void transferSolutions()
			//judge feasible
			//2.每个机组人员的下一航段的起飞机场必须和上一航段的到达机场一致（执勤时间内，每个飞行员的上一个航班的ArrvStn和下一个航班的DptrStn相同）；
			bool judgeArriveDeparture(const AAP_struct& cur)const;
			//4.每趟航班最多乘机人数 MaxDH = 5。
			bool judgeMaxDH(const AAP_struct& cur)const {
				for (auto& it : cur.m_flight_membersid) {
					if (it.size() > m_MaxDH)return false;
				}
				return true;
			}


			// objective
			//:range[0,1], maximum:尽可能多的航班满足机组配置（C1F1）; 
			Real numSatisfiedFlights(const AAP_struct & sol)   const;
			int  numSatisfiedFlightsInt(const AAP_struct& sol) const;
			//:range[0,1], minimum: 尽可能少的总体乘机次数(少摆渡)； 
			Real numFerryFlights(const AAP_struct& sol)const;
			int numFerryFlightsInt(const AAP_struct& sol)const;
			//:range[0,1], minimum: 尽可能少使用替补资格（有正机长资格就尽可能当正机长）。
			Real numSubstituteNum(const AAP_struct& sol)const;
			int numSubstituteNumInt(const AAP_struct& sol)const;



			
			//hard  constrant
            //minimum 1.每个机组人员初始从基地出发并最终回到基地（任务环开始的机场和结束机场相同，即相同Base）；
			Real numNotReturnMem(const AAP_struct& sol)const;
			//minimum 3.每个机组人员相邻两个航段（相邻的两次起飞任务）之间的连接时间不小于 MinCT分钟（MinCT=40分钟）;
			Real totalLessTimeConnectTrail(const AAP_struct& sol)const;

			//soft constrains

			// construct a solution by 
			//SolState 
			void initSolState(SolState& cur)const;
			void addStepSolState(SolState& cur,int next,int pos)const;
			void delStepSolState(SolState& cur)const;
			void initMemState(SolState& cur, int mid)const;
			/*
			* 根据硬约束剪枝
			* 每趟航班最多乘机人数 MaxDH = 5
			每个机组人员初始从基地出发
		    每个机组人员相邻两个航段（相邻的两次起飞任务）之间的连接时间不小于 MinCT分钟（MinCT=40分钟）;
			*/
			void getFeasibleFlightForMember(const SolState& cur, int mid,std::vector<int>& feasible)const;
			//从基地出发并最终回到基地
			bool judgeiFeasibleTaskLoop(const SolState& cur,int mid)const;
			void insertFeasibleTaskLoop( SolState& cur, int mid, int &numFinishedFlight)const;

			bool selectNext(int& next, int& pos, SolState& cur, std::vector<bool>& feasibleNode)const;
			bool dfsSolutions(SolState& cur, std::vector<bool>& feasibleNode, Random *rnd)const;
		public:
			using variable_type = VariableVector<std::vector<int>>;
			using solution_type = Solution<variable_type>;
		
		public:

			//initSolutions with two members
			void initSolutionRandomNetwork(SolutionBase& s, Random *rnd) const;

			void algSolutionRandom2Mem(Random *rnd) const;
			void initSolutionRandom2Mem(SolutionBase& s, Random *rnd) const ;
			void initSolutionRandomRandMem(SolutionBase& s, Random *rnd) const;
			void printSol(std::ostream& out, SolutionBase& s)const;
			void printSolReal(std::ostream& out, SolutionBase& s)const;
			void printSolObj(std::ostream& out, SolutionBase& s)const;
			void calObj(SolutionBase& s)const;
			bool better2Obj(SolutionBase& s1, SolutionBase& s2) const;
			// 通过 problem 继承
			virtual void initializeSolution(SolutionBase& s, Random *rnd) const override;
			virtual bool same(const SolutionBase& s1, const SolutionBase& s2) const override;
			virtual Real variableDistance(const SolutionBase& s1, const SolutionBase& s2) const override;
			virtual Real variableDistance(const VariableBase& s1, const VariableBase& s2) const override;

			void readFile(const std::string& memberFile, const std::string& flightFile);

		protected:
			virtual void initialize_()override;

			// 通过 Problem 继承
			virtual void evaluate_(SolutionBase& s) override;

		protected:
			inline Real diffCurTime(const time_t& t1, const time_t& t2)const {
				return std::difftime(t1, t2) / m_time_span;
			}
		};
		// function 
		extern std::vector<std::string> split(const std::string& input,
			const std::string& regex);
		extern void read_time(std::tm& t, const std::string& time_s, const std::string& format);
		extern void read_time(time_t& t, const std::string& time_s, const std::string& format);
		// for test
		extern void show_time();
	}

	using ComOP_APP = APP::airport_planning_problem;
	}
	

#endif