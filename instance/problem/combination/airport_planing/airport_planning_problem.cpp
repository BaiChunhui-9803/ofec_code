#include "airport_planning_problem.h"
#include<fstream>
#include<vector>
using namespace ofec;
using namespace APP;
using namespace std;

inline void AAP_struct::resetInfo() {
	for (auto& it : m_flight_membersid) it.clear();
	for (auto& it : m_possible_flight) it = false;
	for (auto& it : m_member_flightseq) it.clear();
}
void AAP_struct::resize(int fNum, int mNum) {
	m_flight_membersid.resize(fNum);
	m_possible_flight.resize(fNum);
	m_member_flightseq.resize(mNum);
}


void airport_planning_problem::networkInfo() {
	m_totalEdges = 0;
	for (auto& it : m_network) {
		m_totalEdges += it.m_nexts.size();
	}
}


void airport_planning_problem::initNetwork(){
	m_network.resize(int(m_numFlight + 1));
	for (auto& it : m_network) it.clear();
	for (int idx(0); idx < m_network.size(); ++idx) {
		m_network[idx].m_node_id = idx;
	}
	m_startNetId = m_numFlight;
	m_endNetId = -1;

//	m_flights.push_back(Flight());
//	m_flights.back()

	for (int fid(0); fid < m_numFlight; ++fid) {
		if (connectStart(fid)) {
			m_network[m_startNetId].m_nexts.push_back(fid);
			m_network[fid].m_pres.push_back(m_startNetId);
	    }
		if (connectEnd(fid)) {
			m_network[fid].m_nexts.push_back(m_endNetId);
		}
		for (int fto(0); fto < m_numFlight; ++fto) {
			if (connectionFeasible(fid, fto)) {
			  
				m_network[fid].m_nexts.push_back(fto);
				m_network[fto].m_pres.push_back(fid);
			}
		}
	}
}

bool ofec::APP::airport_planning_problem::connectionFeasible(int fromFid, int toFid)const
{
    
	return (m_flights[fromFid].m_ArrvStnId == m_flights[toFid].m_DptrStnId)
		&& (int(m_flights[fromFid].m_ArrvDateTimeSpan)+int(m_MinCT)<= int(m_flights[toFid].m_DptrDateTimeSpan));
}



void airport_planning_problem::algSolutionRandom2Mem(Random *rnd) const {
	solution_type s;
	auto& x = s.variable().vect();
	x.resize(m_numFlight);
	for (int fId(0); fId < m_numFlight; ++fId) {
		auto& memberSeq(x[fId]);
		memberSeq.resize(m_flights[fId].m_CompTotNum);
	}
}
void airport_planning_problem::initSolutionRandomRandMem(SolutionBase& s, Random *rnd) const {
	auto& x = dynamic_cast<solution_type&>(s).variable().vect();
	x.resize(m_numFlight);
	for (int fId(0); fId < m_numFlight; ++fId) {
		auto& memberSeq(x[fId]);
		memberSeq.resize(
			rnd->uniform.nextNonStd<int>(m_flights[fId].m_CompTotNum, m_MaxDH));
		int captianNum(m_C_FO.first.size()), fONum(m_C_FO.second.size());
		for (int idx(0); idx < m_flights[fId].m_CompNum.first; ++idx) {
			int captianIdx =rnd->uniform.nextNonStd<int>(0, captianNum);
			int add(0);
			int captianId = m_C_FO.first[captianIdx + add];
			for (int jdx(0); jdx < idx; ++jdx) {
				if (captianId >= memberSeq[jdx]) ++add;
				captianId = m_C_FO.first[captianIdx + add];
			}
			memberSeq[idx] = captianId;
			if (m_members[memberSeq[idx]].m_FirstOfficer) --fONum;
			--captianNum;
		}
		for (int idx(m_flights[fId].m_CompNum.first); idx < m_flights[fId].m_CompTotNum; ++idx) {
			int captianIdx = rnd->uniform.nextNonStd<int>(0, fONum);
			int add(0);

			int captianId = m_C_FO.second[captianIdx + add];
			for (int jdx(0); jdx < idx; ++jdx) {
				if (captianId >= memberSeq[jdx] && m_members[memberSeq[jdx]].m_FirstOfficer) ++add;
				captianId = m_C_FO.second[captianIdx + add];
			}
			memberSeq[idx] = captianId ;
			if (m_members[memberSeq[idx]].m_Captain) --captianNum;
			--fONum;
		}
		for (int idx(m_flights[fId].m_CompTotNum); idx < memberSeq.size(); ++idx) {
			int captianIdx = rnd->uniform.nextNonStd<int>(0, fONum + captianNum);
			if (captianIdx < captianNum) {
				int add(0);
				int captianId = m_C_FO.first[captianIdx + add];
				for (int jdx(0); jdx < idx; ++jdx) {
					if (captianId >= memberSeq[jdx] && m_members[memberSeq[jdx]].m_Captain) ++add;
					captianId = m_C_FO.first[captianIdx + add];
				}
				memberSeq[idx] = captianId;
			}
			else {
				captianIdx -= captianNum;
				int add(0);

				int captianId = m_C_FO.second[captianIdx + add];
				for (int jdx(0); jdx < idx; ++jdx) {
					if (captianId >= memberSeq[jdx] && m_members[memberSeq[jdx]].m_FirstOfficer) ++add;
					captianId = m_C_FO.second[captianIdx + add];
				}
				memberSeq[idx] = captianId;
			}
			
			if (m_members[memberSeq[idx]].m_Captain) --captianNum;
			if (m_members[memberSeq[idx]].m_FirstOfficer) --fONum;
		}


	}
}
void airport_planning_problem::initMemState(SolState& curSol, int mid)const {
	auto& cur = curSol.m_memberState[mid];
	auto& curState = cur.m_state;

	curState.m_curPositionId = m_members[mid].m_BaseId;
	curState.m_networkStateId = m_startNetId;
	curState.m_arriveTime = 0;
	curState.m_active = true;
	curState.m_mid = mid;
	cur.m_loopFlights.clear();
}

void airport_planning_problem::initSolState(SolState& cur)const {
	cur.m_memberState.resize(m_numMember);
	for (int mid(0); mid < cur.m_memberState.size(); mid++) {
		initMemState(cur, mid);
	}
	cur.m_flightState.resize(m_numFlight);
	for (int fid(0); fid < cur.m_flightState.size(); fid++) {
		cur.m_flightState[fid].m_totalPos = m_flights[fid].m_totalPos;
	}

	cur.m_flight_members.resize(m_numFlight);
	for (auto& it : cur.m_flight_members) it.clear();
	cur.m_unfeasibleFlightNum = m_numFlight;
	cur.m_feasibleMemberNum = m_numMember;
}

void airport_planning_problem::addStepSolState(SolState& curSol,int next,int pos)const {
     auto& cur(curSol.m_curMemberLoop);
	auto& curState(cur.m_state);
	curState.m_networkStateId= next;
	curState.m_curPositionId = m_flights[next].m_ArrvStnId;
	curState.m_arriveTime = m_flights[next].m_ArrvDateTimeSpan;
	cur.m_loop.push_back(make_pair(next, pos));
}

void airport_planning_problem::delStepSolState(SolState& curSol)const {

    
	auto& cur(curSol.m_curMemberLoop);
	auto& curState(cur.m_state);
	int mid = curState.m_mid;
	cur.m_loop.pop_back();
	if (cur.m_loop.empty()) {
		curState = curSol.m_memberState[mid].m_state;
	}
	else {
		int next = cur.m_loop.back().first;
		curState.m_networkStateId = next;
		curState.m_curPositionId = m_flights[next].m_ArrvStnId;
		curState.m_arriveTime = m_flights[next].m_ArrvDateTimeSpan;
	}
}

void airport_planning_problem::initSolutionRandomNetwork(SolutionBase& s, Random *rnd) const {
	SolState cur_sols;
	initSolState(cur_sols);
	std::vector<int> memberSeq (m_numMember);
	for (int mid(0); mid < memberSeq.size(); ++mid) {
		memberSeq[mid] = mid;
	}
	rnd->uniform.shuffle(memberSeq.begin(), memberSeq.end());
	int memberSeqId(0);
	int memberId = 0;

	int totalLoop(2000);
	int sucLoopNum(0);
	while (cur_sols.m_unfeasibleFlightNum && cur_sols.m_feasibleMemberNum&&totalLoop--) {
		//initMemState(cur_sols.m_curMemberState);
		while (!cur_sols.m_memberState[memberSeq[memberSeqId]].m_state.m_active) {
			++memberSeqId;
		}
		memberId = memberSeq[memberSeqId];
		cur_sols.m_curMemberLoop.m_state = cur_sols.m_memberState[memberId].m_state;
		 // construct a solutions
		{
			std::vector<int> feasible;
			while (true) {
				int next_step(-1);
				getFeasibleFlightForMember(cur_sols, memberId,feasible);
				if (!feasible.empty()) {
					next_step= feasible[rnd->uniform.nextNonStd<int>(0,feasible.size())];
				}
				if (next_step == -1)
				{
					break;
				}
				// select position 
				int pos_idx(0);
				{
					int total(0);
					int totalCO(0);
					std::array<bool, 3> fff;
					fff.fill(false);
					for (int idx = 0; idx < 3; ++idx) {
						if (cur_sols.m_flightState[next_step].m_totalPos[idx] 
						&& m_members[memberId].m_position[idx]) {
							++total;
							fff[idx] = true;
							if (idx <2) ++totalCO;
							break;
						}
					}
					if (totalCO) {
						for (int idx = 0; idx < 2; ++idx) {
							if (fff[idx]) {
								--total;
								if (total == 0) {
									pos_idx = idx;
									break;
								}
							}
						}
					}
					else {
						for (int idx = 0; idx < 3; ++idx) {
							if (fff[idx]) {
								--total;
								if (total == 0) {
									pos_idx = idx;
									break;
								}
							}
						}
					}

				}
				addStepSolState(cur_sols, next_step, pos_idx);
			}
		}
		if (judgeiFeasibleTaskLoop(cur_sols, memberId)) {
			int numFinishedFlight(0);
			insertFeasibleTaskLoop(cur_sols, memberId, numFinishedFlight);
			cur_sols.m_unfeasibleFlightNum -= numFinishedFlight;
			++sucLoopNum;
		}
		//if (!cur_sols.m_memberState[memberId].m_loopFlights.empty()&&m_network[cur_sols.m_memberState[memberId].m_loopFlights.back().first]
		//	.m_nexts.size() == 1) {
		//	cur_sols.m_memberState[memberId].m_state.m_active = false;
		//	--cur_sols.m_feasibleMemberNum;
		//}
		++memberSeqId;
		if (memberSeqId = memberSeq.size()) {
			memberSeqId = 0;
			rnd->uniform.shuffle(memberSeq.begin(), memberSeq.end());
		}
	}
	transferSolutions(s, cur_sols);


}
void airport_planning_problem::getFeasibleFlightForMember(const SolState& cur, int mid,std::vector<int>& feasible)const {
	feasible.clear();
	int departId(cur.m_curMemberLoop.m_state.m_curPositionId);

	for (auto& nextId : m_network[cur.m_curMemberLoop.m_state.m_networkStateId].m_nexts)
	{
		
		if (nextId >= 0 
		&& m_flights[nextId].m_DptrStnId == departId
		&& int(m_flights[nextId].m_DptrDateTimeSpan) >= (cur.m_curMemberLoop.m_state.m_arriveTime + m_MinCT)
		) {
			for (int idx = 0; idx < 3; ++idx) {
				if (cur.m_flightState[nextId].m_totalPos[idx] && m_members[mid].m_position[idx]) {
					feasible.push_back(nextId);
					break;
				}
			}
		}
	}
}


bool airport_planning_problem::judgeiFeasibleTaskLoop(const SolState& cur,int mid)const {
	auto& curState(cur.m_curMemberLoop);
	return (!curState.m_loop.empty()
	&& m_members[mid].m_BaseId == m_flights[curState.m_loop.back().first].m_ArrvStnId);
}


void airport_planning_problem::insertFeasibleTaskLoop( SolState& cur, int mid, int& numFinishedFlight)const {

	numFinishedFlight = 0;
	for (const auto& fid : cur.m_curMemberLoop.m_loop) {
		--cur.m_flightState[fid.first].m_totalPos[fid.second];
		const auto& numPos(cur.m_flightState[fid.first].m_totalPos);
		if (numPos[0] + numPos[1] == 0) {
			++numFinishedFlight;
		}
	}
	cur.m_memberState[mid].m_loopFlights.emplace_back(cur.m_curMemberLoop.m_loop);
}



bool ofec::APP::airport_planning_problem::dfsSolutions(SolState& cur, std::vector<bool>& feasibleNode, Random *rnd) const
{
	return false;
}


//
//int ofec::APP::airport_planning_problem::selectNext(int&next, int&pos, SolState& cur, std::vector<bool>& feasibleNode)const
//{
//	int next_step(-1);
//	getFeasibleFlightForMember(cur_sols, memberId, feasible);
//	if (!feasible.empty()) {
//		next_step = feasible[rnd->uniform.nextNonStd<int>(0, feasible.size())];
//	}
//	if (next_step == -1)
//	{
//		break;
//	}
//	// select position 
//	int pos_idx(0);
//	{
//		int total(0);
//		int totalCO(0);
//		std::array<bool, 3> fff;
//		fff.fill(false);
//		for (int idx = 0; idx < 3; ++idx) {
//			if (cur_sols.m_flightState[next_step].m_totalPos[idx]
//				&& m_members[memberId].m_position[idx]) {
//				++total;
//				fff[idx] = true;
//				if (idx < 2) ++totalCO;
//				break;
//			}
//		}
//		if (totalCO) {
//			for (int idx = 0; idx < 2; ++idx) {
//				if (fff[idx]) {
//					--total;
//					if (total == 0) {
//						pos_idx = idx;
//						break;
//					}
//				}
//			}
//		}
//		else {
//			for (int idx = 0; idx < 3; ++idx) {
//				if (fff[idx]) {
//					--total;
//					if (total == 0) {
//						pos_idx = idx;
//						break;
//					}
//				}
//			}
//		}
//	}
//	return 0;
//}



void airport_planning_problem::initializeSolution(SolutionBase& s, Random *rnd) const
{
     
    //std::vector<int> cur
}

void ofec::APP::airport_planning_problem::printSol(ostream& out, SolutionBase& s) const
{
	auto& x = dynamic_cast<solution_type&>(s).variable().vect();
	for (int fid(0); fid < x.size(); ++fid) {
		out << fid << ",\t";
		for (auto& mid : x[fid]) {
			out << mid << ",\t";
		}
		out << endl;
	}
	out << endl;
}



void ofec::APP::airport_planning_problem::printSolReal(ostream& out, SolutionBase& s) const
{
	auto& sx = dynamic_cast<solution_type&>(s).variable();
	AAP_struct cur;
	transferSolutions(cur, sx);
	auto& x = cur.m_flight_membersid;
	out << "航班 id,\t" << "旅客　ｉｄ" << std::endl;
	for (int fid(0); fid < x.size(); ++fid) {
		out << fid << ",\t";
		for (auto& mid : x[fid]) {
			out << mid << ",\t";
		}
		out << endl;
	}
	out << endl;
}


void ofec::APP::airport_planning_problem::printSolObj(ostream& out, SolutionBase& s) const
{
	auto& x = dynamic_cast<solution_type&>(s).variable();
	AAP_struct cur;
	transferSolutions(cur,x);


	std::string flag = judgeArriveDeparture(cur) == true ? "true" : "false";

    out << "每个机组人员的下一航段的起飞机场必须和上一航段的到达机场一致,\t" << flag << std::endl;
//	flag = judgeMaxDH(cur) == true ? 1 : 0;
	flag = judgeMaxDH(cur) == true ? "true" : "false";
	out << "每趟航班不超过最多乘机人数,\t" << flag << std::endl;
	out << "满足机组配置的航班数,\t" << numSatisfiedFlightsInt(cur) << std::endl;
	out << "总体乘机次数,\t" << numFerryFlightsInt(cur) << std::endl;
	out << "替补资格数,\t" << numSubstituteNumInt(cur) << std::endl;
	out << "每个机组人员初始从基地出发并最终回到基地的比例,\t" << numNotReturnMem(cur) << std::endl;
	out << "每个机组人员相邻两个航段（相邻的两次起飞任务）之间的连接时间不小于 MinCT分钟 比例,\t" << totalLessTimeConnectTrail(cur) << std::endl;



	//for (int fid(0); fid < x.size(); ++fid) {
	//	out << fid << ",\t";
	//	for (auto& mid : x[fid]) {
	//		out << mid << ",\t";
	//	}
	//	out << endl;
	//}
	//out << endl;
}


void airport_planning_problem::calObj(SolutionBase& s)const {
	auto& ss = dynamic_cast<solution_type&>(s);
	auto& x = dynamic_cast<solution_type&>(s).variable();
	AAP_struct cur;
	transferSolutions(cur, x);
	ss.objective().resize(2);
	ss.objective()[1] = 1.0 - numSatisfiedFlights(cur) + numFerryFlights(cur) + numSubstituteNum(cur);
	ss.objective()[0] = totalLessTimeConnectTrail(cur) + numNotReturnMem(cur);
}

bool airport_planning_problem::better2Obj(SolutionBase& s1, SolutionBase& s2) const {
	auto& ss1 = dynamic_cast<solution_type&>(s1);
	auto& ss2 = dynamic_cast<solution_type&>(s2);
	
	if (ss1.objective()[0] == ss2.objective()[0]) {
		return ss1.objective()[1] < ss2.objective()[1];
	}
	else return ss1.objective()[0] < ss2.objective()[0];
}

void ofec::APP::airport_planning_problem::readFile(const std::string& memberFileName, const std::string& flightFileName)
{
	{
		fstream memberFile(m_filepath + memberFileName);
		std::string input;
		bool first = true;
		m_members.clear();
		int m_member_id = 0;
		Member cur_mem;
		while (memberFile) {
			std::getline(memberFile, input);
			std::string comma = ",";
			auto res_str = APP::split(input, comma);
			if (first) {
				first = false;
			    continue;
			}
		
			if (res_str.size() >= 2) {
				cur_mem.m_EmpId = m_member_id;
				cur_mem.m_EmpNo = res_str[0];
				if (res_str[1] == "Y") {
					cur_mem.m_Captain = false;
				}
				else cur_mem.m_Captain = true;

				if (res_str[2] == "Y") {
					cur_mem.m_FirstOfficer = false;
				}
				else cur_mem.m_FirstOfficer = true;

				if (res_str[3] == "Y") {
					cur_mem.m_Deadhead = false;
				}
				else cur_mem.m_Deadhead = true;

				cur_mem.m_Base = res_str[4];
				cur_mem.m_DutyCostPerHour =   stoi(res_str[5]);
				cur_mem.m_ParingCostPerHour = stoi(res_str[6]);
				m_members.emplace_back(cur_mem);
				++m_member_id;
			}
		}

	}

	{
		fstream flightFile(m_filepath + flightFileName);
		std::cout << m_filepath + flightFileName;
		std::string file_input;
		bool first = true;
		m_flights.clear();
		int id = 0;
		Flight cur_flight;
		while (flightFile) {
			std::getline(flightFile, file_input);
			std::string comma = ",";
			auto res_str = APP::split(file_input, comma);
			if (first) {
				first = false;
				continue;
			}
			if (res_str.size() >= 2) {


				cur_flight.m_FltId = id;
				cur_flight.m_FltNum = res_str[0];
	
				{
					auto datetime = res_str[1] + "-" + res_str[2];
					struct std::tm cur_time {};
					read_time(cur_time, datetime, m_fileformat);
					cur_flight.m_DptrDateTime = mktime(&cur_time);
					cur_flight.m_strDptrDateTime = datetime;
				}

				cur_flight.m_DptrStn = res_str[3];

				{
					auto datetime = res_str[4] + "-" + res_str[5];
					struct std::tm cur_time {};
					read_time(cur_time, datetime, m_fileformat);
					cur_flight.m_ArrvDateTime = mktime(&cur_time);
					cur_flight.m_strArrvDateTime = datetime;
				}

				cur_flight.m_ArrvStn = res_str[6];


				{
					std::string rex = "[CF]";
					cur_flight.m_Comp = res_str[7];
					auto comp_info = APP::split(res_str[7], rex);
					cur_flight.m_CompNum.first = stoi(comp_info[1]);
					cur_flight.m_CompNum.second = stoi(comp_info[2]);
					cur_flight.m_CompTotNum = cur_flight.m_CompNum.first + cur_flight.m_CompNum.second;
				}
				m_flights.emplace_back(cur_flight);
				++id;
			}
		}
	}
}

void airport_planning_problem::initialize_()
{

	Problem::initialize_();
	m_MinCT = 40;
	m_MaxBlk = 600;
	//执勤时长最多不超过 MaxDP = 720分钟
	m_MaxDP = 720;
	//相邻执勤之间的休息时间不少于 MinRest = 660分钟
	m_MinRest = 660;
	//每趟航班最多乘机人数 MaxDH = 5
	m_MaxDH = 5;
	//排班周期单个机组人员任务环总时长不超过 MaxTAFB = 14400分钟
	m_MaxTAFB = 14400;
	//连续执勤天数不超过MaxSuccOn = 4天
	m_MaxSuccOn = 4 * 12;
	//相邻两个任务环之间至少有MinVacDay = 2天休息
	m_MinVacDay = 2 * 12;


	/// </summary>



	readFile("Data A-Crew.csv", "Data A-Flight.csv");
	m_startTimeInfo = "2021-08-11";
	read_time(m_startTime, m_startTimeInfo, "%Y-%m-%d");


	m_numFlight = m_flights.size();
	m_numMember = m_members.size();



	m_airport_name2id.clear();
	m_airport_nameNum = 0;

	for (auto& it : m_flights) {
		it.m_DptrStnId = insertAirportName(it.m_DptrStn);
		it.m_ArrvStnId = insertAirportName(it.m_ArrvStn);
	}

	for (auto& it : m_members) {
		it.m_BaseId = insertAirportName(it.m_Base);
	}

	for (auto& it : m_members) {
		it.m_position = { it.m_Captain,it.m_FirstOfficer,it.m_Deadhead };
	}

	for (auto& it : m_flights) {
		it.m_DptrDateTimeSpan = diffCurTime(it.m_DptrDateTime, m_startTime);
		it.m_ArrvDateTimeSpan = diffCurTime(it.m_ArrvDateTime, m_startTime);
		it.m_totalPos[0] = it.m_CompNum.first;
		it.m_totalPos[1] = it.m_CompNum.second;
		it.m_totalPos[2] = m_MaxDH -it.m_CompTotNum;
		it.m_startPos.fill(0);
		for (int posId(1); posId < it.m_startPos.size(); ++posId) {
			it.m_startPos[posId] = it.m_startPos[posId - 1] + it.m_totalPos[posId - 1];
		}
	}


	m_flights_order.resize(m_numFlight);
	for (int idx(0); idx < m_numFlight; ++idx) {
		m_flights_order[idx] = idx;
	}
	std::sort(m_flights_order.begin(), m_flights_order.end(), [&](int a, int b) {
		return m_flights[a].m_DptrDateTimeSpan < m_flights[b].m_DptrDateTimeSpan;
	});

	m_baseFid.clear();
	for (auto& it : m_members) {
		m_baseFid.insert(it.m_BaseId);
	}
	

	{
		m_C_FO.first.clear();
		m_C_FO.second.clear();

		m_set_C_FO.first.clear();
		m_set_C_FO.second.clear();
		for (int mid(0); mid < m_numMember; ++mid) {
			if (m_members[mid].m_Captain) {
				m_C_FO.first.push_back(mid);
				m_set_C_FO.first.insert(mid);
			}

			if (m_members[mid].m_FirstOfficer) {
				m_C_FO.second.push_back(mid);
				m_set_C_FO.second.insert(mid);
			}
		}
	}


	initNetwork();
}

void ofec::APP::airport_planning_problem::transferSolutions(AAP_struct& cur, const VariableBase& s1)const
{

	cur.resize(m_numFlight, m_numMember);
	cur.resetInfo();
	auto& x1 = dynamic_cast<const variable_type&>(s1).vect();
	cur.m_flight_membersid = x1;
	std::vector<MemberTempState> curMemberState(m_numMember);
	for (int mid(0); mid < curMemberState.size(); mid++) {
		curMemberState[mid].m_curPositionId = m_members[mid].m_BaseId;
		curMemberState[mid].m_arriveTime  = 0;
	}
	
	for (auto& fId : m_flights_order) {
		bool suc(true);
		for (auto& mid : cur.m_flight_membersid[fId]) {
			if (int(m_flights[fId].m_DptrDateTimeSpan) < int(curMemberState[mid].m_arriveTime)
				|| m_flights[fId].m_DptrStnId != curMemberState[mid].m_curPositionId) {
				suc = false;
				break;
			}
		}
		cur.m_possible_flight[fId] = suc;
		if (suc) {
			for (auto& mid : cur.m_flight_membersid[fId]) {
				curMemberState[mid].m_arriveTime = m_flights[fId].m_ArrvDateTimeSpan;
				curMemberState[mid].m_curPositionId= m_flights[fId].m_ArrvStnId;
				cur.m_member_flightseq[mid].push_back(fId);
			}
		}
		else {
			cur.m_flight_membersid[fId].clear();
		}
	}
}


void ofec::APP::airport_planning_problem::transferSolutions(SolutionBase& s, SolState& cur_sol)const
{
	auto& x = dynamic_cast<solution_type&>(s).variable().vect();
	x.resize(m_numFlight);
	bool notFeasible(true);
	//int totalLoop(0);
	//for (int mid(0); mid < m_numMember; ++mid) {
	//	totalLoop += cur_sol.m_memberState[mid].m_loopFlights.size();
	//}
	while (notFeasible) {
		notFeasible = false;
		for (int mid(0); mid < m_numMember; ++mid) {
			for (auto&loop : cur_sol.m_memberState[mid].m_loopFlights) {
				bool loopFeasible(true);
				for (auto& fid : loop) {
					// 某些飞机不满足飞行条件
					if (!cur_sol.m_flightState[fid.first].feasible()) {
						loopFeasible = false;
					}
				 }
				if (!loopFeasible) {
					notFeasible = true;
					for (auto& fid : loop) {
						++cur_sol.m_flightState[fid.first].m_totalPos[fid.second];
					}
					loop.clear();
					//--totalLoop;
					//std::cout << "totalLoop\t" << totalLoop << std::endl;
				 }

			}
		}
	}


	for (int fid(0); fid < m_numFlight; ++fid) {
		x[fid].resize(m_MaxDH-cur_sol.m_flightState[fid].sumTotalPos());
		for (auto& it2 : x[fid]) it2 = -1;
	}
	for (int mid(0); mid < m_numMember; ++mid) {

		for (auto loop : cur_sol.m_memberState[mid].m_loopFlights) {
			bool loopFeasible(true);
			for (auto& fid : loop) {
				int posId(m_flights[fid.first].m_startPos[fid.second]);
				while (x[fid.first][posId] != -1) {
					++posId;
				}
				x[fid.first][posId] = mid;
			}
		}
	}



}



bool ofec::APP::airport_planning_problem::judgeArriveDeparture(const AAP_struct& cur)const
{
	for (int mid(0); mid < cur.m_member_flightseq.size(); ++mid) {
		for (size_t seqId(1); seqId < cur.m_member_flightseq[mid].size(); ++seqId) {
			if (m_flights[cur.m_member_flightseq[mid][seqId - 1]].m_ArrvStnId !=
				m_flights[cur.m_member_flightseq[mid][seqId]].m_DptrStnId) return false;
		}
	}
	return true;
}


// objective
Real ofec::APP::airport_planning_problem::numSatisfiedFlights(const AAP_struct& sol)const {
	return static_cast<Real>(numSatisfiedFlightsInt(sol)) / sol.m_possible_flight.size();
}

int airport_planning_problem::numSatisfiedFlightsInt(const AAP_struct& sol)const {
	int totalNum(0);
	for (const auto& it : sol.m_possible_flight) {
		if (it) ++totalNum;
	}
	return totalNum;
}
// objective
Real ofec::APP::airport_planning_problem::numFerryFlights(const AAP_struct& sol)const {
	Real totalNum(0);
	for (int fId(0); fId < sol.m_flight_membersid.size(); ++fId) {
		totalNum += 
		static_cast<Real>(std::max<Real>(0,sol.m_flight_membersid[fId].size() - m_flights[fId].m_CompTotNum))
		/ static_cast<Real>(m_MaxDH - m_flights[fId].m_CompTotNum);
	}
	return totalNum / static_cast<Real>(sol.m_possible_flight.size());
}

int airport_planning_problem::numFerryFlightsInt(const AAP_struct& sol)const {
	int totalNum(0);
	for (int fId(0); fId < sol.m_flight_membersid.size(); ++fId) {
		totalNum +=
			static_cast<int>(std::max<int>(0, 
			sol.m_flight_membersid[fId].size() - m_flights[fId].m_CompTotNum));
	}
	return totalNum;
}


Real ofec::APP::airport_planning_problem::numSubstituteNum(const AAP_struct& sol)const
{
	Real totalNum(0);
	int totalFight(0);
	int curNum(0);
	int total = 0;
	for (int fId(0); fId < sol.m_flight_membersid.size(); ++fId) {
		curNum = 0;
		for (int midx(m_flights[fId].m_CompNum.first); midx < sol.m_flight_membersid[fId].size(); ++midx) {
			if (m_members[sol.m_flight_membersid[fId][midx]].m_Captain) ++curNum;
		 }
		total = sol.m_flight_membersid[fId].size();
		if (total) {
			++totalFight;
			totalNum += static_cast<Real>(curNum) / static_cast<Real>(total);
		}

	}
	return totalNum / static_cast<Real>(sol.m_flight_membersid.size());
	//return totalNum / static_cast<Real>(totalFight);
}
int airport_planning_problem::numSubstituteNumInt(const AAP_struct& sol)const {
	int totalNum(0);
	int totalFight(0);
	int curNum(0);
	int total = 0;
	for (int fId(0); fId < sol.m_flight_membersid.size(); ++fId) {
		curNum = 0;
		for (size_t midx(m_flights[fId].m_CompNum.first); midx < sol.m_flight_membersid[fId].size(); ++midx) {
			if (m_members[sol.m_flight_membersid[fId][midx]].m_Captain) ++curNum;
		}
		total = sol.m_flight_membersid[fId].size();
		if (total) {
			totalNum += static_cast<int>(curNum);
		}
	}
	return totalNum;
}
Real ofec::APP::airport_planning_problem::numNotReturnMem(const AAP_struct& sol) const
{
	Real cur(0);
	for (size_t mid(0); mid < sol.m_member_flightseq.size(); ++mid) {
		if (sol.m_member_flightseq[mid].empty()) ;
		else if (m_members[mid].m_BaseId != m_flights[sol.m_member_flightseq[mid].back()].m_ArrvStnId) cur += 1.0;
	}
	return cur/static_cast<Real>(m_numMember);
}


Real ofec::APP::airport_planning_problem::totalLessTimeConnectTrail(const AAP_struct& cur)const {
    
	Real totalTime = m_numFlight*static_cast<Real>(m_MaxDH-1)*m_MinCT;
	Real result(0);
	std::vector<std::vector<Real>> connectTime(m_numMember);
	for (size_t mid(0); mid < cur.m_member_flightseq.size(); ++mid) {
		connectTime[mid].resize(cur.m_member_flightseq[mid].size());
		for (size_t seqId(1); seqId < cur.m_member_flightseq[mid].size(); ++seqId) {
			result += (- std::min<Real>(0,
				m_flights[cur.m_member_flightseq[mid][seqId]].m_DptrDateTimeSpan
				- m_flights[cur.m_member_flightseq[mid][seqId - 1]].m_ArrvDateTimeSpan- m_MinCT));
			connectTime[mid][seqId] = m_flights[cur.m_member_flightseq[mid][seqId]].m_DptrDateTimeSpan
				- m_flights[cur.m_member_flightseq[mid][seqId - 1]].m_ArrvDateTimeSpan;
		}
	}
	return result / totalTime;
}


void ofec::APP::airport_planning_problem::evaluate_(SolutionBase& s)
{
    
}

bool airport_planning_problem::same(const SolutionBase& s1, const SolutionBase& s2) const
{
	return false;
}

Real airport_planning_problem::variableDistance(const SolutionBase& s1, const SolutionBase& s2) const
{
	return Real();
}

Real airport_planning_problem::variableDistance(const VariableBase& s1, const VariableBase& s2) const
{
	return Real();
}




std::vector<std::string> ofec::APP::split(const std::string& input, const std::string& regex)
{
	std::regex re(regex);
	std::sregex_token_iterator first{ input.begin(), input.end(), re, -1 }, last;
	return { first, last };
}

void ofec::APP::read_time(std::tm& t, const std::string& time_s, const std::string& format) {
	std::istringstream ss(time_s);
	//ss.imbue(std::locale("de_DE.utf-8"));
	ss >> std::get_time(&t, format.c_str());

}

void ofec::APP::show_time()
{
	time_t tm1 = std::time(nullptr);
	clock_t cl1 = std::clock();
	char* pc = std::ctime(&tm1);

	std::cout << "time:\t" << tm1 << std::endl;
	std::cout << "clock:\t" << cl1 << std::endl;
	std::cout << "ctime:\t" << pc << std::endl;

	char buffer[80];
	wchar_t wcbuffer[80];
	// 本地时间
	std::cout << "local " << std::endl;
	std::tm* plocal = std::localtime(&tm1);
	time_t mk1 = std::mktime(plocal);
	std::cout << "mktime:\t" << mk1 << std::endl;
	std::strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", plocal);
	std::cout << "strftime:\t" << buffer << std::endl;
	std::wcsftime(wcbuffer, 80, L"%Y-%m-%d %H:%M:%S", plocal);
	std::wcout << L"wcsftime:\t" << wcbuffer << std::endl;
	char* asc = std::asctime(plocal);
	std::cout << "asctime:\t" << asc << std::endl;


	// 通用时间
	std::cout << "Universal Coordinated" << std::endl;
	std::tm* ptime = std::gmtime(&tm1);
	time_t mk2 = std::mktime(ptime);
	std::cout << "mktime:\t" << mk2 << std::endl;
	std::strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", ptime);
	std::cout << "strftime:\t" << buffer << std::endl;
	std::wcsftime(wcbuffer, 80, L"%Y-%m-%d %H:%M:%S", ptime);
	std::wcout << L"wcsftime:\t" << wcbuffer << std::endl;
	char* asctime = std::asctime(ptime);
	std::cout << "asctime:\t" << asctime << std::endl;

	// c++17
	// windows下timespec不在std里
	timespec ts;
	timespec_get(&ts, TIME_UTC);
	std::strftime(buffer, sizeof(buffer), "%D %T", std::gmtime(&ts.tv_sec));
	std::cout << "Current time:\t" << buffer << '.' << ts.tv_nsec << " UTC" << std::endl;

	time_t tm2 = std::time(nullptr);
	double d = std::difftime(tm2, tm1);
	std::cout << "Wall time passed:\t " << d << " s." << std::endl;

}

void ofec::APP::read_time(time_t& t, const std::string& time_s, const std::string& format)
{
	struct std::tm cur_time {};
	read_time(cur_time, time_s, format);
	t = mktime(&cur_time);

}


