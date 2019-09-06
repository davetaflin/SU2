 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
 #if defined TECIOMPI
#include "mpi.h"
 #endif
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
#include "AuxData_s.h"
#include "basicTypes.h"
#include "CodeContract.h"
#include "FieldData_s.h"
#include "fileio.h"
#include "GhostInfo_s.h"
#include "IJK.h"
#include "NodeMap_s.h"
#include "NodeToElemMap_s.h"
namespace tecplot { namespace tecioszl { class ___3970; struct Zone_s { typedef boost::shared_ptr<Zone_s> Ptr; typedef std::map<___3933::___4636, Ptr> ZoneMap; typedef std::pair<___1172, ___2227> ___4607; struct ___1275 { ___372 ___2488; std::vector<___4607> ___2678; ___1275() {} void writeToFile(tecplot::___3933::FileWriterInterface& outputFile, bool ___4480) const { writeScalar(outputFile, ___2488, ___4480); writeScalar(outputFile, (uint64_t)___2678.size(), ___4480); BOOST_FOREACH(___4607 const& zoneCell, ___2678) { writeScalar(outputFile, zoneCell.first, ___4480); writeScalar(outputFile, zoneCell.second, ___4480); } } uint64_t sizeInFile(bool ___4480) const { uint64_t sizeInFile = 0; sizeInFile += scalarSizeInFile(___2488, ___4480); sizeInFile += scalarSizeInFile((uint64_t)___2678.size(), ___4480); BOOST_FOREACH(___4607 const& zoneCell, ___2678) { sizeInFile += scalarSizeInFile(zoneCell.first, ___4480); sizeInFile += scalarSizeInFile(zoneCell.second, ___4480); } return sizeInFile; } ___1275(tecplot::___3933::___1399& inputFile, bool readASCII) { readScalar(inputFile, ___2488, readASCII); uint64_t length; readScalar(inputFile, length, readASCII); ___2678.resize((size_t)length); for(uint64_t i = 0; i < length; ++i) { readScalar(inputFile, ___2678[i].first, readASCII); readScalar(inputFile, ___2678[i].second, readASCII); } } }; typedef std::pair<int32_t, int64_t> ___458; typedef std::map<___458, ___1275> ___1276; std::string ___2683; ZoneType_e ___2684; tecplot::___3933::___1844 m_partitionOffset; tecplot::___3933::___1844 ___2682; double ___2621; int32_t ___2622; ___1172 ___2614; int64_t ___2503; FaceNeighborMode_e ___2458; int64_t ___2651; int64_t ___2501; int64_t ___2650; std::vector<FieldDataType_e> ___2460; std::vector<int> m_passiveVars; std::vector<ValueLocation_e> ___2670; std::vector<___1172> m_shareVarFromZone; bool m_allVarsAreShared; ___1172 m_shareConnectivityFromZone; GhostInfo_s m_ghostNodeInfo; GhostInfo_s m_ghostCellInfo; std::vector<boost::shared_ptr<___1362> > ___2496; std::vector<boost::shared_ptr<___1362> > ___2400; size_t ___2397; boost::shared_ptr<___2730> ___2497; boost::shared_ptr<___2743> m_nodeToElemMap; ___1276 ___2457; boost::shared_ptr<AuxData_s> ___2345; ZoneMap m_partitionMap; std::vector<int> m_partitionOwners; std::vector<int64_t> m_oldNodeNumbers; Zone_s( ___3970 const& tecioData, std::string const& ___4690, ZoneType_e ___4692, int64_t iMin, int64_t jMin, int64_t kMin, int64_t ___1909, int64_t ___2116, int64_t ___2161, double ___3640, int32_t ___3785, ___1172 ___2974, int64_t ___2802, FaceNeighborMode_e ___1284, int64_t ___4192,
int64_t ___2786, int64_t ___4188, std::vector<FieldDataType_e> const& ___1372, std::vector<int> const& passiveVarVector, std::vector<ValueLocation_e> const& ___4326, std::vector<___1172> const& shareVarFromZoneVector, ___1172 ___3549); Zone_s( ___3970 const& tecioData, Zone_s const* partitionParent, int64_t iMin, int64_t jMin, int64_t kMin, int64_t ___1909, int64_t ___2116, int64_t ___2161); void ___965(tecplot::___3933::___4352 ___4368); void deriveCCValues(tecplot::___3933::___4352 ___4368); void filterData(
 #if defined TECIOMPI
MPI_Comm comm, int mainProcess,
 #endif
int32_t ___4336, ___2479 const& minMax); void writeToFile(tecplot::___3933::FileWriterInterface& outputFile, bool ___4480) const; uint64_t sizeInFile(bool ___4480) const; static boost::shared_ptr<Zone_s> makePtr(tecplot::___3933::___1399& inputFile, bool readASCII); private: Zone_s(); }; }}
