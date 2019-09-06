#include "Zone_s.h"
#include "ThirdPartyHeadersBegin.h"
#include <algorithm>
#include <sstream>
#include <boost/make_shared.hpp>
#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "FieldData_s.h"
 #if defined TECIOMPI
#include "MPIUtil.h"
 #else
#include "JobControl_s.h"
 #endif
#include "TecioData.h"
using namespace tecplot::___3933; namespace tecplot { namespace tecioszl { namespace {
 #define MIN_NODES_FOR_MULTITHREAD 100000 
 #define MIN_CELLS_FOR_MULTITHREAD 50000 
void getZoneCounts(int64_t& nodeCount, int64_t& cellCount, ___1844 const& zoneSize, ZoneType_e ___4692) { if (___4692 == ___4704) { nodeCount = (int64_t)zoneSize.i() * zoneSize.___2105() * zoneSize.___2134(); cellCount = (int64_t)zoneSize.i() * zoneSize.___2105() * std::max((int64_t)1, (int64_t)(zoneSize.___2134() - 1)); } else { nodeCount = (int64_t)zoneSize.i(); cellCount = (int64_t)zoneSize.___2105(); } } } Zone_s::Zone_s( ___3970 const& tecioData, std::string const& ___4690, ZoneType_e ___4692, int64_t iMin, int64_t jMin, int64_t kMin, int64_t ___1909, int64_t ___2116, int64_t ___2161, double ___3640, int32_t ___3785, ___1172 ___2974, int64_t ___2802, FaceNeighborMode_e ___1284, int64_t ___4192, int64_t ___2786, int64_t ___4188, std::vector<FieldDataType_e> const& ___1372, std::vector<int> const& passiveVarVector, std::vector<ValueLocation_e> const& ___4326, std::vector<___1172> const& shareVarFromZoneVector, ___1172 ___3549) : ___2683(___4690) , ___2684(___4692) , m_partitionOffset((___81)(iMin - 1), (___81)(jMin - 1), (___81)(kMin - 1)) , ___2682((___81)(___1909 - iMin + 1), (___81)(___2116 - jMin + 1), (___81)(___2161 - kMin + 1)) , ___2621(___3640) , ___2622(___3785) , ___2614(___2974) , ___2503(___2802) , ___2458(___1284) , ___2651(___4192) , ___2501(___2786) , ___2650(___4188) , ___2460(___1372) , m_passiveVars(passiveVarVector) , ___2670(___4326) , m_shareVarFromZone(shareVarFromZoneVector) , m_shareConnectivityFromZone(___3549) , ___2397(static_cast<size_t>(-1)) , ___2345(new AuxData_s) { REQUIRE(___4690.size() > 0); REQUIRE(VALID_ENUM(___4692, ZoneType_e)); REQUIRE(iMin > 0); REQUIRE(jMin > 0); REQUIRE(kMin > 0); REQUIRE(0 < ___1909); REQUIRE(0 < ___2116); REQUIRE(IMPLICATION(___4692 == ___4704 || ___4692 == ___4698 || ___4692 == ___4699, 0 < ___2161)); REQUIRE("solutionTime can be anything"); REQUIRE(0 <= ___3785); REQUIRE(0 <= ___2974); REQUIRE(0 <= ___2802); REQUIRE(___2670.size() > 0); REQUIRE(VALID_ENUM(___1284, FaceNeighborMode_e)); REQUIRE(IMPLICATION(___4692 == ___4699, 0 <= ___4192)); REQUIRE(IMPLICATION(___4692 == ___4699 || ___4692 == ___4698, 0 <= ___2786)); REQUIRE(IMPLICATION(___4692 == ___4699 || ___4692 == ___4698, 0 <= ___4188)); int64_t nodeCount; int64_t cellCount; getZoneCounts(nodeCount, cellCount, ___2682, ___2684); if (___2684 == ___4696) { ___2682.___3497(2); } else if (___2684 == ___4702)
{ ___2682.___3497(3); } else if (___2684 == ___4700 || ___2684 == ___4701) { ___2682.___3497(4); } else if (___2684 == ___4695) { ___2682.___3497(8); } try { ___2496.resize(___2670.size()); } catch (std::bad_alloc const&) { std::cerr << "Out of memory while storing " << ___2670.size() << " zone value locations." << std::endl; throw; } try { ___2400.resize(___2670.size()); } catch (std::bad_alloc const&) { std::cerr << "Out of memory while storing " << ___2670.size() << " zone derived variable structures." << std::endl; throw; } std::set<___3493> zoneSet; try { zoneSet = tecioData.zoneSet(); } catch (std::bad_alloc const&) { std::cerr << "Out of memory while copying set of defined zones." << std::endl; throw; } ___1172 numVarsShared = 0; try { for(size_t ___4291 = 0; ___4291 < ___2670.size(); ++___4291) { if (m_passiveVars[___4291] != 0) { ___2496[___4291] = TypedFieldDataFactory().make(___2460[___4291]); ___2400[___4291] = TypedFieldDataFactory().make(___2460[___4291]); } else if (m_shareVarFromZone[___4291] != 0) { numVarsShared++; if (zoneSet.find(static_cast<___3493>(m_shareVarFromZone[___4291])) == zoneSet.end()) { std::ostringstream ___2892; ___2892 << "Invalid zone specified for variable sharing.\n" "Specified non-existent zone " << m_shareVarFromZone[___4291] << " for variable " << ___4291 + 1; throw ___3970::Error(___2892.str().c_str()); } else { Zone_s* zonePtr = tecioData.zonePtr(m_shareVarFromZone[___4291]); ___478(VALID_REF(zonePtr)); ___2496[___4291] = zonePtr->___2496[___4291]; ___2400[___4291] = zonePtr->___2400[___4291]; } } else { if (___2397 == static_cast<size_t>(-1)) ___2397 = ___4291; ___2496[___4291] = TypedFieldDataFactory().make(___2460[___4291]); ___2400[___4291] = TypedFieldDataFactory().make(___2460[___4291]); if (___2670[___4291] == ___4330) { ___2496[___4291]->___2669 = ___4330; ___2496[___4291]->___2668 = nodeCount; ___2400[___4291]->___2669 = ___4328; ___2400[___4291]->___2668 = cellCount; } else { ___2496[___4291]->___2669 = ___4328; ___2496[___4291]->___2668 = cellCount; ___2400[___4291]->___2669 = ___4330; ___2400[___4291]->___2668 = nodeCount; } } } } catch (std::bad_alloc const&) { std::cerr << "Out of memory while storing " << ___2670.size() << " zone variable structures." << std::endl; throw; } m_allVarsAreShared = numVarsShared == static_cast<___1172>(___2670.size()); if (m_shareConnectivityFromZone) { if (zoneSet.find(static_cast<___3493>(m_shareConnectivityFromZone)) == zoneSet.end()) { try { std::ostringstream ostream; ostream << "Invalid zone specified for nodemap sharing.\n" "Specified non-existent zone " << m_shareConnectivityFromZone;
throw ___3970::Error(ostream.str().c_str()); } catch(std::bad_alloc const&) { std::cerr << "Out of memory while attempting to report:" << std::endl << "Invalid zone specified for nodemap sharing." << std::endl << "Specified non-existent zone " << m_shareConnectivityFromZone << "." << std::endl; throw; } } else { Zone_s* zonePtr = tecioData.zonePtr(m_shareConnectivityFromZone); ___2497 = zonePtr->___2497; } } else if (___2684 == ___4704) { } else { try { ___2497 = ___2730::makePtr(___2682.___1669(), ___2116, ___2682.___1670()); } catch(std::bad_alloc const&) { std::cerr << "Out of memory while creating zone node map." << std::endl; throw; } } } Zone_s::Zone_s( ___3970 const& tecioData, Zone_s const* partitionParent, int64_t iMin, int64_t jMin, int64_t kMin, int64_t ___1909, int64_t ___2116, int64_t ___2161) { REQUIRE(VALID_REF(partitionParent)); REQUIRE(iMin > 0); REQUIRE(jMin > 0); REQUIRE(kMin > 0); REQUIRE(___1909      <= partitionParent->___2682.i()); REQUIRE(___2116 <= partitionParent->___2682.___2105()); REQUIRE(___2161    <= partitionParent->___2682.___2134()); REQUIRE(___1909 > iMin); REQUIRE(___2116 > jMin); REQUIRE(___2161 > kMin); *this = Zone_s(tecioData, partitionParent->___2683, partitionParent->___2684, iMin, jMin, kMin, ___1909, ___2116, ___2161, partitionParent->___2621, partitionParent->___2622, partitionParent->___2614, 0, partitionParent->___2458, 0, 0, 0, partitionParent->___2460, partitionParent->m_passiveVars, partitionParent->___2670, partitionParent->m_shareVarFromZone, partitionParent->m_shareConnectivityFromZone); } namespace { struct FENodalValueDerivationData { ___2718 begin; ___2718 end; ___1362* nativeFieldData; ___1362* derivedFieldData; ___2743* nodeToElemMap; double minVal; double maxVal; FENodalValueDerivationData(___2718 begin, ___2718 end, ___1362* nativeFieldData, ___1362* derivedFieldData, ___2743* nodeToElemMap) : begin(begin) , end(end) , nativeFieldData(nativeFieldData) , derivedFieldData(derivedFieldData) , nodeToElemMap(nodeToElemMap) , minVal(std::numeric_limits<double>::max()) , maxVal(-std::numeric_limits<double>::max()) {} }; inline double getValueFromTypedRawPtr(void* rawPtr, int64_t index, FieldDataType_e type) { REQUIRE(VALID_REF(rawPtr)); REQUIRE(0 <= index); REQUIRE(VALID_ENUM(type, FieldDataType_e)); double ___4314; switch(type) { case FieldDataType_Float: ___4314 = static_cast<double>(((float*)(rawPtr))[index]); break; case FieldDataType_Double: ___4314 = ((double*)(rawPtr))[index]; break; case FieldDataType_Int32: ___4314 = static_cast<double>(((int32_t*)(rawPtr))[index]); break; case FieldDataType_Int16: ___4314 = static_cast<double>(((int16_t*)(rawPtr))[index]); break; case FieldDataType_Byte: ___4314 = static_cast<double>(((uint8_t*)(rawPtr))[index]); break; case ___1365:
___4314 = static_cast<double>((((uint8_t*)(rawPtr))[index / 8] >> (index % 8)) & 1); break; case ___1368: ___4314 = static_cast<double>(((int64_t*)(rawPtr))[index]); break; default: ___478(___1305); ___4314 = 0.0; break; } return ___4314; } inline void setValueToTypedRawPtr(void* rawPtr, int64_t index, FieldDataType_e type, double ___4314) { REQUIRE(VALID_REF(rawPtr)); REQUIRE(0 <= index); REQUIRE(VALID_ENUM(type, FieldDataType_e)); switch(type) { case FieldDataType_Float: ((float*)(rawPtr))[index] = static_cast<float>(___4314); break; case FieldDataType_Double: ((double*)(rawPtr))[index] = ___4314; break; case FieldDataType_Int32: ((int32_t*)(rawPtr))[index] = static_cast<int32_t>(___4314); break; case FieldDataType_Int16: ((int16_t*)(rawPtr))[index] = static_cast<int16_t>(___4314); break; case FieldDataType_Byte: ((uint8_t*)(rawPtr))[index] = static_cast<uint8_t>(___4314); break; case ___1365: if (___4314 < 1.0) { ___4314 = 0.0; ((uint8_t *)(rawPtr))[index / 8] &= ~(static_cast<uint8_t>(01) << (index % 8)); } else { ___4314 = 1.0; ((uint8_t *)(rawPtr))[index / 8] |= static_cast<uint8_t>(01) << (index % 8); } break; case ___1368: ((int64_t*)(rawPtr))[index] = static_cast<int64_t>(___4314); break; default: ___478(___1305); break; } } void deriveRangeOfNodalValues(___90 threadData) { FENodalValueDerivationData* derivationData = reinterpret_cast<FENodalValueDerivationData*>(threadData); REQUIRE(derivationData->derivedFieldData->___2459 == derivationData->nativeFieldData->___2459); void* rawCCPtr = derivationData->nativeFieldData->getRawPointer(); void* rawNodalPtr = derivationData->derivedFieldData->getRawPointer(); FieldDataType_e ___1363 = derivationData->nativeFieldData->___2459; for(___2718 ___2709 = derivationData->begin; ___2709 < derivationData->end; ++___2709) { double ___4314 = 0.0; ___465 startingIndex = derivationData->nodeToElemMap->m_elemIndex[___2709]; ___465 howManyCells = derivationData->nodeToElemMap->cellCountForNode(___2709); for(___465 whichCell = 0; whichCell < howManyCells; ++whichCell) { ___465 ___449 = derivationData->nodeToElemMap->m_elem[startingIndex + whichCell]; ___4314 += getValueFromTypedRawPtr(rawCCPtr, ___449, ___1363); } ___4314 /= howManyCells; setValueToTypedRawPtr(rawNodalPtr, ___2709, ___1363, ___4314); derivationData->minVal = std::min(derivationData->minVal, ___4314); derivationData->maxVal = std::max(derivationData->maxVal, ___4314); } } } void Zone_s::___965(___4352 ___4368) { REQUIRE(0 < ___4368 && ___4368 <= static_cast<___4352>(___2496.size())); REQUIRE(___2496[___4368 - 1]->___2669 == ___4328); ___1362::Ptr nativeFieldData = ___2496[___4368 - 1]; ___1362::Ptr derivedFieldData = ___2400[___4368 - 1]; derivedFieldData->assignValues(derivedFieldData->___2668, 0.0); if (___2684 == ___4704) { boost::scoped_array<uint16_t> divisor(new uint16_t[derivedFieldData->___2668]);
for(size_t i = 0; i < derivedFieldData->___2668; ++i) divisor[i] = 0; int64_t const ___461 = std::max(int64_t(1), (int64_t)(___2682.i() - 1)); int64_t const ___466 = std::max(int64_t(1), (int64_t)(___2682.___2105() - 1)); int64_t const ___467 = std::max(int64_t(1), (int64_t)(___2682.___2134() - 1)); int64_t const nodeIMax = ___2682.i(); int64_t const nodeJMax = ___2682.___2105(); int64_t const nodeKMax = ___2682.___2134(); std::vector<int64_t> nodes; for (int64_t i = 0; i < ___461; ++i) { for (int64_t ___2105 = 0; ___2105 < ___466; ++___2105) { for (int64_t ___2134 = 0; ___2134 < ___467; ++___2134) { int64_t const index = (___2134 * nodeJMax + ___2105) * nodeIMax + i; nodes.resize(0); nodes.push_back(index); if (nodeIMax > 1) { nodes.push_back(index + 1); } if (nodeJMax > 1) { size_t count = nodes.size(); for(size_t n = 0; n < count; ++n) { nodes.push_back(nodes[n] + nodeIMax); } } if (nodeKMax > 1) { size_t count = nodes.size(); for(size_t n = 0; n < count; ++n) { nodes.push_back(nodes[n] + nodeIMax * nodeJMax); } } for(size_t n = 0; n < nodes.size(); ++n) { derivedFieldData->___3504(nodes[n], derivedFieldData->___1780(nodes[n]) + nativeFieldData->___1780(index)); ++divisor[nodes[n]]; } } } } for(size_t n = 0; n < derivedFieldData->storedValueCount(); ++n) derivedFieldData->___3504((___81)n, derivedFieldData->___1780(n) / (double)divisor[n]); } else { ___478(___2682.___1669() > 0); if (!m_nodeToElemMap) m_nodeToElemMap.reset(new ___2743(*___2497, ___2682.___1670()));
 #if !defined TECIOMPI
int numThreads = 1; if (m_nodeToElemMap->m_nodeCount >= MIN_NODES_FOR_MULTITHREAD) numThreads = std::min(___2122::___2827(), static_cast<int>((m_nodeToElemMap->m_nodeCount - 1) / MIN_NODES_FOR_MULTITHREAD + 1)); if (numThreads == 1) {
 #endif
FENodalValueDerivationData derivationData(0, ___2682.___1670(), nativeFieldData.get(), derivedFieldData.get(), m_nodeToElemMap.get()); deriveRangeOfNodalValues((___90)&derivationData); derivedFieldData->___3499(derivationData.minVal, derivationData.maxVal);
 #if !defined TECIOMPI
} else { std::vector<boost::shared_ptr<FENodalValueDerivationData> > nodalDerivationData; for(int i = 0; i < numThreads; ++i) { ___2718 const begin = static_cast<___2718>((size_t)m_nodeToElemMap->m_nodeCount * i / numThreads); ___2718 const end = static_cast<___2718>((size_t)m_nodeToElemMap->m_nodeCount * (i + 1) / numThreads); nodalDerivationData.push_back(boost::make_shared<FENodalValueDerivationData>( begin, end, nativeFieldData.get(), derivedFieldData.get(), m_nodeToElemMap.get())); } ___2122 ___2119; for(int i = 0; i < numThreads; ++i) ___2119.addJob(deriveRangeOfNodalValues, reinterpret_cast<___90>(nodalDerivationData[i].get())); ___2119.wait(); double minVal = std::numeric_limits<double>::max(); double maxVal = -std::numeric_limits<double>::max(); for(int i = 0; i < numThreads; ++i) { minVal = std::min(minVal, nodalDerivationData[i]->minVal); maxVal = std::max(maxVal, nodalDerivationData[i]->maxVal); } derivedFieldData->___3499(minVal, maxVal); }
 #endif
} } namespace { struct FECCValueDerivationData { int64_t begin; int64_t end; ___1362* nativeFieldData; ___1362* derivedFieldData; ___2730* ___2723; double minVal; double maxVal; FECCValueDerivationData(int64_t begin, int64_t end, ___1362* nativeFieldData, ___1362* derivedFieldData, ___2730* ___2723) : begin(begin) , end(end) , nativeFieldData(nativeFieldData) , derivedFieldData(derivedFieldData) , ___2723(___2723) , minVal(std::numeric_limits<double>::max()) , maxVal(-std::numeric_limits<double>::max()) {} }; void deriveRangeOfCCValues(___90 threadData) { FECCValueDerivationData* derivationData = reinterpret_cast<FECCValueDerivationData*>(threadData); REQUIRE(derivationData->derivedFieldData->___2459 == derivationData->nativeFieldData->___2459); void* rawNodalPtr = derivationData->nativeFieldData->getRawPointer(); void* rawCCPtr = derivationData->derivedFieldData->getRawPointer(); FieldDataType_e ___1363 = derivationData->nativeFieldData->___2459; for(int64_t ___449 = derivationData->begin; ___449 < derivationData->end; ++___449) { double ___4314 = 0.0; ___682 ___2789 = static_cast<___682>(derivationData->___2723->___2500); for(___682 ___681 = 0; ___681 < ___2789; ++___681) { int64_t ___2709 = derivationData->___2723->___4314(___449 * ___2789 + ___681); ___4314 += getValueFromTypedRawPtr(rawNodalPtr, ___2709, ___1363); } ___4314 /= ___2789; setValueToTypedRawPtr(rawCCPtr, ___449, ___1363, ___4314); derivationData->minVal = std::min(derivationData->minVal, ___4314); derivationData->maxVal = std::max(derivationData->maxVal, ___4314); } } } void Zone_s::deriveCCValues(___4352 ___4368) { REQUIRE(0 < ___4368 && ___4368 <= static_cast<___4352>(___2496.size())); REQUIRE(___2496[___4368 - 1]->___2669 == ___4330); ___1362::Ptr nativeFieldData = ___2496[___4368 - 1]; ___1362::Ptr derivedFieldData = ___2400[___4368 - 1]; derivedFieldData->assignValues(derivedFieldData->___2668, 0.0); if (___2684 == ___4704) { int64_t const ___461 = std::max(int64_t(1), (int64_t)(___2682.i() - 1)); int64_t const ___466 = std::max(int64_t(1), (int64_t)(___2682.___2105() - 1)); int64_t const ___467 = std::max(int64_t(1), (int64_t)(___2682.___2134() - 1)); int64_t const nodeIMax = ___2682.i(); int64_t const nodeJMax = ___2682.___2105(); int64_t const nodeKMax = ___2682.___2134(); std::vector<int64_t> nodes; for (int64_t i = 0; i < ___461; ++i) { for (int64_t ___2105 = 0; ___2105 < ___466; ++___2105) { for (int64_t ___2134 = 0; ___2134 < ___467; ++___2134) { int64_t const index = (___2134 * nodeJMax + ___2105) * nodeIMax + i; nodes.resize(0); nodes.push_back(index); if (nodeIMax > 1) { nodes.push_back(index + 1); } if (nodeJMax > 1) { size_t count = nodes.size(); for(size_t n = 0; n < count; ++n) { nodes.push_back(nodes[n] + nodeIMax); } } if (nodeKMax > 1) { size_t count = nodes.size(); for(size_t n = 0; n < count; ++n) { nodes.push_back(nodes[n] + nodeIMax * nodeJMax);
} } double ___4314 = 0.0; for(size_t n = 0; n < nodes.size(); ++n) ___4314 += nativeFieldData->___1780(nodes[n]); derivedFieldData->___3504(index, ___4314 / nodes.size()); } } } } else { ___478(___2682.___1669() > 0);
 #if !defined TECIOMPI
int numThreads = 1; if (___2497->___2392 >= MIN_CELLS_FOR_MULTITHREAD) numThreads = std::min(___2122::___2827(), static_cast<int>((___2497->___2392 - 1) / MIN_CELLS_FOR_MULTITHREAD + 1)); if (numThreads == 1) {
 #endif
FECCValueDerivationData derivationData(0, ___2497->___2392, nativeFieldData.get(), derivedFieldData.get(), ___2497.get()); deriveRangeOfCCValues((___90)&derivationData); derivedFieldData->___3499(derivationData.minVal, derivationData.maxVal);
 #if !defined TECIOMPI
} else { std::vector<boost::shared_ptr<FECCValueDerivationData> > derivationData; for (int i = 0; i < numThreads; ++i) { ___465 const begin = static_cast<___465>((size_t)___2497->___2392 * i / numThreads); ___465 const end = static_cast<___465>((size_t)___2497->___2392 * (i + 1) / numThreads); derivationData.push_back(boost::make_shared<FECCValueDerivationData>( begin, end, nativeFieldData.get(), derivedFieldData.get(), ___2497.get())); } ___2122 ___2119; for (int i = 0; i < numThreads; ++i) ___2119.addJob(deriveRangeOfCCValues, reinterpret_cast<___90>(derivationData[i].get())); ___2119.wait(); double minVal = std::numeric_limits<double>::max(); double maxVal = -std::numeric_limits<double>::max(); for (int i = 0; i < numThreads; ++i) { minVal = std::min(minVal, derivationData[i]->minVal); maxVal = std::max(maxVal, derivationData[i]->maxVal); } derivedFieldData->___3499(minVal, maxVal); }
 #endif
} } namespace { int64_t positionInArray(std::vector<int64_t> const& array, int64_t ___2085) { REQUIRE("array is sorted."); REQUIRE(___2085 >= 0); std::vector<int64_t>::const_iterator it = std::lower_bound(array.begin(), array.end(), ___2085); if (it != array.end() && *it == ___2085) return static_cast<int64_t>(std::distance(array.begin(), it)); else return -1; } void filterGhostInfo(GhostInfo_s& ghostInfo, std::set<int64_t> const& selectedItems, std::vector<int64_t> const& oldItemNumbers) { REQUIRE(selectedItems.size() == oldItemNumbers.size()); if (!ghostInfo.m_items.empty()) { GhostInfo_s newGhostInfo; size_t numSelectedGhostItems = 0; for (size_t i = 0; i < ghostInfo.m_items.size(); ++i) if (selectedItems.find(ghostInfo.m_items[i] - 1) != selectedItems.end()) ++numSelectedGhostItems; newGhostInfo.m_items.reserve(numSelectedGhostItems); newGhostInfo.m_neighbors.reserve(numSelectedGhostItems); newGhostInfo.m_neighborItems.reserve(numSelectedGhostItems); for (size_t i = 0; i < ghostInfo.m_items.size(); ++i) { if (selectedItems.find(ghostInfo.m_items[i] - 1) != selectedItems.end()) { newGhostInfo.m_items.push_back(positionInArray(oldItemNumbers, ghostInfo.m_items[i] - 1) + 1); newGhostInfo.m_neighbors.push_back(ghostInfo.m_neighbors[i]); newGhostInfo.m_neighborItems.push_back(ghostInfo.m_neighborItems[i]); } } newGhostInfo.m_items.swap(ghostInfo.m_items); newGhostInfo.m_neighbors.swap(ghostInfo.m_neighbors); newGhostInfo.m_neighborItems.swap(ghostInfo.m_neighborItems); } }
 #if defined TECIOMPI
void renumberGhostNeighborNodes( MPI_Comm comm, std::vector<int> const& partitionOwners, Zone_s::ZoneMap& partitionMap) { if (teciompi::everyRankAppearsOnce(comm, partitionOwners)) { Zone_s::Ptr partitionPtr = partitionMap.begin()->second; int commSize; MPI_Comm_size(comm, &commSize); std::vector<int32_t> ghostNodeCountToRank(commSize, 0); for (size_t i = 0; i < partitionPtr->m_ghostNodeInfo.m_neighbors.size(); ++i) { int32_t neighbor = partitionPtr->m_ghostNodeInfo.m_neighbors[i]; int owningRank = partitionOwners[neighbor - 1]; ++ghostNodeCountToRank[owningRank]; } std::vector<int32_t> ghostNodeCountFromRank(commSize); MPI_Alltoall(&ghostNodeCountToRank[0], 1, MPI_INT32_T, &ghostNodeCountFromRank[0], 1, MPI_INT32_T, comm); std::vector<int> sendDisplacements(commSize, 0); for (int i = 1; i < commSize; ++i) sendDisplacements[i] = sendDisplacements[i - 1] + ghostNodeCountToRank[i - 1]; int totalSendSize = std::max(1, sendDisplacements[commSize - 1] + ghostNodeCountToRank[commSize - 1]); std::vector<int64_t> ghostNodesToSend(totalSendSize); std::vector<int> index(sendDisplacements); for (size_t i = 0; i < partitionPtr->m_ghostNodeInfo.m_neighbors.size(); ++i) { int32_t neighbor = partitionPtr->m_ghostNodeInfo.m_neighbors[i]; int owningRank = partitionOwners[neighbor - 1]; ghostNodesToSend[index[owningRank]] = partitionPtr->m_ghostNodeInfo.m_neighborItems[i]; ++index[owningRank]; } std::vector<int> receiveDisplacements(commSize, 0); for (int i = 1; i < commSize; ++i) receiveDisplacements[i] = receiveDisplacements[i - 1] + ghostNodeCountFromRank[i - 1]; int totalReceiveSize = std::max(1, receiveDisplacements[commSize - 1] + ghostNodeCountFromRank[commSize - 1]); std::vector<int64_t> ghostNodesToReceive(totalReceiveSize); MPI_Alltoallv(&ghostNodesToSend[0], &ghostNodeCountToRank[0], &sendDisplacements[0], MPI_INT64_T, &ghostNodesToReceive[0], &ghostNodeCountFromRank[0], &receiveDisplacements[0], MPI_INT64_T, comm); for(int i = 0; i < totalReceiveSize; ++i) ghostNodesToReceive[i] = positionInArray(partitionPtr->m_oldNodeNumbers, ghostNodesToReceive[i] - 1) + 1; MPI_Alltoallv(&ghostNodesToReceive[0], &ghostNodeCountFromRank[0], &receiveDisplacements[0], MPI_INT64_T, &ghostNodesToSend[0], &ghostNodeCountToRank[0], &sendDisplacements[0], MPI_INT64_T, comm); index = sendDisplacements; for (size_t i = 0; i < partitionPtr->m_ghostNodeInfo.m_neighbors.size(); ++i) { int32_t neighbor = partitionPtr->m_ghostNodeInfo.m_neighbors[i]; int owningRank = partitionOwners[neighbor - 1]; partitionPtr->m_ghostNodeInfo.m_neighborItems[i] = ghostNodesToSend[index[owningRank]]; ++index[owningRank]; } } else { throw std::runtime_error("Not implemented."); } }
 #else
void renumberGhostNeighborNodes(Zone_s::ZoneMap& partitionMap) { BOOST_FOREACH(Zone_s::ZoneMap::value_type& valuePair, partitionMap) { Zone_s::Ptr& partitionPtr = valuePair.second; for (size_t i = 0; i < partitionPtr->m_ghostNodeInfo.m_items.size(); ++i) { int32_t neighbor   = partitionPtr->m_ghostNodeInfo.m_neighbors[i]; int64_t oldNodeNum = partitionPtr->m_ghostNodeInfo.m_neighborItems[i]; int64_t newNodeNum = positionInArray(partitionMap[neighbor - 1]->m_oldNodeNumbers, oldNodeNum - 1) + 1; partitionPtr->m_ghostNodeInfo.m_neighborItems[i] = newNodeNum; } } }
 #endif
} void Zone_s::filterData(
 #if defined TECIOMPI
MPI_Comm comm, int mainProcess,
 #endif
int32_t ___4336, ___2479 const& minMax) { REQUIRE(___4336 > 0); REQUIRE(minMax.___2067()); if (___2684 != ___4695 && ___2684 != ___4701) { throw std::runtime_error("Not implemented."); } if (m_partitionMap.empty()) { boost::shared_ptr<___1362> ___1351 = ___2496[___4336 - 1]; ___478(___1351->___2669 == ___4330); std::set<int64_t> selectedCells; std::set<int64_t> ghostCellSet(m_ghostCellInfo.m_items.begin(), m_ghostCellInfo.m_items.end()); for (int64_t ___449 = 0; ___449 < ___2497->___2392; ++___449) { ___2479 cellMinMax; for (int32_t ___681 = 0; ___681 < ___2497->___2500; ++___681) { int64_t index = ___681 + ___449 * ___2497->___2500; int64_t ___2709 = ___2497->___4314(index); double ___4314 = ___1351->___1780(___2709); cellMinMax.include(___4314); } if (minMax.intersects(cellMinMax) && ghostCellSet.find(___449 + 1) == ghostCellSet.end()) selectedCells.insert(selectedCells.end(), ___449); } if (selectedCells.empty()) { for (int64_t ___449 = 0; ___449 < ___2497->___2392; ++___449) if (ghostCellSet.find(___449 + 1) == ghostCellSet.end()) { selectedCells.insert(___449); break; } } ___478(!selectedCells.empty()); std::vector<int64_t> const oldCellNumbers(selectedCells.begin(), selectedCells.end()); std::set<int64_t> selectedNodes; BOOST_FOREACH(int64_t const ___449, selectedCells) { for (int32_t ___681 = 0; ___681 < ___2497->___2500; ++___681) { int64_t const index = ___681 + ___449 * ___2497->___2500; int64_t const ___2709 = ___2497->___4314(index); selectedNodes.insert(___2709); } } m_oldNodeNumbers.assign(selectedNodes.begin(), selectedNodes.end()); int64_t const newNumCells = static_cast<int64_t>(oldCellNumbers.size()); int64_t const newNumNodes = static_cast<int64_t>(m_oldNodeNumbers.size()); boost::shared_ptr<___2730> filteredNodeMap = ___2730::makePtr(___2497->___2500, newNumCells, m_oldNodeNumbers[newNumNodes - 1]); for (int64_t newCell = 0; newCell < newNumCells; ++newCell) { for (int32_t ___681 = 0; ___681 < ___2497->___2500; ++___681) { int64_t const index = ___681 + oldCellNumbers[newCell] * ___2497->___2500; int64_t const oldNode = ___2497->___4314(index); int64_t const newNode = positionInArray(m_oldNodeNumbers, oldNode); filteredNodeMap->appendValue(newNode); } } ___2497 = filteredNodeMap; TypedFieldDataFactory factory; for (size_t i = 0; i < ___2496.size(); ++i) { ___1362::Ptr oldFieldData = ___2496[i]; if (oldFieldData->storedValueCount() == 0) continue; ___1362::Ptr newFieldData = factory.make(oldFieldData->___2459); if (oldFieldData->___2669 == ___4328) { throw std::runtime_error("Not implemented.");
 #if 0
newFieldData->reserveValues(newNumCells); for (int64_t n = 0; n < newNumCells; ++n) { int64_t const ___449 = oldCellNumbers[n]; double const ___4314 = oldFieldData->___1780(___449); newFieldData->appendValue(___4314); }
 #endif
} else { newFieldData->reserveValues(newNumNodes); for (int64_t n = 0; n < newNumNodes; ++n) { int64_t const ___2709 = m_oldNodeNumbers[n]; double const ___4314 = oldFieldData->___1780(___2709); newFieldData->appendValue(___4314); } } ___2496[i] = newFieldData; } filterGhostInfo(m_ghostNodeInfo, selectedNodes, m_oldNodeNumbers); filterGhostInfo(m_ghostCellInfo, selectedCells, oldCellNumbers); ___2682.setFENumCells(newNumCells); ___2682.setFENumNodes(newNumNodes); } else { BOOST_FOREACH(ZoneMap::value_type & valuePair, m_partitionMap) { Zone_s::Ptr partitionPtr = valuePair.second; partitionPtr->filterData(
 #if defined TECIOMPI
comm, mainProcess,
 #endif
___4336, minMax); } renumberGhostNeighborNodes(
 #if defined TECIOMPI
comm, m_partitionOwners,
 #endif
m_partitionMap); ___1844 newZoneSize(0, 0, ___2682.___1669()); BOOST_FOREACH(Zone_s::ZoneMap::value_type & valuePair, m_partitionMap) { Zone_s::Ptr& partitionPtr = valuePair.second; size_t writeIndex = 0; for (size_t readIndex = 0; readIndex < partitionPtr->m_ghostNodeInfo.m_items.size(); ++readIndex) { if (partitionPtr->m_ghostNodeInfo.m_neighborItems[readIndex] > 0) { partitionPtr->m_ghostNodeInfo.m_items[writeIndex]         = partitionPtr->m_ghostNodeInfo.m_items[readIndex]; partitionPtr->m_ghostNodeInfo.m_neighbors[writeIndex]     = partitionPtr->m_ghostNodeInfo.m_neighbors[readIndex]; partitionPtr->m_ghostNodeInfo.m_neighborItems[writeIndex] = partitionPtr->m_ghostNodeInfo.m_neighborItems[readIndex]; ++writeIndex; } } partitionPtr->m_ghostNodeInfo.m_items.resize(writeIndex); partitionPtr->m_ghostNodeInfo.m_neighbors.resize(writeIndex); partitionPtr->m_ghostNodeInfo.m_neighborItems.resize(writeIndex); newZoneSize.setFENumNodes(newZoneSize.___1670() + partitionPtr->___2682.___1670() - partitionPtr->m_ghostNodeInfo.m_items.size()); newZoneSize.setFENumCells(newZoneSize.___1668() + partitionPtr->___2682.___1668() - partitionPtr->m_ghostCellInfo.m_items.size()); }
 #if defined TECIOMPI
if (teciompi::everyRankAppearsOnce(comm, m_partitionOwners)) { int64_t localZoneSize[2] = { newZoneSize.___1670(), newZoneSize.___1668() }; int64_t sumZoneSize[2]; MPI_Reduce(localZoneSize, sumZoneSize, 2, MPI_INT64_T, MPI_SUM, mainProcess, comm); newZoneSize.setFENumNodes(sumZoneSize[0]); newZoneSize.setFENumCells(sumZoneSize[1]); } else { throw std::runtime_error("Not implemented."); }
 #endif
___2682 = newZoneSize; } } Zone_s::Zone_s() {} void Zone_s::writeToFile(FileWriterInterface& outputFile, bool ___4480) const { ___4544(outputFile, ___2683, ___4480); writeScalar(outputFile, (uint32_t)___2684, ___4480); writeScalar(outputFile, m_partitionOffset.i(), ___4480); writeScalar(outputFile, m_partitionOffset.___2105(), ___4480); writeScalar(outputFile, m_partitionOffset.___2134(), ___4480); writeScalar(outputFile, ___2682.i(), ___4480); writeScalar(outputFile, ___2682.___2105(), ___4480); writeScalar(outputFile, ___2682.___2134(), ___4480); writeScalar(outputFile, ___2621, ___4480); writeScalar(outputFile, ___2622, ___4480); writeScalar(outputFile, ___2614, ___4480); writeScalar(outputFile, ___2503, ___4480); writeScalar(outputFile, (uint32_t)___2458, ___4480); writeScalar(outputFile, ___2651, ___4480); writeScalar(outputFile, ___2501, ___4480); writeScalar(outputFile, ___2650, ___4480); std::vector<uint32_t>tempFieldDataTypes(___2460.begin(), ___2460.end()); writeVector(outputFile, tempFieldDataTypes, ___4480); writeVector(outputFile, m_passiveVars, ___4480); std::vector<uint32_t> tempValueLocations(___2670.begin(), ___2670.end()); writeVector(outputFile, tempValueLocations, ___4480); writeVector(outputFile, m_shareVarFromZone, ___4480); writeScalar(outputFile, m_shareConnectivityFromZone, ___4480); m_ghostNodeInfo.writeToFile(outputFile, ___4480); m_ghostCellInfo.writeToFile(outputFile, ___4480); writeVectorOfPtrs(outputFile, ___2496, ___4480); writeScalar(outputFile, (uint64_t)___2397, ___4480); if (___2684 != ___4704) ___2497->writeToFile(outputFile, ___4480); writeMapOfPairsToObjects(outputFile, ___2457, ___4480); ___2345->writeToFile(outputFile, ___4480); writeMapOfScalarsToPtrs(outputFile, m_partitionMap, ___4480); writeVector(outputFile, m_partitionOwners, ___4480); writeVector(outputFile, m_oldNodeNumbers, ___4480); } uint64_t Zone_s::sizeInFile(bool ___4480) const { uint64_t sizeInFile = 0; sizeInFile += stringSizeInFile(___2683, ___4480); sizeInFile += scalarSizeInFile((uint32_t)___2684, ___4480); sizeInFile += scalarSizeInFile(m_partitionOffset.i(), ___4480); sizeInFile += scalarSizeInFile(m_partitionOffset.___2105(), ___4480); sizeInFile += scalarSizeInFile(m_partitionOffset.___2134(), ___4480); sizeInFile += scalarSizeInFile(___2682.i(), ___4480); sizeInFile += scalarSizeInFile(___2682.___2105(), ___4480); sizeInFile += scalarSizeInFile(___2682.___2134(), ___4480); sizeInFile += scalarSizeInFile(___2621, ___4480); sizeInFile += scalarSizeInFile(___2622, ___4480); sizeInFile += scalarSizeInFile(___2614, ___4480); sizeInFile += scalarSizeInFile(___2503, ___4480); sizeInFile += scalarSizeInFile((uint32_t)___2458, ___4480);
sizeInFile += scalarSizeInFile(___2651, ___4480); sizeInFile += scalarSizeInFile(___2501, ___4480); sizeInFile += scalarSizeInFile(___2650, ___4480); std::vector<uint32_t> tempFieldDataTypes(___2460.size()); sizeInFile += vectorSizeInFile(tempFieldDataTypes, ___4480); sizeInFile += vectorSizeInFile(m_passiveVars, ___4480); std::vector<uint32_t> tempValueLocations(___2670.size()); sizeInFile += vectorSizeInFile(tempValueLocations, ___4480); sizeInFile += vectorSizeInFile(m_shareVarFromZone, ___4480); sizeInFile += scalarSizeInFile(m_shareConnectivityFromZone, ___4480); sizeInFile += m_ghostNodeInfo.sizeInFile(___4480); sizeInFile += m_ghostCellInfo.sizeInFile(___4480); sizeInFile += vectorOfPtrsSizeInFile(___2496, ___4480); sizeInFile += scalarSizeInFile((uint64_t)___2397, ___4480); if (___2684 != ___4704) sizeInFile += ___2497->sizeInFile(___4480); sizeInFile += mapOfPairsToObjectsSizeInFile(___2457, ___4480); sizeInFile += ___2345->sizeInFile(___4480); sizeInFile += mapOfScalarsToPtrsSizeInFile(m_partitionMap, ___4480); sizeInFile += vectorSizeInFile(m_partitionOwners, ___4480); sizeInFile += vectorSizeInFile(m_oldNodeNumbers, ___4480); return sizeInFile; } boost::shared_ptr<Zone_s> Zone_s::makePtr(___1399& inputFile, bool readASCII) { Zone_s* newZone = new Zone_s; readString(inputFile, newZone->___2683, readASCII); READ_ENUM(newZone->___2684, ZoneType_e, inputFile, readASCII); ___81 i, ___2105, ___2134; readScalar(inputFile, i, readASCII); readScalar(inputFile, ___2105, readASCII); readScalar(inputFile, ___2134, readASCII); newZone->m_partitionOffset = ___1844(i, ___2105, ___2134); readScalar(inputFile, i, readASCII); readScalar(inputFile, ___2105, readASCII); readScalar(inputFile, ___2134, readASCII); newZone->___2682 = ___1844(i, ___2105, ___2134); readScalar(inputFile, newZone->___2621, readASCII); readScalar(inputFile, newZone->___2622, readASCII); readScalar(inputFile, newZone->___2614, readASCII); readScalar(inputFile, newZone->___2503, readASCII); READ_ENUM(newZone->___2458, FaceNeighborMode_e, inputFile, readASCII); readScalar(inputFile, newZone->___2651, readASCII); readScalar(inputFile, newZone->___2501, readASCII); readScalar(inputFile, newZone->___2650, readASCII); READ_ENUM_VECTOR(newZone->___2460, FieldDataType_e, inputFile, readASCII); readVector(inputFile, newZone->m_passiveVars, readASCII); READ_ENUM_VECTOR(newZone->___2670, ValueLocation_e, inputFile, readASCII); readVector(inputFile, newZone->m_shareVarFromZone, readASCII); readScalar(inputFile, newZone->m_shareConnectivityFromZone, readASCII); newZone->m_ghostNodeInfo = GhostInfo_s(inputFile, readASCII); newZone->m_ghostCellInfo = GhostInfo_s(inputFile, readASCII); readVectorOfPtrs(inputFile, newZone->___2496, readASCII); newZone->___2400.resize(newZone->___2496.size());
readScalar(inputFile, newZone->___2397, readASCII); if (newZone->___2684 != ___4704) newZone->___2497 = ___2730::makePtr(inputFile, readASCII); readMapOfPairsToObjects(inputFile, newZone->___2457, readASCII); newZone->___2345 = AuxData_s::makePtr(inputFile, readASCII); readMapOfScalarsToPtrs(inputFile, newZone->m_partitionMap, readASCII); readVector(inputFile, newZone->m_partitionOwners, readASCII); readVector(inputFile, newZone->m_oldNodeNumbers, readASCII); int64_t nodeCount; int64_t cellCount; getZoneCounts(nodeCount, cellCount, newZone->___2682, newZone->___2684); for(size_t ___4291 = 0; ___4291 < newZone->___2400.size(); ++___4291) { newZone->___2400[___4291] = TypedFieldDataFactory().make(newZone->___2460[___4291]); if (newZone->___2670[___4291] == ___4330) { newZone->___2400[___4291]->___2669 = ___4328; newZone->___2400[___4291]->___2668 = cellCount; } else { newZone->___2400[___4291]->___2669 = ___4330; newZone->___2400[___4291]->___2668 = nodeCount; } } return boost::shared_ptr<Zone_s>(newZone); } }}
