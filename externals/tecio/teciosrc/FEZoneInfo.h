 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <algorithm>
#include <utility>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "FESubzonePartitionerInterface.h"
#include "fileStuff.h"
#include "SzlFileLoader.h"
namespace tecplot { namespace ___3933 { class ___1350 { UNCOPYABLE_CLASS(___1350); public: ___1350( ___465                                      ___2781, ___2718                                      ___2821, ___682                                    ___2789, boost::shared_ptr<FESubzonePartitionerInterface> subzonePartitioner, ___2090::___2980                         ___2977 = ___2090::UNKNOWN_PARTITION, ___465                                      numGhostCells = 0, ___2718                                      numGhostNodes = 0) : m_numCells(___2781) , m_numNodes(___2821) , m_numCorners(___2789) , m_subzonePartitioner(subzonePartitioner) , m_partition(___2977) , m_numGhostCells(numGhostCells) , m_numGhostNodes(numGhostNodes) , m_numReferencedPartitions(0) {} ~___1350() {} typedef boost::unordered_map<___2090::SubzoneOffset_t, boost::unordered_set<PartitionSubzone> > NeighborCszRefMap; typedef boost::unordered_map<___2090::SubzoneOffset_t, std::vector<___2479> > NszMinMaxMap; static ___2090::SubzoneOffset_t calcNumSubzones(___81 ___2812) { REQUIRE(___2812>0); REQUIRE(___2090::MAX_ITEM_OFFSET+1==DEFAULT_SUBZONE_MAX_FE_SIZE); ___2090::SubzoneOffset_t const numSubzones = ___2090::SubzoneOffset_t( 1 + ( (___2812-1) >> ___2090::ItemOffsetBitSize ) ); ENSURE(numSubzones > 0 && numSubzones <= ___2812); return numSubzones; } ___2090::___2980 getPartition() const { return m_partition; } ___465   ___1766() const { return m_numCells; } ___2718   ___1768() const { return m_numNodes; } ___682 ___1767() const { return m_numCorners; } ___465   getNumGhostCells() const { return m_numGhostCells; } ___2718   getNumGhostNodes() const { return m_numGhostNodes; } inline ___2090::SubzoneOffset_t ___2783() const { ___2090::SubzoneOffset_t const ___2783 = m_subzonePartitioner->___2783(); ENSURE(___2783 == 0 || ___2783 == calcNumSubzones(m_numCells - m_numGhostCells)); return ___2783; } ___2090::ItemOffset_t ___2782(___2090::SubzoneOffset_t ___469) const { return m_subzonePartitioner->___2782(___469); } ___465 ___4608(___2090 szCoordinate) const { return m_subzonePartitioner->___4608(szCoordinate); } inline ___2090::SubzoneOffset_t ___2823() const { ___2090::SubzoneOffset_t const ___2823 = m_subzonePartitioner->___2823(); ENSURE(___2823 == 0 || ___2823 == calcNumSubzones(m_numNodes - m_numGhostNodes)); return ___2823; } ___2090::ItemOffset_t ___2822(___2090::SubzoneOffset_t ___2734) const { return m_subzonePartitioner->___2822(___2734);} ___2718 ___4657(___2090 szCoordinate) const { return m_subzonePartitioner->___4657(szCoordinate); } void addNeighborNodeCoordinate(___2718 ___4656, ___2090 szCoordinate) { m_subzonePartitioner->setNodeSubzoneCoordinate(___4656, szCoordinate); }
inline void addNeighborCszInfo( tecplot::___2090::___2980 neighbor, std::vector<___2718> const& nodes, std::vector<boost::unordered_set<___2090::SubzoneOffset_t> > const& cellSubzones, std::vector<std::vector<___2479> > const& varMinMaxes) { REQUIRE(getPartition() != ___2090::INVALID_PARTITION); REQUIRE(neighbor != getPartition()); REQUIRE(!nodes.empty()); REQUIRE(nodes.size() == cellSubzones.size()); REQUIRE(cellSubzones.size() == varMinMaxes.size()); if (m_numReferencedPartitions == 0) { m_numReferencedPartitions = 2; m_referencedPartitions.alloc(m_numReferencedPartitions); m_referencedPartitions[0] = getPartition(); m_referencedPartitions[1] = neighbor; } else { std::set<tecplot::___2090::___2980> partitions(m_referencedPartitions.begin(), m_referencedPartitions.end(m_numReferencedPartitions)); if (partitions.insert(neighbor).second) { ++m_numReferencedPartitions; m_referencedPartitions.___937(); m_referencedPartitions.alloc(m_numReferencedPartitions); int i = 0; BOOST_FOREACH(tecplot::___2090::___2980 ___2977, partitions) m_referencedPartitions[i++] = ___2977; } } for(size_t i = 0; i < nodes.size(); ++i) { ___2090 szCoordinate = ___3924(nodes[i]); ___2090::SubzoneOffset_t ___2757 = szCoordinate.subzoneOffset(); BOOST_FOREACH(___2090::SubzoneOffset_t ___469, cellSubzones[i]) m_neighborCszRefs[___2757].insert(PartitionSubzone(neighbor, ___469)); std::vector<___2479>& minMaxes = m_nszMinMaxes[___2757]; if (minMaxes.empty()) minMaxes.resize(varMinMaxes[i].size()); for(size_t ___2105 = 0; ___2105 < varMinMaxes[i].size(); ++___2105) minMaxes[___2105].include(varMinMaxes[i][___2105]); } } inline void addNeighborCszRef(___2090::SubzoneOffset_t ___2757, ___2090 cellCoordinate) { m_neighborCszRefs[___2757].insert(PartitionSubzone(cellCoordinate.___2977(), cellCoordinate.subzoneOffset())); } NeighborCszRefMap const& getNeighborCszRefs() const { return m_neighborCszRefs; } NszMinMaxMap const& getNszMinMaxes() const { return m_nszMinMaxes; } inline void resetNeighborInfo() { m_neighborCszRefs.clear(); m_nszMinMaxes.clear(); m_numReferencedPartitions = 0; m_referencedPartitions.___937(); } inline ___2090 ___3924(___2718 ___4656) const { REQUIRE(0 <= ___4656 && ___4656 < m_numNodes); return m_subzonePartitioner->___3924(___4656); } inline ___2090 szCoordinateAtZoneCell(___2718 zoneCell) const { REQUIRE(0 <= zoneCell && zoneCell < m_numCells); return m_subzonePartitioner->szCoordinateAtZoneCell(zoneCell); } inline void addReferencedPartitions(std::set<tecplot::___2090::___2980> const& partitionSet) { REQUIRE(getPartition() != ___2090::INVALID_PARTITION); REQUIRE(!partitionSet.empty()); std::set<tecplot::___2090::___2980> partitions(partitionSet); for(tecplot::___2090::___2980 i = 0; i < m_numReferencedPartitions; ++i) partitions.insert(m_referencedPartitions[i]); partitions.insert(getPartition()); if (m_numReferencedPartitions != static_cast<tecplot::___2090::___2980>(partitions.size()))
{ m_referencedPartitions.___937(); m_referencedPartitions.alloc((uint64_t)partitions.size()); int i = 0; BOOST_FOREACH(tecplot::___2090::___2980 ___2977, partitions) { m_referencedPartitions[i] = ___2977; ++i; } m_numReferencedPartitions = i; } } tecplot::___2090::___2980 getNumReferencedPartitions() const { return m_numReferencedPartitions; } PartitionArray const& getReferencedPartitions() const { return m_referencedPartitions; } private: ___465 const                                m_numCells; ___2718 const                                m_numNodes; ___682 const                              m_numCorners; boost::shared_ptr<FESubzonePartitionerInterface> m_subzonePartitioner; ___2090::___2980                         m_partition; ___465 const                                m_numGhostCells; ___2718 const                                m_numGhostNodes; NeighborCszRefMap                                m_neighborCszRefs; NszMinMaxMap                                     m_nszMinMaxes; tecplot::___2090::___2980                m_numReferencedPartitions; PartitionArray                                   m_referencedPartitions; }; }}
