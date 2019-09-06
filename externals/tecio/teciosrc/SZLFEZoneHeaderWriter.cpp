#include "SZLFEZoneHeaderWriter.h"
#include "ThirdPartyHeadersBegin.h"
#include <boost/assign.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "AltTecUtil.h"
#include "FEZoneInfo.h"
#include "fileStuff.h"
#include "ItemSetIterator.h"
#include "SzlFileLoader.h"
#include "writeValueArray.h"
namespace tecplot { namespace ___3933 { SZLFEZoneHeaderWriter::SZLFEZoneHeaderWriter( ItemSetIterator&    varIter, ___4636         zone, ___4636         ___341, ___37&         ___36, ___1350 const&   ___1349, bool                ___4499, ___1392 const& varFileLocs, ___1392 const& ___687, ___1392 const& ___2758, UInt16Array const&  ___2834, UInt16Array const&  numRefCellSubzones, UInt8Array const&   cszIncludesPartitionOffsetsBitArray, UInt8Array const&   nszIncludesPartitionOffsetsBitArray) : ZoneHeaderWriterAbstract(varIter, zone, ___341, ___36) , m_feZoneInfo(___1349) , m_writeConnectivity(___4499) , ___2673(varFileLocs) , m_cszConnectivityFileLocs(___687) , m_nszConnectivityFileLocs(___2758) , m_numRefNodeSubzones(___2834) , m_numRefCellSubzones(numRefCellSubzones) , m_cszIncludesPartitionOffsetsBitArray(cszIncludesPartitionOffsetsBitArray) , m_nszIncludesPartitionOffsetsBitArray(nszIncludesPartitionOffsetsBitArray) {} SZLFEZoneHeaderWriter::~SZLFEZoneHeaderWriter() {} uint64_t SZLFEZoneHeaderWriter::sizeInFile(bool ___2002) const { uint64_t ___3358 = zoneHeaderTagsSizeInFile(9, ___2002); if (m_writeConnectivity) { size_t const ___2783 = static_cast<size_t>(m_feZoneInfo.___2783()); ___3358 += arraySizeInFile<uint64_t, true /* ___2025 */>(___2783, ___2002); size_t const ___2823 = static_cast<size_t>(m_feZoneInfo.___2823()); ___3358 += arraySizeInFile<uint64_t, true /* ___2025 */>(___2823, ___2002); ___3358 += arraySizeInFile<uint16_t, false>(___2783, ___2002); ___3358 += arraySizeInFile<uint16_t, false>(___2823, ___2002); ___3358 += arraySizeInFile<uint8_t, true /* ___2025 */>(numBytesForNumBits(___2783), ___2002); ___3358 += arraySizeInFile<uint8_t, true /* ___2025 */>(numBytesForNumBits(___2823), ___2002); } size_t const numVarsToWrite = static_cast<size_t>(m_varIter.___2812()); ___3358 += arraySizeInFile<uint64_t, true /* ___2025 */>(numVarsToWrite, ___2002); size_t const numReferencedPartitions = static_cast<size_t>(m_feZoneInfo.getNumReferencedPartitions()); if (numReferencedPartitions > 0) ___3358 += arraySizeInFile<uint32_t, false>(numReferencedPartitions, ___2002); return ___3358; } ___372 SZLFEZoneHeaderWriter::write(FileWriterInterface& fileWriter) const { REQUIRE(fileWriter.___2041()); ___4352 const numVarsToWrite = m_varIter.___2812(); REQUIRE(___2673.size() == uint64_t(numVarsToWrite)); REQUIRE(EQUIVALENCE(m_writeConnectivity, m_cszConnectivityFileLocs.size() > 0)); REQUIRE(IMPLICATION(m_writeConnectivity, m_cszConnectivityFileLocs.size() == static_cast<uint64_t>(m_feZoneInfo.___2783()))); REQUIRE(EQUIVALENCE(m_writeConnectivity, m_nszConnectivityFileLocs.size() > 0)); REQUIRE(IMPLICATION(m_writeConnectivity, m_nszConnectivityFileLocs.size() == static_cast<uint64_t>(m_feZoneInfo.___2823())));
___2090::SubzoneOffset_t const ___2783 = m_feZoneInfo.___2783(); ___2090::SubzoneOffset_t const ___2823 = m_feZoneInfo.___2823(); ___372 ___2039 = ___4226; try { ___3945 ___3944 = boost::assign::map_list_of<uint16_t, uint64_t> (___4342, ___330) (CSZ_CONNECT_FILE_LOC_TAG, ___330) (NSZ_CONNECT_FILE_LOC_TAG, ___330) (NUM_REF_PARTITIONS_TAG, uint64_t(m_feZoneInfo.getNumReferencedPartitions())) (REF_PARTITIONS_TAG, ___330) (NUM_REF_NODE_SUBZONES_TAG, ___330) (NUM_REF_CELL_SUBZONES_TAG, ___330) (CELL_SUBZONE_INCLUDES_PARTITIONS_TAG, ___330) (NODE_SUBZONE_INCLUDES_PARTITIONS_TAG, ___330); ___1393 headerFileLoc = fileWriter.fileLoc(); ___2039 = ___4565(fileWriter, ___3944); ___4636 const fileZone = ___2677 - m_baseZone; if (___2039 && m_writeConnectivity) { ___3944[CSZ_CONNECT_FILE_LOC_TAG] = fileWriter.fileLoc(); ___2039 = ___4563<uint64_t, true, 0>(fileWriter, appendZoneSuffix(CSZ_CONNECT_FILE_LOC_DESCRIPTION).c_str(), fileZone, ___2783, &m_cszConnectivityFileLocs[0]); } if (___2039 && m_writeConnectivity) { ___3944[NSZ_CONNECT_FILE_LOC_TAG] = fileWriter.fileLoc(); ___2039 = ___4563<uint64_t, true, 0>(fileWriter, appendZoneSuffix(NSZ_CONNECT_FILE_LOC_DESCRIPTION).c_str(), fileZone, ___2823, &m_nszConnectivityFileLocs[0]); } if (___2039) { ___3944[___4342] = fileWriter.fileLoc(); ___2039 = ___4563<uint64_t, true, 0>(fileWriter, appendZoneSuffix(VAR_FILE_LOC_DESCRIPTION).c_str(), fileZone, numVarsToWrite, &___2673[0]); } if (___2039 && m_feZoneInfo.getNumReferencedPartitions() > 0) { ___3944[REF_PARTITIONS_TAG] = fileWriter.fileLoc(); ___2039 = ___4563<uint32_t, false, 0>(fileWriter, appendZoneSuffix(REF_PARTITIONS_DESCRIPTION).c_str(), fileZone, static_cast<size_t>(m_feZoneInfo.getNumReferencedPartitions()), &m_feZoneInfo.getReferencedPartitions()[0]); } if (___2039 && m_writeConnectivity) { ___3944[NUM_REF_NODE_SUBZONES_TAG] = fileWriter.fileLoc(); ___2039 = ___4563<uint16_t, false, 0>(fileWriter, appendZoneSuffix(NUM_REF_NODE_SUBZONES_DESCRIPTION).c_str(), fileZone, ___2783, &m_numRefNodeSubzones[0]); } if (___2039 && m_writeConnectivity) { ___3944[NUM_REF_CELL_SUBZONES_TAG] = fileWriter.fileLoc(); ___2039 = ___4563<uint16_t, false, 0>(fileWriter, appendZoneSuffix(NUM_REF_CELL_SUBZONES_DESCRIPTION).c_str(), fileZone, ___2823, &m_numRefCellSubzones[0]); } if (___2039 && m_writeConnectivity) { ___3944[CELL_SUBZONE_INCLUDES_PARTITIONS_TAG] = fileWriter.fileLoc(); ___2039 = ___4563<uint8_t, true, 0>(fileWriter, appendZoneSuffix(CELL_SUBZONE_INCLUDES_PARTITIONS_DESCRIPTION).c_str(), fileZone, numBytesForNumBits(___2783), &m_cszIncludesPartitionOffsetsBitArray[0]); } if (___2039 && m_writeConnectivity) { ___3944[NODE_SUBZONE_INCLUDES_PARTITIONS_TAG] = fileWriter.fileLoc(); ___2039 = ___4563<uint8_t, true, 0>(fileWriter, appendZoneSuffix(NODE_SUBZONE_INCLUDES_PARTITIONS_DESCRIPTION).c_str(),
fileZone, numBytesForNumBits(___2823), &m_nszIncludesPartitionOffsetsBitArray[0]); } ___1393 endFileLoc = fileWriter.fileLoc(); ___2039 = ___2039 && fileWriter.___3459(headerFileLoc) && ___4565(fileWriter, ___3944) && fileWriter.___3459(endFileLoc); } catch(std::bad_alloc const&) { ___2039 = ___1186("Out of memory while writing zone %d header.", ___2677 + 1); } catch(...) { ___2039 = ___1186("Unrecoverable error while writing zone %d header.", ___2677 + 1); } ENSURE(VALID_BOOLEAN(___2039)); return ___2039; } }}
