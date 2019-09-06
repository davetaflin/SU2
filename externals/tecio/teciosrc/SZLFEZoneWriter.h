 #pragma once
#include "ZoneWriterAbstract.h"
#include "FieldData.h"
#include "NodeMap.h"
#include "SZLFEZoneHeaderWriter.h"
namespace tecplot { namespace ___3933 { class ___1339; class ___1350; class ItemSetIterator; class SZLFEZoneWriter : public ___4709 { public: SZLFEZoneWriter( ItemSetIterator&                           varIter, ___4636                                zone, ___4636                                ___341, std::vector<___372> const&              ___4564, ___372                                  ___4499, ___37&                                ___36, boost::shared_ptr<___1350 const> const& zoneInfo); virtual ~SZLFEZoneWriter(); static uint64_t fieldDataSubzoneHeaderFileSize(bool ___2002); static uint64_t cszConnectivityHeaderFileSize(bool ___2002); protected: SZLFEZoneHeaderWriter m_headerWriter; boost::shared_ptr<___1350 const> m_feZoneInfo; void setZoneNumberLabel(std::string const& zoneNumberLabel); private: std::string m_zoneNumberLabel; ___1392 ___2673; ___1392 m_cszConnectivityFileLocs; ___1392 m_nszConnectivityFileLocs; UInt16Array m_numRefNodeSubzones; UInt16Array m_numRefCellSubzones; UInt8Array m_cszIncludesPartitionOffsetsBitArray; UInt8Array m_nszIncludesPartitionOffsetsBitArray; virtual uint64_t zoneConnectivityFileSize(bool ___2002); virtual uint64_t zoneDataFileSize(bool ___2002); virtual uint64_t zoneHeaderFileSize(bool ___2002); virtual ___372 writeZoneConnectivity(FileWriterInterface& szpltFile); virtual ___372 writeZoneData(FileWriterInterface& szpltFile); virtual ___372 writeZoneHeader(FileWriterInterface& szpltFile); ___372 calculateVarSubzoneMinMaxes( ___4352   ___4336, ___2481& nszDataMinMax, ___2481& cszDataMinMax); ___372 ___4512( FileWriterInterface&         file, ValueLocation_e              ___4326, ___4352                   ___4336, ___2090::SubzoneOffset_t ___3880); template <typename T, bool isBitArray> ___372 ___4496( FileWriterInterface&         szpltFile, ___1352 const&             ___1351, ___2090::SubzoneOffset_t ___469); template <typename T, bool isBitArray> ___372 ___4531( FileWriterInterface&         szpltFile, ___1352 const&             ___1351, ___2090::SubzoneOffset_t ___2734); template <typename T, bool isBitArray> uint64_t cellSubzoneFieldDataFileSize(bool ___2002, ___2090::SubzoneOffset_t ___469) const; template <typename T, bool isBitArray> uint64_t nodeSubzoneFieldDataFileSize(bool ___2002, ___2090::SubzoneOffset_t ___2734) const; template <typename T, bool isBitArray> uint64_t subzoneFieldDataFileSize(bool ___2002, ___2090::SubzoneOffset_t ___3880, ValueLocation_e ___4326) const; template <typename T, bool isBitArray> ___372 writeVariable( FileWriterInterface&     szpltFile, ___4352 const         ___4336, ___2481 const&       nszDataMinMax, ___2481 const&       cszDataMinMax); template <typename T, bool isBitArray> uint64_t variableFileSize(bool ___2002, ValueLocation_e ___4326);
___372 ___4501( FileWriterInterface&         file, ___2090::SubzoneOffset_t ___469); ___372 ___4500( FileWriterInterface&           file, ___1339 const& compressor, bool                           outputPartitionIndices); uint64_t cszConnectivityDataFileSize( bool ___2002, size_t totalNumCellCorners, size_t numRefNszs, ___2090::___2980 numRefPartitions, bool outputPartitionIndices); ___372 writeCszConnectivity( FileWriterInterface&      szpltFile, ___2729                ___2723, PartitionSubzoneSetArray& nszRefCszSets); uint64_t cszConnectivityFileSize(bool ___2002, PartitionSubzoneSetArray& nszRefPtnCszSets); ___372 ___4537( FileWriterInterface&         file, ___2090::SubzoneOffset_t ___2734); uint64_t nszConnectivityHeaderFileSize(bool ___2002); ___372 ___4535( FileWriterInterface&      file, PartitionSubzoneSetArray& nszRefPtnCszSets); uint64_t nszConnectivityFileSize(bool ___2002, PartitionSubzoneSetArray& nszRefPtnCszSets); PartitionSubzoneSet getRefPtnCszs(___2729 ___2723, ___2090::SubzoneOffset_t ___2734); }; }}
