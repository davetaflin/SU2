#include "stdafx.h"
#include "MASTER.h"
 #define DATASETMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "CHARTYPE.h"
#include "STRUTIL.h"
#include "AUXDATA.h"
#include "ARRLIST.h"
#include "STRLIST.h"
#include "ALLOC.h"
#include "SET.h"
#include "DATASET.h"
#include "FILESTREAM.h"
#include "DATASET0.h"
#include "stringformat.h"
#include <float.h>
void ___490(___4683 *___4677) { REQUIRE(VALID_REF(___4677)); if (___4677->___2686 != NULL) ___1530(___4677->___2686, "ZoneSpec name"); if (___4677->___230 != NULL) ___236(&___4677->___230); ___3541(___4677); } void ___4679(___4683 **___4677) { REQUIRE(VALID_REF(___4677)); REQUIRE(VALID_REF(*___4677) || *___4677 == NULL); if (*___4677 != NULL) { ___490(*___4677); ___1531(*___4677, "ZoneSpec structure"); *___4677 = NULL; } ENSURE(*___4677 == NULL); } ___372 ___4681(void       *___2098, ___90 ___494) { ___4683 **ZoneSpecRef = (___4683 **)___2098; REQUIRE(VALID_REF(ZoneSpecRef)); REQUIRE(VALID_REF(*ZoneSpecRef) || *ZoneSpecRef == NULL); ___4278(___494); if (*ZoneSpecRef != NULL) ___4679(ZoneSpecRef); ENSURE(*ZoneSpecRef == NULL); return ___4226; } void ___3541(___4683 *___4677) { REQUIRE(VALID_REF(___4677)); ___4677->___2686                         = NULL; ___4677->___4263                     = ___1991; ___4677->___2975                   = ___333; ___4677->___3786                     = ___3788; ___4677->___3641                 = 0.0; ___4677->___2830                      = 0; ___4677->___2831                      = 0; ___4677->___2832                      = 0; ___4677->___4236                         = ___4704; ___4677->___2064                    = ___1305; ___4677->___4647.___3174 = ___2708; ___4677->___4647.___2027 = ___4226; ___4677->___230                      = NULL; ___4677->___419             = ___4226; ___4677->___1440                 = ___1290; ___4677->___228 = ___4226; ___4677->___2805      = 0; ___4677->___2799 = 0; ___4677->___2800 = 0; } ___4683 *___4678(void) { ___4683 *___3359; ___3359 = (___4683 *)ALLOC_ITEM(___4683, "ZoneSpec structure"); if (___3359 != NULL) ___3541(___3359); ENSURE(___3359 == NULL || VALID_REF(___3359)); return ___3359; } ___2227 ___4668(___134 ASSERT_ONLY(___4667), ___2227    ___693, ___2227    ___3354, ___90   ___494) { ___2227 ___3359; REQUIRE(ArrayListIsValid(___4667)); REQUIRE((___3354 == 0 && ___693 == 0) || ___3354 > ___693); REQUIRE(___693 <= ___2382); ___4278(___494); if (___3354 <= ___2382) { if (___3354 != 0 && ___693 == 0) { ___3359 = ___3354; } else { const ___2227 DEFAULT_CAPACITY = 32; ___2227       BlockSize = MAX(DEFAULT_CAPACITY, ___693 / 2); if (___3354 == 0) ___3359 = DEFAULT_CAPACITY; else ___3359 = ((___3354 - 1) / BlockSize + 1) * BlockSize; if (___3359 > ___2382) ___3359 = ___2382; } } else ___3359 = 0; ENSURE(___3359 == 0 || ___3359 >= ___3354);
ENSURE(___3359 <= ___2382); return ___3359; }
