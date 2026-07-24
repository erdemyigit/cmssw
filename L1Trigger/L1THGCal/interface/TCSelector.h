#ifndef L1Trigger_L1THGCal_TCselector_H
#define L1Trigger_L1THGCal_TCselector_H

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"

typedef StringCutObjectSelector<HGCalTriggerDetId> TCSelector;

#endif

