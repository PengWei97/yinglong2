// Added function: When misorientation is below a certain threshold, grains merge

#pragma once

#include "GrainTrackerBase.h"

class GTMergeMisori : public GrainTrackerBase
{
public:
  static InputParameters validParams();

  GrainTrackerMerge(const InputParameters & parameters);
};