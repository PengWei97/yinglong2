// Added function: When misorientation is below a certain threshold, grains merge

#pragma once

#include "GrainTrackerGG.h"

class GrainTrackerGG2 : public GrainTrackerGG
{
public:
  static InputParameters validParams();

  GrainTrackerGG2(const InputParameters & parameters);
  virtual ~GrainTrackerGG2();
};