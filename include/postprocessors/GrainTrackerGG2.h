// Added function: When misorientation is below a certain threshold, grains merge

#pragma once

#include "GrainTrackerGG.h"
#include "EulerAngleProvider.h"  // used to provide Euler angle by weipeng
#include "CalculateMisorientationAngle.h" // calculated misorientation angle based on Euler angles by weipeng

class GrainTrackerGG2 : public GrainTrackerGG
{
public:
  static InputParameters validParams();

  GrainTrackerGG2(const InputParameters & parameters);
  virtual ~GrainTrackerGG2();
protected:
  virtual void trackGrains();

  virtual void remapGrains();

  // Get the adjacent grain id, used to calculate the GB misorientation angle
  unsigned int getTopoRelaGrainID(const FeatureData & grain_i);  

  // re-merge grains due to misorientation angle from euler angles calculation
  void mergeGrainsBasedMisorientation();  

  // True if two grains are determined to perform a merge operation
  bool _remerge_grains;  

  const EulerAngleProvider & _euler;
  misoriAngle_isTwining _s_misoriTwin;

  std::vector<unsigned int> _inactive_grains_id;
};