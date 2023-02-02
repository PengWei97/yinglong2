#include "GrainTrackerGG2.h"

registerMooseObject("yinglongApp", GrainTrackerGG2);

InputParameters
GrainTrackerGG2::validParams()
{
  InputParameters params = GrainTrackerGG::validParams();
  params.addClassDescription("Grain Tracker object for running reduced order parameter simulations "
                             "without grain coalescence.");

  return params;
}

GrainTrackerGG2::GrainTrackerGG2(const InputParameters & parameters)
  : GrainTrackerGG(parameters)
{
}

GrainTrackerGG2::~GrainTrackerGG2() {}