#include "GTMergeMisori.h"

registerMooseObject("yinglongApp", GTMergeMisori);

InputParameters
GTMergeMisori::validParams()
{
  InputParameters params = GTMergeMisori::validParams();
  params.addClassDescription("Grain Tracker derived object for merging of grains based on misorientation angle.");

  return params;
}

GTMergeMisori::GTMergeMisori(const InputParameters & parameters)
  : GrainTrackerBase(parameters)
{
}