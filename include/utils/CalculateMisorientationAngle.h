#pragma once

#include "Moose.h"
#include "EulerAngles.h"
#include "EulerAngleProvider.h"

typedef Eigen::Quaternion<Real> quatReal;

// set twin type
enum class TwinType {NONE, TT1, ST1};

struct MisorientationAngleData{Real _misor = 0.0; bool _is_twin = false; TwinType _twin_type = TwinType::NONE;};

/**
 * This is a real version of the misorientation angle between grains
 */
class CalculateMisorientationAngle
{
public:
  // function 1: input Euler1 and Euler2, output s
  static MisorientationAngleData calculateMisorientaion(EulerAngles & Euler1, EulerAngles & Euler2, MisorientationAngleData & s, const std::string & CrystalType = "hcp");  

  // function 2: Obtaining the key orientation using quaternion, including twinning, CS, SS
  static std::vector<quatReal> getKeyQuat(const std::string & QuatType, const std::string & CrystalKype = "hcp");

  // function 3.1: computer the scalar dot product using quaternion
  static Real dotQuaternion(const quatReal & o1, const quatReal & o2, 
                     const std::vector<quatReal> & qcs, 
                     const std::vector<quatReal> & qss);

  // function 3.2: computes inv(o1) .* o2 usig quaternion
  static quatReal itimesQuaternion(const quatReal & q1, const quatReal & q2);

  // function 3.3: computes outer inner product between two quaternions
  static Real dotOuterQuaternion(const quatReal & rot1, const std::vector<quatReal> & rot2);

  // function 3.4: X*Y is the matrix product of X and Y. ~twice~
  static Real mtimes2Quaternion(const quatReal & q1, const std::vector<quatReal> & q2, const quatReal & qTwin);  
};