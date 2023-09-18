// Nicolo Grilli
// ADSC Singapore
// 3 October 2020

#pragma once

#include "GeneralUserObject.h"

/**
 * Read properties from file - grain or element
 * Input file syntax: prop1 prop2 etc. See test.
 * For grain level, voronoi tesellation with random grain centers are generated;
 * Element center points used for assigning properties
 * Usable for generated mesh
 * using read_type = indexgrain
 * it is possible to read grains as physical volumes assigned in GMSH
 */

class GrainPropertyReadFile : public GeneralUserObject
{
public:
  static InputParameters validParams();

  GrainPropertyReadFile(const InputParameters & parameters);
  virtual ~GrainPropertyReadFile() {}

  virtual void initialize() {}
  virtual void execute() {}
  virtual void finalize() {}

  /**
   * This function reads element data from file
   */
  void readElementData();

  /**
   * This function Read grain data from file
   */
  virtual void readGrainData();

  /**
   * This function generates grain center point
   * Presently random generated
   */
  virtual void initGrainCenterPoints();

  /**
   * This function assign property data to elements
   */
  Real getData(const Elem *, unsigned int) const;

  /**
   * This function assign properties to element read from file with element based properties
   */
  Real getElementData(const Elem *, unsigned int) const;

  /**
   * This function assign properties to element read from file with grain  based properties
   * Grain distribution in the RVE can be Periodic or non-periodic (default)
   */
  Real getGrainData(const Elem *, unsigned int) const;

  /**
   * This function assign properties to element read from file with grain based properties
   * Grains are assigned based on the physical volumes in GMSH generated mesh
   */
  Real getIndexGrainData(const Elem *, unsigned int) const;

  /**
   * This function calculates minimum distance between 2 points
   * considering periodicity of the simulation volume
   */
  Real minPeriodicDistance(Point, Point) const;

protected:
  ///Name of file containing property values
  std::string _prop_file_name;
  ///Store property values read from file
  std::vector<Real> _data;
  ///Number of properties in a row
  unsigned int _nprop;
  ///Number of grains (for property read based on grains)
  unsigned int _ngrain;
  ///Type of read - element or grain or indexgrain
  MooseEnum _read_type;
  ///Random seed - used for generating grain centers
  unsigned int _rand_seed;
  ///Type of grain structure - non-periodic default
  MooseEnum _rve_type;

  MooseMesh & _mesh;
  std::vector<Point> _center;

private:
  unsigned int _nelem;
  Point _top_right;
  Point _bottom_left;
  Point _range;
  Real _max_range;
};
