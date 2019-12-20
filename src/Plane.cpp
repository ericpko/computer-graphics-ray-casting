#include "Plane.h"
#include "Ray.h"

bool Plane::intersect(
  const Ray & ray, const double min_t, double & t, Eigen::Vector3d & n) const
{
  ////////////////////////////////////////////////////////////////////////////
  bool intersect = false;

  // Set up some vectors to represent the plane and line
  Eigen::Vector3d p = Plane::point;                   // A specific point
  Eigen::Vector3d e = ray.origin;
  Eigen::Vector3d v = ray.direction;

  /**
   * First check if the direction of the line is perpendicular to the normal of
   * the plane. If this is the case, then the line does not inersect with plane
   */
  if (v.dot(Plane::normal) == 0) {
    return false;
  }

  // If we make it here, then we know that the ray intersects the plane at some
  // point t
  t = (Plane::normal.dot(p - e)) / (Plane::normal.dot(v)); // See ipad for calc
  if (t >= min_t) {
    n = Plane::normal;
    intersect = true;
  }

  return intersect;
  ////////////////////////////////////////////////////////////////////////////
}


/**
 * NOTE: I just found this out, but since we're in the namespace of Plane,
 * we can use the keyword this as a pointer
 * eg) this->point, this->normal, etc.
 * https://www.tutorialspoint.com/cplusplus/cpp_this_pointer.htm
 */
