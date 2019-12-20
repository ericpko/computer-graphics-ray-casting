#include "TriangleSoup.h"
#include "Ray.h"
// Hint
#include "first_hit.h"

bool TriangleSoup::intersect(
  const Ray & ray, const double min_t, double & t, Eigen::Vector3d & n) const
{
  ////////////////////////////////////////////////////////////////////////////
  /**
   * NOTE: We don't need to iterate through triangles since first_hit does
   * this for us.
   */
  int hit_id;
  return first_hit(ray, min_t, TriangleSoup::triangles, hit_id, t, n);
  ////////////////////////////////////////////////////////////////////////////
}
