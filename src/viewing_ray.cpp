#include "viewing_ray.h"

void viewing_ray(
  const Camera & camera,
  const int i,
  const int j,
  const int width,
  const int height,
  Ray & ray)
{
  ////////////////////////////////////////////////////////////////////////////
  /**
   * We need a transformation T from (i, j) to the (u, v) plane
   * i.e. T(i, j) = (u, v)
   * hence, u = u(i, j) and v = v(i, j)
   * The textbook provides these transformations:
   * u = l + (r - l)(i + 0.5) / width
   * v = b + (t - b)(j + 0.5) / height
   *
   * We have two coordinates: The u, v, w, Camera Frame coordinates
   * and the image representation frame defined in the textbook.
   * u, v, w, camera frame coordinates vs pixel coordinates.
   * Remember that our defined coordinate system for (i, j) is defined
   * in the book starting from bottom left -> cartesian coordinate system
   *
   * Note that the image dimensions (i.e. l, r, b, t) are defined from the
   * camera frame coordinates
   *
   * formula r(t) = e + tv      We need to find direction vector v
   */

  // NOTE: The camera width and height is the dimensions of the camera frame.

  // Get the scalar (weights) for the uvw coordinates
  double u = -(camera.width / 2.0) + (camera.width * (j + 0.5) / width);
  double v = (camera.height / 2.0) - (camera.height * (i + 0.5) / height);
  double w = -camera.d;

  // Find the point <P> of (i, j) in the uvw camera frame coordinates
  // This is a projection onto the image plane
  Eigen::Vector3d P = camera.e + u * camera.u + v * camera.v + w * camera.w;

  // Set the viewing ray origin and direction
  ray.origin = camera.e;
  ray.direction = P - camera.e;
  ////////////////////////////////////////////////////////////////////////////
}



/**
 * NOTE: below is my original transformation, which is not necessarily wrong,
 * HOWEVER, these are the transformations from the textbook, which is mapping
 * (i, j) = (x, y) in cartesian coordinates. In reality (in c++), however,
 * (i, j) = (row, col). Hence, we have to switch i and j. FURTHER, in c++
 * we start in the upper left corner and work our way across and down, not
 * up, so we have to remove one of the negative signs!
 * This caused all my images to be sideways, and I coudn't figure this out
 * for the longest time since everything was working perfectly fine.
 */

  // uvw_point[0] = -(camera.width / 2.0) + (camera.width * (i + 0.5) / width);    // u
  // uvw_point[1] = -(camera.height / 2.0) + (camera.height * (j + 0.5) / height); // v
  // uvw_point[2] = camera.e[2] + -camera.d;                                       // w
