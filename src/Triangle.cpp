#include "Triangle.h"
#include "Ray.h"
#include "Eigen/Dense"

bool Triangle::intersect(
  const Ray & ray, const double min_t, double & t, Eigen::Vector3d & n) const
{
  ////////////////////////////////////////////////////////////////////////////
  /**
   * NOTE: The following solution is found on page 79 in the textbook.
   * This is finally the fastest compile time of all my solutions below...
   * Just save every variable.
   */
  bool intersect = false;

  // Set up some variables
  Eigen::Vector3d P = std::get<0>(Triangle::corners);       // Point P
  Eigen::Vector3d Q = std::get<1>(Triangle::corners);       // Point Q
  Eigen::Vector3d R = std::get<2>(Triangle::corners);       // Point R
  Eigen::Vector3d PQ = Q - P;                               // direction vector
  Eigen::Vector3d PR = R - P;                               // direction vector
  Eigen::Vector3d v = ray.direction;
  Eigen::Vector3d O = ray.origin;

  // Get the normal of the two vectors in the plane of the triangle
  Eigen::Vector3d normal = PQ.cross(PR);

  // Check if the ray intersects with the plane made by the triangle
  if (v.dot(normal) == 0) {
    return false;
  }

  // If we make it here, then we know that the ray intersects the plane at some
  // point t
  t = (normal.dot(P - O)) / (normal.dot(v));
  if (t < min_t) {
    return false;
  }

  /**
   * Now we need to check if the intersection point is in the triangle
   * pg 79
   */
  // Eigen::Matrix3d A;
  // A << -PQ, -PR, v;
  // Eigen::Vector3d y = P - O;           // Ax = y, where x = <beta, gamma, t>

  // Set up the columns of matrix A
  // Column vector (P - Q)
  double a = -PQ[0];
  double b = -PQ[1];
  double c = -PQ[2];
  // Column vector (P - R)
  double d = -PR[0];
  double e = -PR[1];
  double f = -PR[2];
  // Column vector d (direction vector of ray)
  double g = v[0];
  double h = v[1];
  double i = v[2];

  // Column vector y                      Ax = y
  double j = P[0] - O[0];
  double k = P[1] - O[1];
  double l = P[2] - O[2];

  /**
   * Rather than use A.determinant(), we can reduce the number of operations
   * by reusing numbers.
   * double M = A.determinant();
   */
  double ei_minus_hf = e * i - h * f;
  double gf_minus_di = g * f - d * i;
  double dh_minus_eg = d * h - e * g;

  double ak_minus_jb = a * k - j * b;
  double jc_minus_al = j * c - a * l;
  double bl_minus_kc = b * l - k * c;

  // Get the determinant of A
  double M = a * ei_minus_hf + b * gf_minus_di + c * dh_minus_eg;

  // Solve for the values of x = <beta, gamma, t> in the system Ax = y
  double beta = (j * ei_minus_hf + k * gf_minus_di + l * dh_minus_eg) / M;
  double gamma = (i * ak_minus_jb + h * jc_minus_al + g * bl_minus_kc) / M;
  // double t_prime = -((f * ak_minus_jb + e * jc_minus_al + d * bl_minus_kc) / M);
  // assert(t == t_prime);                  // works properly

  if (beta >= 0 && gamma >= 0 && (beta + gamma) <= 1) {
    n = normal.normalized();
    intersect = true;
  }

  return intersect;
  ////////////////////////////////////////////////////////////////////////////
}




/**
 * NOTE: Solution 1 (below) **Works**
 * Solution follows the textbook solution, but still super slow to compile
 * the bunny.
 */
// #include "Triangle.h"
// #include "Ray.h"
// #include "Eigen/Dense"

// bool Triangle::intersect(
//   const Ray & ray, const double min_t, double & t, Eigen::Vector3d & n) const
// {
//   ////////////////////////////////////////////////////////////////////////////
//   /**
//    * NOTE: The following solution is found on page 79 in the textbook.
//    */
//   bool intersect = false;

//   // Set up some variables
//   Eigen::Vector3d a = std::get<0>(Triangle::corners);       // Point a
//   Eigen::Vector3d b = std::get<1>(Triangle::corners);       // Point b
//   Eigen::Vector3d c = std::get<2>(Triangle::corners);       // Point c
//   Eigen::Vector3d ab = b - a;                               // direction vector
//   Eigen::Vector3d ac = c - a;                               // direction vector
//   Eigen::Vector3d d = ray.direction;
//   Eigen::Vector3d e = ray.origin;

//   // Get the normal of the two vectors in the plane of the triangle
//   Eigen::Vector3d normal = ab.cross(ac);

//   // Check if the ray intersects with the plane made by the triangle
//   if (d.dot(normal) == 0) {
//     return false;
//   }

//   // If we make it here, then we know that the ray intersects the plane at some
//   // point t
//   t = (normal.dot(a - e)) / (normal.dot(d));
//   if (t < min_t) {
//     return false;
//   }

//   /**
//    * Now we need to check if the intersection point is in the triangle
//    * pg 79
//    */
//   Eigen::Matrix3d A;
//   A << -ab, -ac, d;
//   Eigen::Vector3d y = a - e;              // Ax = y, where x = <beta, gamma, t>

//   /**
//    * Rather than use A.determinant(), we can reduce the number of operations
//    * by reusing numbers.
//    * double M = A.determinant();
//    */
//   double ei_minus_hf = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
//   double gf_minus_di = A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2);
//   double dh_minus_eg = A(0, 1) * A(1, 2) - A(1, 1) * A(0, 2);

//   double ak_minus_jb = A(0, 0) * y[1] - y[0] * A(1, 0);
//   double jc_minus_al = y[0] * A(2, 0) - A(0, 0) * y[2];
//   double bl_minus_kc = A(1, 0) * y[2] - y[1] * A(2, 0);

//   // Get the determinant of A
//   double M = A(0, 0) * ei_minus_hf + A(1, 0) * gf_minus_di + A(2, 0) * dh_minus_eg;

//   // Solve for the values of x = <beta, gamma, t> in the system Ax = y
//   double beta = (y[0] * ei_minus_hf + y[1] * gf_minus_di + y[2] * dh_minus_eg) / M;
//   double gamma = (A(2, 2) * ak_minus_jb + A(1, 2) * jc_minus_al + A(0, 2) * bl_minus_kc) / M;
//   // double t_prime = -((A(2, 1) * ak_minus_jb + A(1, 1) * jc_minus_al + A(0, 1) * bl_minus_kc) / M);
//   // assert(t == t_prime);                  // works properly

//   if (beta >= 0 && gamma >= 0 && (beta + gamma) <= 1) {
//     n = normal.normalized();
//     intersect = true;
//   }

//   return intersect;
//   ////////////////////////////////////////////////////////////////////////////
// }





/**
 * NOTE: Solution 2 (below) **Works**
 * This is actually my prefered solution because the code is nicer by using
 * Eigen to solve the linear system Ax = b. However, this solution takes much
 * long to compile.
 */
// #include "Triangle.h"
// #include "Ray.h"
// #include "Eigen/Dense"

// bool Triangle::intersect(
//   const Ray & ray, const double min_t, double & t, Eigen::Vector3d & n) const
// {
//   ////////////////////////////////////////////////////////////////////////////
//   bool intersect = false;

//   // Set up some variables
//   Eigen::Vector3d P = std::get<1>(Triangle::corners);       // Point P
//   Eigen::Vector3d Q = std::get<0>(Triangle::corners);       // Point Q
//   Eigen::Vector3d R = std::get<2>(Triangle::corners);       // Point R
//   Eigen::Vector3d PQ = Q - P;
//   Eigen::Vector3d PR = R - P;
//   Eigen::Vector3d e = ray.origin;
//   Eigen::Vector3d v = ray.direction;

//   // Get the normal of the two vectors in the plane of the triangle
//   Eigen::Vector3d normal = PQ.cross(PR);

//   // Check if the ray intersects with the plane made by the triangle
//   if (v.dot(normal) == 0) {
//     return false;
//   }

//   // If we make it here, then we know that the ray intersects the plane at some
//   // point t
//   t = (normal.dot(P - e)) / (normal.dot(v));
//   if (t < min_t) {
//     return false;
//   }

//   /**
//    * Now we need to check if the intersection point is in the triangle
//    * https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
//    * How did I get these vector values? Look at ipad calculations
//    */
//   Eigen::Matrix3d A;
//   A << PQ, PR, -v; // put the column vectors into matrix A. Same as: A.col(0) = PQ, ...
//   Eigen::Vector3d b = e - P;
//   // Now we are going to solve Ax = b
//   Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);

//   // Check if a solution exists
//   bool a_solution_exists = (A * x).isApprox(b);
//   if (!a_solution_exists) {
//     return false;
//   }

//   // x contains our x[0] = beta, x[1] = gamma, x[2] = t values in that order
//   // NOTE: We only need one of these if, else if statements, but just to be safe...
//   if (x[0] < 0 || x[0] > 1 || x[1] < 0 || x[1] > 1) {
//     return false;

//   } else if (x[0] >= 0 && x[1] >= 0 && (x[0] + x[1]) <= 1) {
//     n = normal.normalized();
//     intersect = true;
//   }

//   // If we make it here, then t == x[2]
//   return intersect;
//   ////////////////////////////////////////////////////////////////////////////
// }

/**
   * NOTE: rounding error
   * If we make it here, then t = x[2]. However, since we are dealing with
   * doubles, there is some difference by the way that Eigen solves the linear
   * system. However, if we round to two decimal places, then t == x[2]
   */
  // assert(t == x[2]);
  // double t_round = ceilf(t * 100) / 100;
  // double x_idx2_round = ceilf(x[2] * 100) / 100;
  // assert(t_round == x_idx2_round);                  // This works


















/**
   * NOTE: Solution 3
   * I know the textbooks offers the solution by using baycentric
   * coordinates, but this direction (below) seems more intuitive to me and I also
   * just wanted to see it through.
   * There are a couple of ways to check if a point p_0 is in a triangle:
   * 1) Barycentric coordinates
   * 2) Sum the angle between the point and every vertex. It will only
   * equal 2 * pi iff p_0 is in the triangle
   * 3) Sum the area's of each vertex of the triangle with p_0 and the sum
   * will equal the area of the triangle PRQ iff p_0 is in the triangle.
   * Below I have provided one solution where I used area's of triangles.
   */


// #include "Triangle.h"
// #include "Ray.h"

// bool Triangle::intersect(
//   const Ray & ray, const double min_t, double & t, Eigen::Vector3d & n) const
// {
//   ////////////////////////////////////////////////////////////////////////////
//   bool intersect = false;

//   // Set up some variables
//   Eigen::Vector3d P = std::get<1>(Triangle::corners);       // Point P
//   Eigen::Vector3d Q = std::get<0>(Triangle::corners);       // Point Q
//   Eigen::Vector3d R = std::get<2>(Triangle::corners);       // Point R
//   Eigen::Vector3d PQ = Q - P;
//   Eigen::Vector3d PR = R - P;
//   Eigen::Vector3d e = ray.origin;
//   Eigen::Vector3d v = ray.direction;

//   // Get the normal of the two vectors in the plane of the triangle
//   Eigen::Vector3d normal = PQ.cross(PR);

//   // Check if the ray intersects with the plane made by the triangle
//   if (v.dot(normal) == 0) {
//     return false;
//   }

//   // If we make it here, then we know that the ray intersects the plane at some
//   // point t
//   t = (normal.dot(P - e)) / (normal.dot(v));
//   if (t < min_t) {
//     return false;
//   }

//   /**
//    * Now we need to check if the intersection point is in the triangle
//    * First find the point of intersection.
//    */
//   Eigen::Vector3d p1 = e + t * v;                 // point of intersection

//   /**
//    * Check if p1 is in the triangle.
//    * NOTE: remember that the norm of a vector gives the area of a paralellogram,
//    * so 1/2 of that is the area of the triangle.
//    */
//   double A = (1 / 2) * normal.norm();             // Area of PQR
//   Eigen::Vector3d pP = P - p1;
//   Eigen::Vector3d pQ = Q - p1;
//   Eigen::Vector3d pR = R - p1;
//   double A1 = (1 / 2) * pP.cross(pQ).norm();      // Area of pPQ
//   double A2 = (1 / 2) * pP.cross(pR).norm();      // Area of pPR
//   double A3 = (1 / 2) * pQ.cross(pR).norm();      // Area of pQR

//   if (A1 + A2 + A3 == A) {
//     // Then p is in the triangle PQR
//     intersect = true;
//     n = normal.normalized();
//   }


//   return intersect;
//   ////////////////////////////////////////////////////////////////////////////
// }
