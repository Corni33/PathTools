/**
  PathTools.h
  Collection of datatypes and templated functions for processing 2D points and paths 

  @version 0.2 2018/03/27
*/

#pragma once

#include <vector>
#include <cmath>
#include <iostream>

#include <Eigen/Core>


namespace PathTools
{

/**
  Point 
*/
template <typename T>
struct Point
{
  T x, y;

  Point() : x(static_cast<T>(0)), y(static_cast<T>(0)) {}
  Point(T x_, T y_) : x(x_), y(y_) {}  
};


/**
  Point on a Path
*/
template <typename T>
struct PathPoint
{
  size_t ind;
  T u;

  PathPoint() : ind(0), u(static_cast<T>(0)) {}
  PathPoint(size_t ind_, T u_) : ind(ind_), u(u_) {}
};

/**
  Path 
*/
template <typename T>
class Path
{
public:
  void add_point(const T x, const T y);
  void add_point(const Point<T> p); 

  const std::vector<T>& x_points() { return m_x; }
  const std::vector<T>& y_points() { return m_y; }

  /**
  Caclulates the x/y-coordinates of a "PathPoint" on this path
    @param PathPoint  
    @return Point containing the coordinates of the PathPoint
  */
  Point<T>      coordinates(PathPoint<T> const& pp);

  /**
  Caclulates the closest "PathPoint" to Point p
    @param p: Point to which the closest "PathPoint" will be calculated
    @return closest "PathPoint"
  */
  PathPoint<T>  closest_point(Point<T> const& p);

  /**
  Moves a desired distance along the path and returns the resulting PathPoint
    @param start_point: PathPoint from which the movement will start
    @param dist_desired: Desired moving distance along the path
    @return PathPoint that is reached after the movement
  */
  PathPoint<T>  move_along(PathPoint<T> const& pp_start, T const dist_desired);

private:
  std::vector<T> m_x, m_y, m_vx, m_vy, m_dst_sq, m_p_length, m_u;

};

/**
  Calculates the squared distance between two 2D points
  @param p First point
  @param q Second point
  @return squared distance between p and q
*/
template <typename T>
T dist_2d_squared(Point<T> const& p, Point<T> const& q);

template <typename T>
T dist_2d_squared(const T x0, const T y0, const T x1, const T y1);

/**
Calculates the distance between two 2D points
@param p First point
@param q Second point
@return distance between p and q
*/
template <typename T>
T dist_2d(Point<T> const& p, Point<T> const& q);

template <typename T>
T dist_2d(const T x0, const T y0, const T x1, const T y1);


/*****************************************************************************************************
Overloading of some basic operators 
******************************************************************************************************/

template <typename T>
Point<T> operator+ (Point<T> const& lhs, Point<T> const& rhs)
{
  return Point<T>(lhs.x + rhs.x, lhs.y + rhs.y);
}

template <typename T>
Point<T> operator- (Point<T> const& lhs, Point<T> const& rhs)
{
  return Point<T>(lhs.x - rhs.x, lhs.y - rhs.y);
}

template <typename T>
Point<T> operator* (T s, const Point<T>& p) // scalar product
{
  return Point<T>(s * p.x, s * p.y);
}

template <typename T>
Point<T> operator* (const Point<T>& p, T s)  // scalar product
{
  return Point<T>(s * p.x, s * p.y);
}

template <typename T>
std::ostream &operator<< (std::ostream& strm, const Point<T>& p) 
{ 
  return  strm << "(" << p.x << ", " << p.y << ")";
}

template <typename T>
std::ostream &operator<< (std::ostream& strm, const PathPoint<T>& pp)
{
  return  strm << "(ind=" << pp.ind << ", u=" << pp.u << ")";
}

/*****************************************************************************************************
Function implementations
******************************************************************************************************/

template <typename T>
void Path<T>::add_point(const T x, const T y)
{
  if (m_x.size() > 0)
  {
    T last_x = m_x.back();
    T last_y = m_y.back();

    auto ind_last = m_x.size() - 1;

    m_x.push_back(x);
    m_y.push_back(y);

    // vector pointing from last point to current point
    m_vx[ind_last] = x - last_x;
    m_vx.push_back(0);
    m_vy[ind_last] = y - last_y;
    m_vy.push_back(0);

    // squared distance between last point and current point
    m_dst_sq[ind_last] = dist_2d_squared(x, y, last_x, last_y);
    m_dst_sq.push_back(0);

    // distance from beginning of path until current point
    m_p_length[ind_last] += dist_2d(x, y, last_x, last_y);
    m_p_length.push_back(m_p_length.back());
  }
  else
  {
    m_x.push_back(x);
    m_y.push_back(y);
    m_vx.push_back(0);
    m_vy.push_back(0);
    m_dst_sq.push_back(0);
    m_p_length.push_back(0);
  }

  m_u.push_back(0);
}

template <typename T>
void Path<T>::add_point(const Point<T> p)
{ 
  add_point(p.x, p.y);
}

template <typename T>
T dist_2d_squared(Point<T> const& p, Point<T> const& q)
{
  return dist_2d_squared(p.x, p.y, q.x, q.y);
}

template <typename T>
T dist_2d_squared(const T x0, const T y0, const T x1, const T y1)
{
  return (x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0);
}

template <typename T>
T dist_2d(Point<T> const& p, Point<T> const& q)
{
  return std::sqrt(dist_2d_squared(p, q));
}

template <typename T>
T dist_2d(const T x0, const T y0, const T x1, const T y1)
{
  return std::sqrt(dist_2d_squared(x0, y0, x1, y1));
}

/*template <typename T>
void closest_point_on_line_segment(Point<T> const& a, Point<T> const& b, Point<T> const& p, Point<T>& closest_point, T& interp_factor)
{
  Point<T> v = b - a; 

  T dist_squared = dist_2d_squared(a, b);

  if (dist_squared < 0.00001) 
  {
    interp_factor = 0;
    closest_point = a;
    return;
  }
 
  T u = ((p.x - a.x) * v.x + (p.y - a.y) * v.y) / dist_squared;
  
  u = std::min<T>(1, std::max<T>(u, 0)); 

  interp_factor = u;
  closest_point = a + u * v;  
  return;
}*/

template <typename T>
Point<T> Path<T>::coordinates(PathPoint<T> const& pp)
{
  size_t n_points = m_x.size();

  if (n_points == 0 || pp.ind >= n_points)
  {
    return Point<T>(0, 0);
  }

  if (n_points == 1)
  {
    return Point<T>(m_x[0], m_y[0]);
  }

  if (n_points - 1 == pp.ind)
  {
    return Point<T>(m_x[pp.ind], m_y[pp.ind]);
  }

  return Point<T>(m_x[pp.ind] + pp.u * m_vx[pp.ind], m_y[pp.ind] + pp.u * m_vy[pp.ind]);
}

template <typename T>
PathPoint<T> Path<T>::closest_point(Point<T> const& p)
{ 
  using array_map = Eigen::Map<Eigen::Array<T, -1, 1>>;

  // wrap std::vector data into Eigen array maps
  auto sz = m_x.size();
  array_map x       (m_x.data(),      sz, 1);
  array_map y       (m_y.data(),      sz, 1);
  array_map vx      (m_vx.data(),     sz, 1);
  array_map vy      (m_vy.data(),     sz, 1);
  array_map dst_sq  (m_dst_sq.data(), sz, 1);
  array_map u       (m_u.data(),      sz, 1);
  
  // find minimum distance to every line segment of the path
  u = (((p.x - x) * vx + (p.y - y) * vy)).min(dst_sq).max(0);

  int ind;
  (Eigen::square(x + vx * u - p.x) + Eigen::square(y + vy * u - p.y)).minCoeff(&ind);

  return PathPoint<T>(ind, dst_sq[ind] > 0.00001 ? u[ind] / dst_sq[ind] : 0);
}

template <typename T>
PathPoint<T> Path<T>::move_along(PathPoint<T> const& pp_start, T const dist_desired)
{
  PathPoint<T> pp_target = pp_start;

  T length_seg = m_p_length[pp_start.ind + 1] - m_p_length[pp_start.ind];
  T dist_moved = (1.0 - pp_start.u) * length_seg;

  if (dist_moved > dist_desired) 
  {
    pp_target.u += dist_desired / length_seg;
    return pp_target;
  }

  T offset = m_p_length[pp_start.ind + 1] - dist_moved;
  T dist_desired_adapted = dist_desired + offset;

  size_t n_points = m_x.size();
  for (size_t i = pp_start.ind; i != n_points; ++i) 
  {   
    if (m_p_length[i] > dist_desired_adapted)
    {
      T dist_seg_begin = m_p_length[i - 1];
      pp_target.u = (dist_desired_adapted - dist_seg_begin) / (m_p_length[i] - m_p_length[i - 1]);
      pp_target.ind = i - 1;
      return pp_target;
    }    
  }

  pp_target.u = 0;
  pp_target.ind = n_points - 1;
  return pp_target;
}


}