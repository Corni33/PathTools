/**
  PathTools.h
  Collection of datatypes and template functions for processing 2D points and paths 

  @version 0.1 2018/02/17
*/

#pragma once

#include <vector>
#include <cmath>
#include <limits>
#include <iostream>


namespace PathTools
{

/**
  Point 
*/
template <typename T>
struct Point
{
  T x, y;

  Point() : x(0), y(0) {} 
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

  PathPoint() : ind(0), u(0) {}
  PathPoint(size_t ind_, T u_) : ind(ind_), u(u_) {}
};

/**
  Path 
*/
template <typename T>
class Path
{
public:
  std::vector<Point<T>> points;
  std::vector<T> length;

  /**
  Pre-calculates the path length to each point of the path
  */
  void          prepare(void);

  /**
  Caclulates the coordinates of a PathPoint on this path
  @param PathPoint  
  @return Point containing the coordinates of the PathPoint
  */
  Point<T>      coordinates(PathPoint<T> const& pp);

  /**
  Caclulates the closest PathPoint on the path to Point p
  @param p: Point to which the closest PathPoint will be calculated
  @return closest PathPoint
  */
  PathPoint<T>  closest_point(Point<T> const& p);

  /**
  Moves a desired distance along the path and returns the resulting PathPoint
  @param start_point: PathPoint from which the movement will start
  @param dist_desired: Desired moving distance along the path
  @return PathPoint that is reached after the movement
  */
  PathPoint<T>  move_along(PathPoint<T> const& pp_start, T const dist_desired);
};

/**
  Calculates the squared distance between two 2D points
  @param p First point
  @param q Second point
  @return squared distance between p and q
*/
template <typename T>
T dist_2d_squared(Point<T> const& p, Point<T> const& q);

/**
Calculates the distance between two 2D points
@param p First point
@param q Second point
@return distance between p and q
*/
template <typename T>
T dist_2d(Point<T> const& p, Point<T> const& q);

/**
Calculates the index of the path point which is closest to a given point
@param path  
@param p Point, to which the closest path point should be calculated
@return Index of the path point which is closest to p
*/
template <typename T>
size_t closest_point_index(Path<T> const& path, Point<T> const& p);


/*****************************************************************************************************
Overloading of basic operators for data type "Point"
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
T dist_2d_squared(Point<T> const& p, Point<T> const& q)
{
  Point<T> tmp = p - q;
  return tmp.x*tmp.x + tmp.y*tmp.y;
}

template <typename T>
T dist_2d(Point<T> const& p, Point<T> const& q)
{
  return std::sqrt(dist_2d_squared(p, q));
}

template <typename T>
size_t closest_point_index(Path<T> const& path, Point<T> const& p)
{
  size_t closest_point_index = 0, ind = 0;
  T shortest_dist = std::numeric_limits<T>::max();
  T current_dist = T();

  for (Point<T> const& q : path.points)
  {
    current_dist = dist_2d_squared(p, q);
    if (dist_2d_squared(p, q) < shortest_dist)
    {
      shortest_dist = current_dist;
      closest_point_index = ind;
    }
    ++ind;
  }

  return closest_point_index;
}

template <typename T>
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
  
  if (u > 1)
  {
    u = static_cast<T>(1);
  }
  else if (u < 0)
  {
    u = static_cast<T>(0);
  }

  interp_factor = u;
  closest_point = a + u * v;  
  return;
}

template <typename T>
void Path<T>::prepare()
{
  this->length.resize(this->points.size());
  T dist = 0;

  this->length[0] = 0;

  for (size_t i = 0; i != this->points.size() - 1; ++i)
  {
    dist += dist_2d(this->points[i], this->points[i + 1]);
    this->length[i + 1] = dist;
  }

  for (auto const& l : this->length)
  {
    std::cout << "prepare: dist = " << l << std::endl;
  }
  return;
}

template <typename T>
Point<T> Path<T>::coordinates(PathPoint<T> const& pp)
{
  if (this->points.size() == 0 || pp.ind >= this->points.size())
  {
    return Point<T>();
  }

  if (this->points.size() == 1)
  {
    return this->points[0];
  }

  if (this->points.size() - 1 == pp.ind)
  {
    return this->points[pp.ind];
  }

  return this->points[pp.ind] + pp.u * (this->points[pp.ind + 1] - this->points[pp.ind]);
}

template <typename T>
PathPoint<T> Path<T>::closest_point(Point<T> const& p)
{
  PathPoint<T> pp;

  size_t const path_size = this->points.size();

  if (path_size <= 1)
  {    
    return pp; // ind=0, u=0
  }

  size_t ind_closest = closest_point_index(*this, p);
  std::cout << "ind closest: " << ind_closest << std::endl;
  Point<T> p_temp;

  if (ind_closest == 0)
  {
    pp.ind = 0;
    closest_point_on_line_segment(this->points[0], this->points[1], p, p_temp, pp.u);
    return pp;
  }

  if (ind_closest == path_size - 1)
  {    
    closest_point_on_line_segment(this->points[path_size - 2], this->points[path_size - 1], p, p_temp, pp.u);
    if (pp.u >= static_cast<T>(1.0))
    {
      pp.ind = path_size - 1;
    }
    else
    {
      pp.ind = path_size - 2;
    }
    return pp;
  }
  
  T u_next;

  closest_point_on_line_segment(this->points[ind_closest], this->points[ind_closest + 1], p, p_temp, u_next); // next segment
  T dist_next = dist_2d_squared(p, p_temp);

  closest_point_on_line_segment(this->points[ind_closest - 1], this->points[ind_closest], p, p_temp, pp.u); // previous segment
  T dist_prev = dist_2d_squared(p, p_temp);  
  
  if (dist_next < dist_prev)
  {
    pp.u = u_next;
    pp.ind = ind_closest;
  }
  else
  {
    pp.ind = ind_closest - 1;
  }

  return pp;
}

template <typename T>
PathPoint<T> Path<T>::move_along(PathPoint<T> const& pp_start, T const dist_desired)
{
  PathPoint<T> pp_target = pp_start;

  T length_seg = this->length[pp_start.ind + 1] - this->length[pp_start.ind];
  T dist_moved = (static_cast<T>(1.0) - pp_start.u) * length_seg;

  if (dist_moved > dist_desired) 
  {
    pp_target.u += dist_desired / length_seg;
    return pp_target;
  }

  T offset = this->length[pp_start.ind + 1] - dist_moved;
  T dist_desired_adapted = dist_desired + offset;

  for (size_t i = pp_start.ind; i != this->length.size(); ++i)
  {   
    if (this->length[i] > dist_desired_adapted)
    {
      T dist_seg_begin = this->length[i - 1];
      pp_target.u = (dist_desired_adapted - dist_seg_begin) / (this->length[i] - this->length[i - 1]);
      pp_target.ind = i - 1;
      return pp_target;
    }    
  }

  pp_target.u = 0;
  pp_target.ind = this->length.size() - 1;
  return pp_target;
}


}