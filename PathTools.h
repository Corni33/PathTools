/**
  PathTools.h
  Collection of datatypes and template functions for processing 2D points and paths 

  @version 0.1 2018/02/17
*/

#pragma once

#include <vector>
#include <cmath>
#include <limits>

namespace PathTools
{

/**
  2D Point 
*/
template <typename T>
struct Point2d
{
  Point2d() : x(T()), y(T()) {} // default constructor
  Point2d(T x_, T y_) : x(x_), y(y_) {}
  T x, y;
};

/**
  2D Path 
*/
template <typename T>
using Path2d = std::vector<Point2d<T>>;

/**
  Calculates the squared distance between two 2D points

  @param p First point
  @param q Second point
  @return squared distance between p and q
*/
template <typename T>
T dist_2d_squared(Point2d<T> const& p, Point2d<T> const& q);

/**
Calculates the (euclidean) distance between two 2D points

@param p First point
@param q Second point
@return distance between p and q
*/
template <typename T>
T dist_2d(Point2d<T> const& p, Point2d<T> const& q);

/**
Calculates the index of the path point which is closest to a given point

@param path  
@param p Point, to which the closest path point should be calculated
@return Index of the path point which is closest to p
*/
template <typename T>
size_t closest_point_index(Path2d<T> const& path, Point2d<T> const& p);

template <typename T>
T dist_to_line_segment(Point2d<T> const& a, Point2d<T> const& b, Point2d<T> const& p, Point2d<T>& closest_point);

template <typename T>
T closest_point_interpolated(Path2d<T> const& path, Point2d<T> const& p, const size_t closest_point_index, Point2d<T>& closest_point);

template <typename T>
T closest_point_interpolated(Path2d<T> const& path, Point2d<T> const& p, Point2d<T>& closest_point);



/*****************************************************************************************************
Implementation of template functions
******************************************************************************************************/

template <typename T>
Point2d<T> operator+ (Point2d<T> const& lhs, Point2d<T> const& rhs)
{
  return Point2d<T>(lhs.x + rhs.x, lhs.y + rhs.y);
}

template <typename T>
Point2d<T> operator- (Point2d<T> const& lhs, Point2d<T> const& rhs)
{
  return Point2d<T>(lhs.x - rhs.x, lhs.y - rhs.y);
}

template <typename T>
Point2d<T> operator* (T s, const Point2d<T>& p)
{
  return Point2d<T>(s * p.x, s * p.y);
}

template <typename T>
Point2d<T> operator* (const Point2d<T>& p, T s)
{
  return Point2d<T>(s * p.x, s * p.y);
}

template <typename T>
T dist_2d_squared(Point2d<T> const& p, Point2d<T> const& q)
{
  Point2d<T> tmp = p - q;
  return tmp.x*tmp.x + tmp.y*tmp.y;
}

template <typename T>
T dist_2d(Point2d<T> const& p, Point2d<T> const& q)
{
  return std::sqrt(dist_2d_squared(p, q));
}

template <typename T>
size_t closest_point_index(Path2d<T> const& path, Point2d<T> const& p)
{
  size_t closest_point_index = 0, ind = 0;
  T shortest_dist = std::numeric_limits<T>::max();
  T current_dist = T();

  for (Point2d<T> const& q : path)
  {
    current_dist = dist_2d_squared(p, q);
    if (current_dist < shortest_dist)
    {
      shortest_dist = current_dist;
      closest_point_index = ind;
    }
    ++ind;
  }

  return closest_point_index;
}

template <typename T>
T dist_to_line_segment(Point2d<T> const& a, Point2d<T> const& b, Point2d<T> const& p, Point2d<T>& closest_point)
{
  Point2d<T> v = b - a; // v is the vector pointing from a to b

  T dist_squared = dist_2d_squared(a, b);

  if (dist_squared < 0.00001) 
  {
    closest_point = a;
    return dist_2d(p, closest_point);
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

  closest_point = a + u * v;

  return dist_2d(p, closest_point);
}

template <typename T>
T closest_point_interpolated(Path2d<T> const& path, Point2d<T> const& p, Point2d<T>& closest_point)
{
  size_t path_size = path.size();

  if (path_size == 0)
  {
    return T();
  }

  size_t ind_closest = closest_point_index(path, p);

  if (path_size == 1) 
  {
    return dist_2d(path[0], p);
  }  

  if (ind_closest == 0)
  {
    return dist_to_line_segment(path[0], path[1], p, closest_point);
  }

  if (ind_closest == path_size - 1)
  {
    return dist_to_line_segment(path[path_size - 2], path[path_size - 1], p, closest_point);
  }

  Point2d<T> closest_point_next; // closest point to next line segment after path[ind_closest]

  T dist_next = dist_to_line_segment(path[ind_closest], path[ind_closest + 1], p, closest_point_next); // next segment
  T dist_prev = dist_to_line_segment(path[ind_closest - 1], path[ind_closest], p, closest_point); // previous segment

  if (dist_next < dist_prev)
  {
    closest_point = closest_point_next;
    return dist_next;
  }
  else
  {
    return dist_prev;
  }
}


/*template <typename T>
T look_ahead_point(Path2d<T> const& path, Point2d<T> const& p, Point2d<T>& closest_point)
{
  // TODO
}*/


}