#include "laser_line_extraction/line_extraction.h"
#include <algorithm>
#include <Eigen/Dense>
#include <iostream>
#include <ros/ros.h>

namespace line_extraction
{

///////////////////////////////////////////////////////////////////////////////
// Constructor / destructor
///////////////////////////////////////////////////////////////////////////////
LineExtraction::LineExtraction()
{
}

LineExtraction::~LineExtraction()
{
}

///////////////////////////////////////////////////////////////////////////////
// Main run function
///////////////////////////////////////////////////////////////////////////////
void LineExtraction::extractLines(std::vector<Line>& lines) 
{
  // Resets
  filtered_indices_ = c_data_.indices;
  lines_.clear();

  // Filter indices
  filterCloseAndFarPoints();
  outlier_indices_.clear();
  filterOutlierPoints();
  outlier_number_ = outlier_indices_.size();

  // // Return no lines if not enough points left
  if (filtered_indices_.size() <= std::max(params_.min_line_points, static_cast<unsigned int>(3)))
  {
    return;
  }

  // // Split indices into lines and filter out short and sparse lines
  split(filtered_indices_);
  filterLines();
  // checkOutliers();

  // Fit each line using least squares and merge colinear lines
  for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it)
  {
    it->leastSqFit();
  }
  
  // If there is more than one line, check if lines should be merged based on the merging criteria
  if (lines_.size() > 1)
  {
    mergeLines();
  }

  lines = lines_;
}

///////////////////////////////////////////////////////////////////////////////
// Data setting
///////////////////////////////////////////////////////////////////////////////
void LineExtraction::setCachedData(const std::vector<double>& bearings,
                                   const std::vector<double>& cos_bearings,
                                   const std::vector<double>& sin_bearings,
                                   const std::vector<unsigned int>& indices)
{
  c_data_.bearings = bearings;
  c_data_.cos_bearings = cos_bearings;
  c_data_.sin_bearings = sin_bearings;
  c_data_.indices = indices;
}

void LineExtraction::setRangeData(const std::vector<double>& ranges)
{
  r_data_.ranges = ranges;
  r_data_.xs.clear();
  r_data_.ys.clear();
  for (std::vector<unsigned int>::const_iterator cit = c_data_.indices.begin(); 
       cit != c_data_.indices.end(); ++cit)
  {
    r_data_.xs.push_back(c_data_.cos_bearings[*cit] * ranges[*cit]);
    r_data_.ys.push_back(c_data_.sin_bearings[*cit] * ranges[*cit]);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Parameter setting
///////////////////////////////////////////////////////////////////////////////
void LineExtraction::setBearingVariance(double value)
{
  params_.bearing_var = value;
}

void LineExtraction::setRangeVariance(double value)
{
  params_.range_var = value;
}

void LineExtraction::setLeastSqAngleThresh(double value)
{
  params_.least_sq_angle_thresh = value;
}

void LineExtraction::setLeastSqRadiusThresh(double value)
{
  params_.least_sq_radius_thresh = value;
}

void LineExtraction::setMaxLineGap(double value)
{
  params_.max_line_gap = value;
}

void LineExtraction::setMinLineLength(double value)
{
  params_.min_line_length = value;
}

void LineExtraction::setMinLinePoints(unsigned int value)
{
  params_.min_line_points = value;
}

void LineExtraction::setMinRange(double value)
{
  params_.min_range = value;
}

void LineExtraction::setMaxRange(double value)
{
  params_.max_range = value;
}

void LineExtraction::setMinSplitDist(double value)
{
  params_.min_split_dist = value;
}

void LineExtraction::setOutlierDist(double value)
{
  params_.outlier_dist = value;
}

///////////////////////////////////////////////////////////////////////////////
// Utility methods
///////////////////////////////////////////////////////////////////////////////
double LineExtraction::chiSquared(const Eigen::Vector2d &dL, const Eigen::Matrix2d &P_1,
                                  const Eigen::Matrix2d &P_2)
{
  return dL.transpose() * (P_1 + P_2).inverse() * dL;
}

double LineExtraction::distBetweenPoints(unsigned int index_1, unsigned int index_2)
{
  return sqrt(pow(r_data_.xs[index_1] - r_data_.xs[index_2], 2) + 
              pow(r_data_.ys[index_1] - r_data_.ys[index_2], 2));
}

///////////////////////////////////////////////////////////////////////////////
// Filtering points
///////////////////////////////////////////////////////////////////////////////
void LineExtraction::filterCloseAndFarPoints()
{
  std::vector<unsigned int> output;
  for (std::vector<unsigned int>::const_iterator cit = filtered_indices_.begin(); 
       cit != filtered_indices_.end(); ++cit)
  {
    const double& range = r_data_.ranges[*cit];
    if (range >= params_.min_range && range <= params_.max_range)
    {
      output.push_back(*cit);
    }
  }
  filtered_indices_ = output;
}

void LineExtraction::filterOutlierPoints()
{
  if (filtered_indices_.size() < 3)
  {
    return;
  }

  std::vector<unsigned int> output;
  unsigned int p_i, p_j, p_k;
  for (std::size_t i = 0; i < filtered_indices_.size(); ++i)
  {

    // Get two closest neighbours

    p_i = filtered_indices_[i];
    if (i == 0) // first point
    {
      p_j = filtered_indices_[i + 1];
      p_k = filtered_indices_[i + 2];
    }
    else if (i == filtered_indices_.size() - 1) // last point
    {
      p_j = filtered_indices_[i - 1];
      p_k = filtered_indices_[i - 2];
    }
    else // middle points
    {
      p_j = filtered_indices_[i - 1];
      p_k = filtered_indices_[i + 1];
    }

    // Check if point is an outlier

    if (fabs(r_data_.ranges[p_i] - r_data_.ranges[p_j]) > params_.outlier_dist &&
        fabs(r_data_.ranges[p_i] - r_data_.ranges[p_k]) > params_.outlier_dist) 
    {
      // Check if it is close to line connecting its neighbours
      std::vector<unsigned int> line_indices;
      line_indices.push_back(p_j);
      line_indices.push_back(p_k);
      Line line(c_data_, r_data_, params_, line_indices);
      line.endpointFit();
      if (line.distToPoint(p_i) > params_.min_split_dist)
      {
        // 
        if (r_data_.ranges[p_i] < 30)
        {
          // ROS_WARN_STREAM(r_data_.ranges[p_i]);
          outlier_indices_.push_back(p_i);
        }
        //
        continue; // point is an outlier
      }
    }

    output.push_back(p_i);
  }

  filtered_indices_ = output;
}

///////////////////////////////////////////////////////////////////////////////
// Filtering and merging lines
///////////////////////////////////////////////////////////////////////////////
void LineExtraction::filterLines()
{
  std::vector<Line> output;
  for (std::vector<Line>::const_iterator cit = lines_.begin(); cit != lines_.end(); ++cit)
  {
    if (cit->length() >= params_.min_line_length && cit->numPoints() >= params_.min_line_points)
    {
      output.push_back(*cit);
    }
    else
    {
      for (int i=0; i < cit->indices_.size(); i++)
      {
        outlier_indices_.push_back(cit->indices_[i]);
      }
    }
  }
  lines_ = output;
}

void LineExtraction::mergeLines()
{
  std::vector<Line> merged_lines;

  for (std::size_t i = 1; i < lines_.size(); ++i)
  {
    // Get L, P_1, P_2 of consecutive lines
    Eigen::Vector2d L_1(lines_[i-1].getRadius(), lines_[i-1].getAngle());
    Eigen::Vector2d L_2(lines_[i].getRadius(), lines_[i].getAngle());
    Eigen::Matrix2d P_1;
    P_1 << lines_[i-1].getCovariance()[0], lines_[i-1].getCovariance()[1],
           lines_[i-1].getCovariance()[2], lines_[i-1].getCovariance()[3];
    Eigen::Matrix2d P_2;
    P_2 << lines_[i].getCovariance()[0], lines_[i].getCovariance()[1],
           lines_[i].getCovariance()[2], lines_[i].getCovariance()[3];

    // Merge lines if chi-squared distance is less than 3
    if (chiSquared(L_1 - L_2, P_1, P_2) < 3)
    {
      // Get merged angle, radius, and covariance
      Eigen::Matrix2d P_m = (P_1.inverse() + P_2.inverse()).inverse();
      Eigen::Vector2d L_m = P_m * (P_1.inverse() * L_1 + P_2.inverse() * L_2);
      // Populate new line with these merged parameters
      boost::array<double, 4> cov;
      cov[0] = P_m(0,0);
      cov[1] = P_m(0,1);
      cov[2] = P_m(1,0);
      cov[3] = P_m(1,1);
      std::vector<unsigned int> indices;
      const std::vector<unsigned int> &ind_1 = lines_[i-1].getIndices();
      const std::vector<unsigned int> &ind_2 = lines_[i].getIndices();
      indices.resize(ind_1.size() + ind_2.size());
      indices.insert(indices.end(), ind_1.begin(), ind_1.end());
      indices.insert(indices.end(), ind_2.begin(), ind_2.end());
      Line merged_line(L_m[1], L_m[0], cov, lines_[i-1].getStart(), lines_[i].getEnd(), indices);
      // Project the new endpoints
      merged_line.projectEndpoints();
      lines_[i] = merged_line;
    }
    else
    {
      merged_lines.push_back(lines_[i-1]);
    }

    if (i == lines_.size() - 1)
    {
      merged_lines.push_back(lines_[i]);
    }
  }
  lines_ = merged_lines;
}

///////////////////////////////////////////////////////////////////////////////
// Splitting points into lines
///////////////////////////////////////////////////////////////////////////////
void LineExtraction::split(const std::vector<unsigned int>& indices)
{
  // Don't split if only a single point (only occurs when orphaned by gap)
  if (indices.size() <= 1)
  {
    return;
  }

  Line line(c_data_, r_data_, params_, indices);
  line.endpointFit();
  double dist_max = 0;
  double gap_max = 0;
  double dist, gap;
  int i_max, i_gap;
  
  // Find the farthest point and largest gap
  for (std::size_t i = 1; i < indices.size() - 1; ++i)
  {
    dist = line.distToPoint(indices[i]);
    if (dist > dist_max)
    {
      dist_max = dist;
      i_max = i;
    }
    gap = distBetweenPoints(indices[i], indices[i+1]);
    if (gap > gap_max)
    {
      gap_max = gap;
      i_gap = i;
    }
  }

  // Check for gaps at endpoints
  double gap_start = distBetweenPoints(indices[0], indices[1]);
  if (gap_start > gap_max)
  {
    gap_max = gap_start;
    i_gap = 1;
  }
  double gap_end = distBetweenPoints(indices.rbegin()[1], indices.rbegin()[0]);
  if (gap_end > gap_max)
  {
    gap_max = gap_end;
    i_gap = indices.size() - 1;
  }

  // Check if line meets requirements or should be split
  if (dist_max < params_.min_split_dist && gap_max < params_.max_line_gap)
  {
    lines_.push_back(line);
  }
  else
  {
    int i_split = dist_max >= params_.min_split_dist ? i_max : i_gap;
    std::vector<unsigned int> first_split(&indices[0], &indices[i_split - 1]);
    std::vector<unsigned int> second_split(&indices[i_split], &indices.back());
    split(first_split);
    split(second_split);
  }

}

void LineExtraction::checkOutliers()
{
  std::vector<double> kickOutList;
  for(int i=0; i < lines_.size(); i++)
  {
     // Go through all lines to see if they can be extended
    if(outlier_indices_.size() == 0)
    {
      continue;
    }
    for(int j=0;  j < outlier_indices_.size(); j++)
    {
      // For each line check with all outliers to see if they are within range
      double pointToLineDist = pointToLinePerpendicularDist(i, outlier_indices_[j]);
      double pointToEndDist = pointToEndOfLineDistance(i, outlier_indices_[j]);

      if (pointToLineDist < 0.05 && 
          pointToEndDist < 0.2)
      {
        // If within tolerances, see if it replaces start or end point
        fitIntoLine(i, outlier_indices_[j]);
        kickOutList.push_back(j);
      }
    }
    for (int k = kickOutList.size() -1; k >= 0; k--)
    {
      if(kickOutList[k] > outlier_indices_.size()-1)
      {
        continue;
      }
      outlier_indices_.erase(outlier_indices_.begin()+kickOutList[k]);
    }
  }
  
  if (outlier_number_ != outlier_indices_.size())
  {
    outlier_number_ = outlier_indices_.size();
    ROS_WARN_STREAM("iteration");
    checkOutliers();
  }
}

double LineExtraction::pointToLinePerpendicularDist(int line_index, int outlier_index)
{
  double a = lines_[line_index].getStart()[1] - lines_[line_index].getEnd()[1];
  double b = lines_[line_index].getEnd()[0] - lines_[line_index].getStart()[0];
  double c = lines_[line_index].getStart()[0]*lines_[line_index].getEnd()[1]-lines_[line_index].getEnd()[0]*lines_[line_index].getStart()[1];

  return ( std::abs(a*r_data_.xs[outlier_index] + b*r_data_.ys[outlier_index] + c) / std::sqrt(std::pow(a, 2) + std::pow(b, 2)) );
}

double LineExtraction::pointToEndOfLineDistance(int line_index, int outlier_index)
{
  double distToStart = std::sqrt(std::pow(lines_[line_index].getStart()[1] - r_data_.ys[outlier_index], 2) 
     + std::pow(lines_[line_index].getStart()[0] - r_data_.xs[outlier_index], 2));
  
  double distToEnd = std::sqrt(std::pow(lines_[line_index].getEnd()[1] - r_data_.ys[outlier_index], 2) 
     + std::pow(lines_[line_index].getEnd()[0] - r_data_.xs[outlier_index], 2));

  return distToStart > distToEnd ? distToEnd : distToStart;
}

void LineExtraction::fitIntoLine(int line_index, int outlier_index)
{
  std::vector<double> lineMidpoint (2);
  lineMidpoint[0] = (lines_[line_index].getStart()[0] + lines_[line_index].getEnd()[0]) / 2;
  lineMidpoint[1] = (lines_[line_index].getStart()[1] + lines_[line_index].getEnd()[1]) / 2;
  
  double midToStartDist = std::sqrt(std::pow(lines_[line_index].getStart()[1] - lineMidpoint[1], 2) 
                          + std::pow(lines_[line_index].getStart()[0] - lineMidpoint[0], 2));
  
  double midToEndDist   = std::sqrt(std::pow(lines_[line_index].getEnd()[1] - lineMidpoint[1], 2) 
                          + std::pow(lines_[line_index].getEnd()[0] - lineMidpoint[0], 2));

  double outlierToStartDist = std::sqrt(std::pow(lines_[line_index].getStart()[1] - r_data_.ys[outlier_index], 2) 
                              + std::pow(lines_[line_index].getStart()[0] - r_data_.xs[outlier_index], 2));

  double outlierToEndDist   = std::sqrt(std::pow(lines_[line_index].getEnd()[1] - r_data_.ys[outlier_index], 2) 
                              + std::pow(lines_[line_index].getEnd()[0] - r_data_.xs[outlier_index], 2));
  
  double outlierToMidDist = std::sqrt(std::pow(r_data_.ys[outlier_index] - lineMidpoint[1], 2) 
                          + std::pow(r_data_.xs[outlier_index] - lineMidpoint[0], 2));

  if (outlierToStartDist < outlierToEndDist)
  {
    if(outlierToMidDist > midToStartDist)
    {
      lines_[line_index].indices_.insert(lines_[line_index].indices_.begin(), outlier_index);
      ROS_WARN_STREAM("start point" << lines_[line_index].getStart()[0] << " " << lines_[line_index].getStart()[1]);
      lines_[line_index].endpointFit();
      ROS_WARN_STREAM("start point" << lines_[line_index].getStart()[0] << " " << lines_[line_index].getStart()[1]);
      ROS_WARN_STREAM("changed start" << lines_[line_index].indices_.size() << " " << line_index);
    }
  }
  else
  {
    if(outlierToMidDist > midToEndDist)
    {
      lines_[line_index].indices_.push_back(outlier_index);
      ROS_WARN_STREAM("end point" << lines_[line_index].getEnd()[0] << " " << lines_[line_index].getEnd()[1]);
      lines_[line_index].endpointFit();
      ROS_WARN_STREAM("end point" << lines_[line_index].getEnd()[0] << " " << lines_[line_index].getEnd()[1]);
      ROS_WARN_STREAM("changed end" << lines_[line_index].indices_.size() << " " << line_index);
    }
  }
}

} // namespace line_extraction
