
#define PCL_NO_PRECOMPILE
#include <pcl/kdtree/kdtree_flann.h>
#include "las_io.h"
#include "computeEigen.h"
#include "computeCovMy.h"

int main (int argc, char** argv)
{
	string dir="D:\\data\\";
	string filename=dir+"LiDAR.las";
 
  pcl::PointCloud<MyPoint>::Ptr cloud (new pcl::PointCloud<MyPoint>);
  getCloudFromLas(filename,cloud,true);
  std::cout << "PointCloud has: " << cloud->points.size () << " data points." << std::endl; //*

  // Creating the KdTree object for the search method of the extraction

  pcl::KdTreeFLANN<MyPoint> kdtree;
  kdtree.setInputCloud (cloud);
  std::ofstream ofs("result.txt");
  for (int i=0;i<cloud->size();i++)
  {
	  std::vector<int> pointIdxSearch;
	  std::vector<float> pointSquaredDistance;
	  kdtree.radiusSearch(i,2.0,pointIdxSearch,pointSquaredDistance);
	  Eigen::Matrix3f cov;
	  Eigen::Vector4f xyz_centroid;
	  computeCovMy(*cloud,pointIdxSearch,cov,xyz_centroid); 	  
	  Eigen::Matrix3f eigenvector;
	  Eigen::Vector3f eigenvalue;
	  computeEigen(cov,eigenvector, eigenvalue);
	  ofs <<(*cloud)[i]<<" "<< eigenvalue.transpose()<<endl;
  } 
  ofs.close();
  system("pause");
  return (0);
}
