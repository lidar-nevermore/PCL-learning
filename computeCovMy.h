
#ifndef COMPUTE_COV_MY_H_
#define COMPUTE_COV_MY_H_

#include <pcl/common/centroid.h>
template <typename PointT> inline unsigned int
	computeCovMy (const pcl::PointCloud<PointT> &cloud,
	const std::vector<int> &indices,						  
	Eigen::Matrix3f &covariance_matrix,
	Eigen::Vector4f &xyz_centroid)
{

	if (cloud.empty()||indices.empty())
		return (0);

	// Initialize to 0
	covariance_matrix.setZero ();
	unsigned point_count;
	float ave_x=0.0f,ave_y=0.0f,ave_z=0.0f;	
	float nor_x=0.0f,nor_y=0.0f,nor_z=0.0f;	
	if (cloud.is_dense)
	{
		point_count = indices.size ();
		for (std::vector<int>::const_iterator iIt = indices.begin (); iIt != indices.end (); ++iIt)
		{						
			ave_x += cloud[*iIt].x;			
			ave_y += cloud[*iIt].y;
			ave_z += cloud[*iIt].z;		

		}
		ave_x/=static_cast<float> (point_count);
		ave_y/=static_cast<float> (point_count);
		ave_z/=static_cast<float> (point_count);	
		for (std::vector<int>::const_iterator iIt = indices.begin (); iIt != indices.end (); ++iIt)
		{	
			nor_x=cloud[*iIt].x-ave_x;
			nor_y=cloud[*iIt].y-ave_y;
			nor_z=cloud[*iIt].z-ave_z;
			covariance_matrix (0, 0) += nor_x*nor_x;
			covariance_matrix (1, 1) += nor_y*nor_y;
			covariance_matrix (2, 2) += nor_z*nor_z;
			covariance_matrix (0, 1) += nor_x*nor_y;
			covariance_matrix (0, 2) += nor_x*nor_z;
			covariance_matrix (1, 2) += nor_y*nor_z;
		}
		covariance_matrix (1, 0) = covariance_matrix (0, 1);
		covariance_matrix (2, 0) = covariance_matrix (0, 2);
		covariance_matrix (2, 1) = covariance_matrix (1, 2);
		covariance_matrix/=static_cast<float> (point_count);
	}
	else
	{
		point_count = 0;
		for (std::vector<int>::const_iterator iIt = indices.begin (); iIt != indices.end (); ++iIt)
		{			
			if (!isFinite (cloud[*iIt]))
				continue;
			++point_count;
			ave_x += cloud[*iIt].x;
			ave_y += cloud[*iIt].y;
			ave_z += cloud[*iIt].z;			
		}
		ave_x/=static_cast<float> (point_count);
		ave_y/=static_cast<float> (point_count);
		ave_z/=static_cast<float> (point_count);
		for (std::vector<int>::const_iterator iIt = indices.begin (); iIt != indices.end (); ++iIt)
		{		
			if (!isFinite (cloud[*iIt]))
				continue;	
			nor_x=cloud[*iIt].x-ave_x;
			nor_y=cloud[*iIt].y-ave_y;
			nor_z=cloud[*iIt].z-ave_z;
			covariance_matrix (0, 0) += nor_x*nor_x;
			covariance_matrix (1, 1) += nor_y*nor_y;
			covariance_matrix (2, 2) += nor_z*nor_z;
			covariance_matrix (0, 1) += nor_x*nor_y;
			covariance_matrix (0, 2) += nor_x*nor_z;
			covariance_matrix (1, 2) += nor_y*nor_z;
		}
		covariance_matrix (1, 0) = covariance_matrix (0, 1);
		covariance_matrix (2, 0) = covariance_matrix (0, 2);
		covariance_matrix (2, 1) = covariance_matrix (1, 2);
		covariance_matrix/=static_cast<float> (point_count);
	}	
	xyz_centroid[0] = ave_x; 
	xyz_centroid[1] = ave_y; 
	xyz_centroid[2] = ave_z;
	xyz_centroid[3] = 0;
	return (point_count);
}

#endif  //#ifndef COMPUTE_COV_MY_H_
