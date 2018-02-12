#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include <queue>

//#include "Region_growing_my.cpp"

using namespace std;

const double M_PI=3.141592654;
string dir="D:/Lidar/LiDAR/regiongrowing/source/";//文件路径
string inputfile=dir+"nor.txt";
string point_nghbr_file=dir+"idx.txt";


inline bool	comparePair (std::pair<float, int> i, std::pair<float, int> j)
{
	return (i.first < j.first);
}
int readpoints(std::string filename,vector<vector<float>>& cloud,vector<vector<float>>& normals)
{
	ifstream ifs(inputfile);
	string str;	
	vector<float> pt(3);
	vector<float> pnor(4);		
	while (getline(ifs,str))
	{		
		stringstream ss(str);		
		while (ss>>pt[0]>>pt[1]>>pt[2]>>pnor[0]>>pnor[1]>>pnor[2]>>pnor[3])	
		{
			cloud.push_back(pt);
			normals.push_back(pnor);
		};
	}
	ifs.close();	
	return cloud.size();
}
bool validatePoint (
	int initial_seed, 
	int point, 
	int nghbr, 
	bool& is_a_seed,
	const vector<vector<float>>& normals_,
	float theta_threshold,
	float curvature_threshold) 
{
	is_a_seed = true;
	float cosine_threshold = cosf (theta_threshold);	

	//check the angle between normals	
	float dot_product = fabsf (normals_[point][0]*normals_[nghbr][0]+normals_[point][1]*normals_[nghbr][1]+normals_[point][2]*normals_[nghbr][2]);
	if (dot_product < cosine_threshold)
	{
		return (false);
	}	
	// check the curvature
	if (normals_[nghbr][3] > curvature_threshold)
	{
		is_a_seed = false;
	}

	return (true);
}
int growRegion (int initial_seed, int segment_number,
	const vector<vector<float>>& normals_,
	const vector<vector<int>>& point_neighbours,
	vector<int>& point_labels_,
	vector<int> & num_pts_in_segment_,
	float theta,
	float curvature)
{
	std::queue<int> seeds;
	seeds.push (initial_seed);
	point_labels_[initial_seed] = segment_number;

	int num_pts_in_segment = 1;

	while (!seeds.empty ())
	{
		int curr_seed;
		curr_seed = seeds.front ();		
		seeds.pop ();		
		size_t i_nghbr = 0;
		while ( i_nghbr < point_neighbours[curr_seed].size () )
		{
			int index = point_neighbours[curr_seed][i_nghbr];
			if (point_labels_[index] != -1)
			{
				i_nghbr++;
				continue;
			}

			bool is_a_seed = false;
			bool belongs_to_segment = validatePoint (initial_seed, curr_seed, index, is_a_seed,normals_,theta,curvature);

			if (belongs_to_segment == false)
			{
				i_nghbr++;
				continue;
			}

			point_labels_[index] = segment_number;
			num_pts_in_segment++;

			if (is_a_seed)
			{
				seeds.push (index);
			}

			i_nghbr++;
		}// next neighbour
	}// next seed

	return (num_pts_in_segment);
}


int getPointNeighbours (int num,std::string filename, vector<vector<int>>& point_neighbours)
{
	point_neighbours.clear();
	std::vector<int> neighbours;		
	ifstream ifs(filename);
	string str;
	point_neighbours.resize (num, neighbours);	
	for (int i_point = 0; i_point < num; i_point++)
	{
		int data=0;
		neighbours.clear ();
		getline(ifs,str);
		stringstream ss(str);
		while (ss>>data)
			neighbours.push_back(data);
		//search_->nearestKSearch (i_point, neighbour_number_, neighbours, distances);
		point_neighbours[i_point].swap (neighbours);
	}	
	ifs.close();
	return point_neighbours.size();
}
void applySmoothRegionGrowingAlgorithm (	
	const vector<vector<float>>& normals_,
	const vector<vector<int>>& point_neighbours,
	vector<int>& point_labels_,
	vector<int> & num_pts_in_segment_,
	float theta,
	float curvature	)
{
	int num_of_pts = static_cast<int> (normals_.size ());
	point_labels_.resize (num_of_pts, -1);

	std::vector< std::pair<float, int> > point_residual;
	std::pair<float, int> pair;
	point_residual.resize (num_of_pts, pair);

	
	for (int i_point = 0; i_point < num_of_pts; i_point++)
	{			
		point_residual[i_point].first = normals_[i_point][3];
		point_residual[i_point].second = i_point;
	}
	std::sort (point_residual.begin (), point_residual.end (), comparePair);
	
	int seed_counter = 0;
	int seed = point_residual[seed_counter].second;

	int segmented_pts_num = 0;
	int number_of_segments = 0;
	while (segmented_pts_num < num_of_pts)
	{
		int pts_in_segment;		
		pts_in_segment = growRegion (seed, number_of_segments,normals_,point_neighbours,point_labels_,num_pts_in_segment_,theta,curvature);
		segmented_pts_num += pts_in_segment;
		num_pts_in_segment_.push_back (pts_in_segment);
		number_of_segments++;

		//find next point that is not segmented yet
		for (int i_seed = seed_counter + 1; i_seed < num_of_pts; i_seed++)
		{
			int index = point_residual[i_seed].second;
			if (point_labels_[index] == -1)
			{
				seed = index;
				seed_counter = i_seed;
				break;
			}
		}
	}
}


void assembleRegions (	
	const vector<vector<int>>& point_neighbours,
	vector<int> & num_pts_in_segment_,
	vector<int>& point_labels_,
	vector<vector<int>>&clusters,
	int min_pts_per_cluster_,
	int max_pts_per_cluster_)
{
	int number_of_segments = static_cast<int> (num_pts_in_segment_.size ());
	int number_of_points = static_cast<int> (point_labels_.size ());

	vector<int> segment;
	vector<vector<int>> clusters_(number_of_segments, segment);
	for (int i_seg = 0; i_seg < number_of_segments; i_seg++)
	{
		clusters_[i_seg].resize ( num_pts_in_segment_[i_seg], 0);
	}
	std::vector<int> counter;
	counter.resize (number_of_segments, 0);
	for (int i_point = 0; i_point < number_of_points; i_point++)
	{
		int segment_index = point_labels_[i_point];
		if (segment_index != -1)
		{
			int point_index = counter[segment_index];
			clusters_[segment_index][point_index] = i_point;
			counter[segment_index] = point_index + 1;
		}
	}
	clusters.resize (clusters_.size ());
	for (int i=0;i<point_labels_.size();i++)
	{
		point_labels_[i]=0;
	}
	auto cluster_iter_input = clusters.begin ();
	int num=0;
	for (auto cluster_iter = clusters_.begin (); cluster_iter != clusters_.end (); cluster_iter++)
	{
		if ((static_cast<int> (cluster_iter->size ()) >= min_pts_per_cluster_) &&
			(static_cast<int> (cluster_iter->size ()) <= max_pts_per_cluster_))
		{
			*cluster_iter_input = *cluster_iter;
			cluster_iter_input++;
			num++;
			for (auto it=cluster_iter->begin();it!=cluster_iter->end();it++)
			{
				point_labels_[*it]=num;
			}
		}
	}
	
	clusters.resize(num);
	
}

int main (int argc, char** argv)
{	
	vector <vector<int>> clusters;
	vector<vector<int>> point_neighbours_;
	vector<int> num_pts_in_segment_;
	vector<int> point_labels_;
	vector<vector<float>> input_;
	vector<vector<float>> normals_;	
	float theta_threshold_=10.0 / 180.0 * M_PI;
	float curvature_threshold_=0.03;
	int min_pts_per_cluster_=50;
	int max_pts_per_cluster_=20000;

	int num_of_pts=readpoints(inputfile,input_,normals_);
	getPointNeighbours(num_of_pts,point_nghbr_file,point_neighbours_);
	applySmoothRegionGrowingAlgorithm(normals_,point_neighbours_,point_labels_,num_pts_in_segment_,theta_threshold_,curvature_threshold_);
	assembleRegions(point_neighbours_,num_pts_in_segment_,point_labels_,clusters,min_pts_per_cluster_,max_pts_per_cluster_);

	ofstream ofs(dir+"result.txt");
	for (int index=0;index<input_.size();index++)
	{
		ofs <<input_[index][0]  <<" "
			<<input_[index][1]  <<" "
			<<input_[index][2]  <<" "			
			<<point_labels_[index]  <<endl;
	}

	ofs.close();
	
	system("pause");
	return (0);
}
