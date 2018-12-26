#pragma once


#include <flann/flann.hpp>
#include <memory>


template<typename T>
class pointCloud
{
public:
	pointCloud()=default;
	pointCloud(const std::vector<std::vector<T>> & input) { setinput(input); };
	~pointCloud() { delete[] data.ptr(); flann_index.reset(); }
	char *readline(FILE *input)
	{
		int len;
		if (fgets(line, max_line_len, input) == NULL)
			return NULL;
		while (strrchr(line, '\n') == NULL)
		{
			max_line_len *= 2;
			line = (char *)realloc(line, max_line_len);
			len = (int)strlen(line);
			if (fgets(line + len, max_line_len - len, input) == NULL)
				break;
		}
		len = (int)strlen(line);
		if (line[len - 1] == '\n')
			line[len - 1] = '\0';
		return line;
	}
	int loadfromfile(const char* filename) 
	{
		dim = 3;
		FILE *fp;
		int num_elements, i;
		char *endptr, *val, *next_val;
		double *line_data, value;
		fopen_s(&fp, filename, "r");
		if (fp == NULL)
		{
			printf("can't open input file %s\n", filename);
			return 0;
		}
		num_points = 0;
		max_line_len = 1024;
		line = (char *)malloc(max_line_len * sizeof(char));
		while (readline(fp) != NULL)
		{
			num_points++;
		}
		rewind(fp);
		data= flann::Matrix<T>(new T[dim*num_points * sizeof(T)], num_points, dim);
		line_data = (double *)malloc(max_line_len / 2 * sizeof(double));
		for (i = 0; i < num_points; i++)
		{
			readline(fp);
			val = strtok_s(line, " ", &next_val); // value1
			if (val == NULL)
			{
				printf("Empty line at line %d\n", i + 1);
				return 0;
			}
			line_data[0] = strtod(val, &endptr);// convert to double
			num_elements = 1;
			//while (1) //==================
			while(num_elements<3)
			{
				val = strtok_s(NULL, " ", &next_val); // value 2:end
				if (val == NULL)
					break;
				errno = 0;
				value = strtod(val, &endptr);// convert to double
				if (endptr == val || errno != 0 || *endptr != '\0')
				{
					printf("Wrong input format at line %d,column %d\n", num_points + 1, num_elements + 1);
					return 0;
				}
				line_data[num_elements] = value;
				num_elements++;
			}
			data[i][0]= line_data[0];
			data[i][1] = line_data[1];
			data[i][2] = line_data[2];
		}
		free(line);
		free(line_data);
		fclose(fp);
		return 1;

	};

	int setinput(const std::vector<std::vector<T>> & input) 
	{
		dim = 3;
		num_points = input.size();
		data = flann::Matrix<T>(new T[dim*num_points * sizeof(T)], num_points, dim);
		for (int i = 0; i < num_points; i++)
		{
			data[i][0] = input[i][0];
			data[i][1] = input[i][1];
			data[i][2] = input[i][2];
		}
		return num_points;
	};
	void initSearch() 
	{
		flann_index.reset(new flann::Index< flann::L2<T>>( data, flann::KDTreeSingleIndexParams(15))); // max 15 points/leaf
		flann_index->buildIndex();
	};
	int radiusSearch(flann::Matrix<T> &query, float radius,
		std::vector<int> &k_indices, std::vector<T> &k_sqr_distances) const
	{
		std::vector<std::vector<int> > indices(1);
		std::vector<std::vector<float> > dists(1);
		int neighbors_in_radius=flann_index->radiusSearch(query, indices, dists, radius, flann::SearchParams(128));
		k_indices = indices[0];
		k_sqr_distances = dists[0];
		return neighbors_in_radius;
	};
	int radiusSearch(int i, float radius,
		std::vector<int> &k_indices, std::vector<T> &k_sqr_distances) const
	{
		flann::Matrix<T> query(data[i], 1, dim);
		return radiusSearch(query, radius, k_indices, k_sqr_distances);
	};
	int radiusSearch(flann::Matrix<T> &queries, T r,
		std::vector<std::vector<int>> &k_indices, std::vector<std::vector<T>> &k_sqr_distances) const
	{
		return flann_index->radiusSearch(queries, k_indices, k_sqr_distances, r, flann::SearchParams(128));
	};
	int knnSearch(flann::Matrix<T> &query, int k,
		std::vector<int> &k_indices, std::vector<T> &k_sqr_distances) const
	{
		if (k > num_points)
			k = num_points;
		flann::Matrix<int> k_indices_mat(&k_indices[0], 1, k);
		flann::Matrix<float> k_sqr_distances_mat(&k_sqr_distances[0], 1, k);
		return flann_index->knnSearch(query, k_indices_mat, k_sqr_distances_mat, k, flann::SearchParams(128));
	};
	int knnSearch(int i, int k,
		std::vector<int> &k_indices, std::vector<T> &k_sqr_distances) const
	{
		flann::Matrix<T> query(data[i], 1, dim);
		return radiusSearch(query, radius, k_indices, k_sqr_distances);
	};
	int knnSearch(flann::Matrix<T> &queries , int k,
		std::vector<std::vector<int>> &k_indices, std::vector<std::vector<T>> &k_sqr_distances) const
	{
		if (k > num_points)
			k = num_points;
		return flann_index->knnSearch(queries, k_indices, k_sqr_distances, k, flann::SearchParams(128));
	};

	char *line;
	int max_line_len;
	int dim;
	int num_points;
	flann::Matrix<T> data;
	std::shared_ptr<flann::Index< flann::L2<T>>> flann_index;
	
};

