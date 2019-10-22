#include <flann/flann.hpp>
#include <flann/io/hdf5.h>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <string>

using namespace flann;
typedef double datatype;
static double taketime()
{
	return (double)(clock()) / CLOCKS_PER_SEC;
}

void exit_with_help()
{
	printf(
		"Usage: radiusDensity [options] input_h5_file density_file\n"
		"options:\n"
		"-h help\n"
		"-p points : set points name in input h5file (default pts)\n"		
		"-q query : set query name in input h5file (default pts)\n"		
		"-c cores : num of cpu cores to use, 0 for auto (default 0)\n"		
		"-r radius : set the radius of search (default 1)\n"
		"-n neighbor file : save radiusSearch neighbor results to file (default NULL)\n"	
		"-d distance file : save radiusSearch distance results to file (default NULL)\n\n"			
	);
	exit(1);
}

template <typename T, class function>
void savefile(const std::vector<std::vector<T>> &vec, const char *filename, const function& func)
{
	T value;
	FILE *fp;	
	fopen_s(&fp, filename, "w");
	if (fp == NULL)
	{
		printf("can't open output file %s\n", filename);
		return;
	}
	for (int i = 0; i < vec.size(); ++i)
	{
		for (int j = 0; j < vec[i].size(); ++j)
		{
			fprintf(fp, "%g ", func(vec[i][j]));
		}
		fprintf(fp, "\n");
	}	
	fclose(fp);
	return;
}
class parameters
{
public:
	parameters() :pointsName("pts"), queryName("pts"), 
		nghbrFilename(), distFilename(), hdf5Filename(), densityFilename(),
		cores(0), radius(1) {};
	std::string pointsName;
	std::string queryName;
	std::string nghbrFilename;	
	std::string distFilename;
	std::string hdf5Filename;
	std::string densityFilename;
	size_t cores;
	float radius;
} param;
void parse_command_line(int argc, char **argv)
{
	int i;
	// parse options
	for (i = 1; i < argc; i++)
	{
		if (argv[i][0] != '-') break;
		if (++i >= argc)
			exit_with_help();
		switch (argv[i - 1][1])
		{
		case 'p':
			param.pointsName = argv[i];
			break;
		case 'q':
			param.queryName = argv[i];
			break;
		case 'n':
			param.nghbrFilename = argv[i];
			break;
		case 'd':
			param.distFilename = argv[i];
			break;
		case 'r':
			param.radius = std::stof(argv[i]);
			break;
		case 'c':
			param.cores = std::stoi(argv[i]);
			break;		
		default:
			fprintf(stderr, "Unknown option: -%c\n", argv[i - 1][1]);
			exit_with_help();
		}
	}
	// determine filenames

	if (i >= argc)
		exit_with_help();	
	param.hdf5Filename = argv[i];

	if (i < argc - 1)
		param.densityFilename=argv[i + 1];
	else
	{
		char *p = strrchr(argv[i], '/');
		if (p == NULL)
			p = argv[i];
		else
			++p;
		param.densityFilename.assign(p);	
	}
}
int main(int argc, char** argv)
{
	parse_command_line(argc, argv);
	double start_time = taketime();
	Matrix<datatype> dataset;
	Matrix<datatype> query;
	load_from_file(dataset, param.hdf5Filename, param.pointsName);
	load_from_file(query, param.hdf5Filename, param.queryName);

	std::cout <<"dataset rows:"<< dataset.rows << " ,cols:" << dataset.cols << " " << std::endl;
	std::cout << "query rows:" << dataset.rows << " ,cols:" << dataset.cols << " " << std::endl;
	std::vector< std::vector<size_t> > indices;
	std::vector<std::vector<datatype> > dists;

	// construct an randomized kd-tree index using 4 kd-trees
	flann::Index<flann::L2_Simple<datatype> > flann_index(dataset, flann::KDTreeSingleIndexParams(64));
	flann_index.buildIndex();
	flann::SearchParams searchPara(128);
	searchPara.cores = param.cores;
	// do a knn search, using 128 checks
	flann_index.radiusSearch(query, indices, dists, param.radius,searchPara);
	printf("compution time: %f sec\n", taketime() - start_time);

	if (!param.nghbrFilename.empty())
		savefile(indices, param.nghbrFilename.c_str(), [](double a) {return a + 1; });
	if (!param.distFilename.empty())
		savefile(dists, param.distFilename.c_str(), [](double a) {return sqrt(a); });

	FILE *fp;
	fopen_s(&fp, param.densityFilename.c_str(), "w");
	for (int i=0;i<indices.size();i++)			
			fprintf(fp, "%d\n", indices[i].size());
	fclose(fp);	

	printf("total time: %f sec\n", taketime() - start_time);
	delete[] dataset.ptr(); 
	delete[] query.ptr();
	exit(1);
	//return 1;//very slow, why?
}
