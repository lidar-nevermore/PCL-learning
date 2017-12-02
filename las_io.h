#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include "laslib/lasreader.hpp"
#include "laslib/laswriter.hpp"
using namespace std;
using namespace pcl;

struct MyPoint
{
	PCL_ADD_POINT4D; 
	PCL_ADD_RGB;
	float userdata;
	int CountOfReturns;
	int ReturnNumber;
	int classification;
	friend ostream &operator<<(ostream &out,MyPoint pt);
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW	
}EIGEN_ALIGN16;
POINT_CLOUD_REGISTER_POINT_STRUCT(MyPoint,
	(float,x,x)
	(float,y,y)
	(float,z,z)
	(uint8_t,r,r)
	(uint8_t,g,g)
	(uint8_t,b,b)
	(uint8_t,a,a)
	(float,rgb,rgb)	
	(uint32_t,rgba,rgba)	
	(float,userdata,userdata)
	(int,CountOfReturns,CountOfReturns)
	(int,ReturnNumber,ReturnNumber)
	(int,classification,classification)	
	)
ostream &operator<<(ostream &out,MyPoint pt)
{
	out<<pt.x<<" "<<pt.y<<" "<<pt.z<<" "
		<<(int)pt.r<<" "<<(int)pt.g<<" "<<(int)pt.b<<" "
		<<pt.userdata<<" "<<pt.CountOfReturns<<" "
		<<pt.ReturnNumber<<" "<<pt.classification;
	return out;
};
/*
// pack r/g/b into rgb
 uint8_t r = 255, g = 0, b = 0;    // Example: Red color
uint32_t rgb = ((uint32_t)r << 16 | (uint32_t)g << 8 | (uint32_t)b);
 p.rgb = *reinterpret_cast<float*>(&rgb);
 \endcode

 To unpack the data into separate values, use:

 \code
 PointXYZRGB p;
 // unpack rgb into r/g/b
 uint32_t rgb = *reinterpret_cast<int*>(&p.rgb);
 uint8_t r = (rgb >> 16) & 0x0000ff;
 uint8_t g = (rgb >> 8)  & 0x0000ff;
 uint8_t b = (rgb)       & 0x0000ff;
*/

int	getCloudFromLas(const string &filename,pcl::PointCloud<MyPoint>::Ptr &m_cloud,bool IsScale )
{
LASreadOpener lasreadopener;
lasreadopener.set_file_name(filename.c_str());
LASreader* lasreader = lasreadopener.open();
if (lasreader==NULL) return -1;
int m_totalcount= lasreader->header.number_of_point_records;
double dx,dy;
if (IsScale)
{
	double minX=lasreader->header.min_x;
	double minY=lasreader->header.min_y;
	double maxX=lasreader->header.max_x;
	double maxY=lasreader->header.max_y;	
	dx=(minX+maxX)/2.0;
	dy=(minY+maxY)/2.0;
}
else
{
	dx=dy=0.0;
}
m_cloud->resize(m_totalcount);
int i=0;
while (lasreader->read_point()) 
{
	(*m_cloud)[i].x= lasreader->point.get_x()-dx;
	(*m_cloud)[i].y= lasreader->point.get_y()-dy;
	(*m_cloud)[i].z=lasreader->point.get_z();
	(*m_cloud)[i].userdata= lasreader->point.user_data;
	(*m_cloud)[i].r=lasreader->point.rgb[0]>>8;	
	(*m_cloud)[i].g=lasreader->point.rgb[1]>>8;	
	(*m_cloud)[i].b=lasreader->point.rgb[2]>>8;
	//(*m_cloud)[i].rgba=((uint32_t)(*m_cloud)[i].r << 16 | (uint32_t)(*m_cloud)[i].g << 8 | (uint32_t)(*m_cloud)[i].b);	
	
	(*m_cloud)[i].ReturnNumber= lasreader->point.return_number;
	(*m_cloud)[i].classification= lasreader->point.classification;
	(*m_cloud)[i].CountOfReturns= lasreader->point.number_of_returns_of_given_pulse;
	++i;
}
lasreader->close();
delete lasreader;
return 0;
}
void saveCloudToLas(const string &filename,pcl::PointCloud<MyPoint>::Ptr &m_cloud )
{
	LASheader lasheader;
	lasheader.x_scale_factor = 0.01;
	lasheader.y_scale_factor = 0.01;
	lasheader.z_scale_factor = 0.01;
	lasheader.x_offset = 0.0;
	lasheader.y_offset = 0.0;
	lasheader.z_offset = 0.0;
	lasheader.point_data_format = 2;
	lasheader.point_data_record_length = 26;	

	LASpoint laspoint;
	laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);

	// open laswriter
	LASwriteOpener laswriteopener;
	laswriteopener.set_file_name(filename.c_str());	
	LASwriter* laswriter = laswriteopener.open(&lasheader);		

	// write points
	size_t sz=m_cloud->size();
	for (size_t i = 0; i < sz; i++)
	{				
		laspoint.set_x((*m_cloud)[i].x);
		laspoint.set_y((*m_cloud)[i].y);
		laspoint.set_z((*m_cloud)[i].z);		
		laspoint.user_data=(*m_cloud)[i].userdata;
		laspoint.rgb[0]=static_cast<unsigned short>((*m_cloud)[i].r)<<8;
		laspoint.rgb[1]=static_cast<unsigned short>((*m_cloud)[i].g)<<8;
		laspoint.rgb[2]=static_cast<unsigned short>((*m_cloud)[i].b)<<8;
		laspoint.classification=(*m_cloud)[i].classification;
		laspoint.return_number=(*m_cloud)[i].ReturnNumber;		
		laspoint.number_of_returns_of_given_pulse=(*m_cloud)[i].CountOfReturns;	

		laswriter->write_point(&laspoint);
		laswriter->update_inventory(&laspoint);
	}
	// update the header
	laswriter->update_header(&lasheader, TRUE);
	// close the writer
	laswriter->close();;
}
