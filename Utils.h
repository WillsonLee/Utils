#ifndef LYXUTILS_H
#define LYXUTILS_H
//说明：一些常用工具函数
//项目-属性-配置属性-C/C++ -预处理器-预处理器定义中要添加_CRT_SECURE_NO_WARNINGS,否则lyxutils::debug::getTimeString函数会出现线程不安全报错
#include <string>
#include <iostream>
#include <time.h>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <thread>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

/**\ self made lyxutils package
  *\encapsulated some convenient functionalities
  *\some namespaces are defined within this utils namespace, for brevity, '#define' directive can be applied to give them shorter alias
  *\eg: #define lal lyxutils::algorithms, then lyxutils::algorithms namespace has the alias lal
  *\or better, you can declare a namespace and assign the lyxutils::some_inner_namespace to it
  *\eg: namespace myIO=lyxutils::io;
  */
namespace lyxutils
{
	namespace debug
	{
		/**\brief get the current time string
		  *\return the current time string in HH:MM:SS format
		  */
		std::string getTimeString();
		/**\brief create a folder named folderPath if it is not existed
		  *\param[in] folderPath the folder path to be created
		  *\return true if creation succeeded or the folder is already existed;false otherwise
		  */
		bool createFolder(const std::string folderPath);
		/**\brief open a log file for writing log
		  *\param[in] fileName the file name of log file
		  *\return writable log file(write in append mode)
		  */
		std::ofstream* openLogFile(std::string fileName);
		/**\brief write message to opened log file
		  *\param[in] log the file stream to write to
		  *\param[in] message the message to write
		  *\param[in] end the message end(default end is "\n",automatically carriage return)
		  */
		void writeLog(std::ofstream &log, std::string message, std::string end = "\n");
		/**\brief time the function f
		  *\param[in] f the function to be timed
		  *\note: the unit of time is millisecond,the function being timed should take no argument
		  *\input parameter or should be adapted using lambda function or bind function
		  *\function can be easily adapted using lambda function, for example,function like:
		  *\int f(params) can be adapted like time=timeIt([&result,params](){result=f(params);});
		  */
		template<class Func>
		double timeIt(Func f)
		{
			double start, stop, duration;
			start = clock();
			f();
			stop = clock();
			duration = stop - start;
			return duration;
		}
		/**\brief print console progress bar
		  *\param[in] value the current progress value(0~100)
		  *\param[in] delimiterL the left delimiter
		  *\param[in] increment the progress symbol
		  *\param[in] delimiterR the right deimiter
		  *\param[in] status a string to declare current status(eg. "processing","reading",etc)
		  *\note: when value is 0 or 100, it will automatically go to next line; 
		  *\progress update within 0-100 will update in the same line
		  *\note the update frequency should not be very 
		  *\for(int i=0;i<total;++i){
		  *\	if(i%(total/100)==0)lyxutils::consoleProgressBar(i*100/total);
		  *\	//some enduring repetitive procedure
		  *\}
		  */
		void consoleProgressBar(int value, std::string delimiterL = "[", std::string increment = ">", std::string delimiterR = "]", std::string status = "");
	}
	/**\brief encapsule some convenient string operation
	  */
	namespace str_utils
	{
		/**\brief split the string with given delimiter
		*\param[in] str the string to be splitted
		*\param[in] delimiter the delimiter(note:can only be single char, if several are given, the delimiter is any one of them)
		*/
		std::vector<std::string> split(const std::string &str, const std::string &delimiter);
		/**\brief split the string with a string
		  */
		std::vector<std::string> split_with_string(const std::string &str, const std::string &delString);
		/**\brief to lower case
		  */
		std::string lower(const std::string &str);
		/**\brief to upper case
		  */
		std::string upper(const std::string &str);
		/**\brief replace old string with new string
		  */
		std::string replace(const std::string &str, const std::string &oldStr, const std::string &newStr);
	}
	/**\brief encapsule some IO operation
	  */
	namespace io
	{
		  /**\brief read csv file
		    *\param[in] fileName the input file name
		    *\param[out] table the 2d vector to store the data
		    *\param[in] sep the separator to separate each column, usually ","
		    *\param[out] report the reading process report(record open file name,time,bad lines, etc)
		    *\param[in] fields the number of fields to read(-1 means use total number of fields in the first line of the file)
		    *\param[in] read_header whether to read first line as header or not
		    *\param[out] headers the headers of the table(if read_header is true)
		    *\param[in] numberOfRedirect number of fields to be redirected(not stored in table),should be less than fields
		    *\param[in] directFunc the function to redirect first numberOfRedirect of fields
		    *\note:parameter of directFunc must have the same dimension as numberOfRedirect
		    *\the last parameter directFuncs can be defined as follows if trying to redirect to point cloud(cloud_in is a global pc pointer)
					void direct(const std::vector<float> &vec) {
						pcl::PointXYZRGB pt = pcl::PointXYZRGB();
						pt.x = vec[0];
						pt.y = vec[1];
						pt.z = vec[2];
						pt.r = (uint8_t)(vec[3]);
						pt.g = (uint8_t)(vec[4]);
						pt.b = (uint8_t)(vec[5]);
						cloud_in->points.push_back(pt);
					}
			*\addition: if the columns of data is less than required, it will fill with 0s; 
			*\			if data column is larger than required, the rest will be abbreviated
		    */
		bool read_csv(const std::string &fileName, std::vector<std::vector<float> > &table, const std::string &sep,
		        std::string &report, std::vector<std::string> &headers, int fields = -1, bool read_header = false,
		        int numberOfRedirect = 0, void(*directFunc)(const std::vector<float>&) = NULL);
		/**\brief write csv file
		  */
		template<class DataType>
		void write_csv(const std::string & fileName, const std::vector<std::vector<DataType>>& table, const std::string & sep)
		{
			std::ofstream ofs(fileName);
			ofs << std::setprecision(8);
			for (int i = 0; i < table[0].size(); ++i) {
				ofs << (*table.begin())[i];
				std::for_each(table.begin() + 1, table.end(), [i, &ofs, sep](const std::vector<DataType> &vec) {ofs << sep << vec[i]; });
				ofs << std::endl;
			}
			if (ofs.is_open()) {
				ofs.close();
			}
		}
		/**\brief write points to csv file
		  */
		template<class PointXY,class allocator=std::allocator<PointXY> >
		void points_to_csv(std::string fileName, std::vector<PointXY, allocator> &points, std::string sep = ",", int precision = 8)
		{
			std::ofstream ofs(fileName);
			ofs << std::setprecision(precision);
			for (int i = 0; i < points.size(); ++i) {
				ofs << points[i].x << sep << points[i].y << std::endl;
			}
			if (ofs.is_open()) {
				ofs.close();
			}
		}
		/**\brief write 3D points to file
		  */
		template<class PointXYZ, class allocator = std::allocator<PointXYZ> >
		void points3D_to_csv(std::string fileName, std::vector<PointXYZ, allocator> &points, std::string sep = ",", int precision = 8)
		{
			std::ofstream ofs(fileName);
			ofs << std::setprecision(precision);
			for (int i = 0; i < points.size(); ++i) {
				ofs << points[i].x << sep << points[i].y << sep << points[i].z << std::endl;
			}
			if (ofs.is_open()) {
				ofs.close();
			}
		}
		/**\brief write 3D colored points to file
		*/
		template<class PointXYZRGB, class allocator = std::allocator<PointXYZRGB> >
		void points3DRGB_to_csv(std::string fileName, std::vector<PointXYZRGB, allocator> &points, std::string sep = ",", int precision = 8)
		{
			std::ofstream ofs(fileName);
			ofs << std::setprecision(precision);
			for (int i = 0; i < points.size(); ++i) {
				ofs << points[i].x << sep << points[i].y << sep << points[i].z << sep
					<< points[i].r << sep << points[i].g << sep << points[i].b << std::endl;
			}
			if (ofs.is_open()) {
				ofs.close();
			}
		}
		/**\brief check file for existence
		  */
		bool fileExists(std::string fileName);
	}
	/**\brief encapsulate some alorithms and Data structures
	  */
	namespace algorithms 
	{
		const double PI = 3.1415926535897932;

		/**\brief compute the average x y z coordinates of points contained in pList
		  *\param[in] pList the vector contains the points
		  *\param[out] average_x,average_y,average_z the computed averages
		  */
		template<class PointXYZ,class allocator=std::allocator<PointXYZ> >
		void computeAverages(const std::vector<PointXYZ,allocator> &pList, double &average_x, double &average_y, double &average_z)
		{
			average_x = 0;
			average_y = 0;
			average_z = 0;
			for (int i = 0; i < pList.size(); ++i) {
				average_x += pList[i].x;
				average_y += pList[i].y;
				average_z += pList[i].z;
			}
			if (pList.size() > 0) {
				average_x = average_x / pList.size();
				average_y = average_y / pList.size();
				average_z = average_z / pList.size();
			}
		}
		/**\brief use least square method to fit plane
		  *\param[in] pList the input point list
		  *\param[out] a,b,c,d the output plane parameters, where a^2+b^2+c^2=1
		  */
		template<class PointXYZ,class allocator=std::allocator<PointXYZ> >
		void leastSquarePlaneFit(const std::vector<PointXYZ,allocator> &pList, double &a, double &b, double &c, double &d)
		{
			Eigen::Matrix3d A;
			double average_x, average_y, average_z;
			double sum_xx = 0, sum_xy = 0, sum_xz = 0;
			double sum_yy = 0, sum_yz = 0;
			double sum_zz = 0;
			computeAverages(pList, average_x, average_y, average_z);
			for (int i = 0; i < pList.size(); ++i) {
				sum_xx += (pList[i].x - average_x)*(pList[i].x - average_x);
				sum_xy += (pList[i].x - average_x)*(pList[i].y - average_y);
				sum_xz += (pList[i].x - average_x)*(pList[i].z - average_z);
				sum_yy += (pList[i].y - average_y)*(pList[i].y - average_y);
				sum_yz += (pList[i].y - average_y)*(pList[i].z - average_z);
				sum_zz += (pList[i].z - average_z)*(pList[i].z - average_z);
			}
			A << sum_xx, sum_xy, sum_xz,
				sum_xy, sum_yy, sum_yz,
				sum_xz, sum_yz, sum_zz;
			Eigen::EigenSolver<Eigen::Matrix3d> es(A);
			Eigen::Matrix3d D = es.pseudoEigenvalueMatrix();
			Eigen::Matrix3d V = es.pseudoEigenvectors();
			int mark = 0;
			double lamda = D(0, 0);
			for (int i = 0; i < 3; ++i) {
				if (D(i, i) < lamda) {
					mark = i;
					lamda = D(i, i);
				}
			}
			a = V(0, mark);
			b = V(1, mark);
			c = V(2, mark);
			d = -a*average_x - b*average_y - c*average_z;
			if (c < 0) {//keep c always positive
				a = -a;
				b = -b;
				c = -c;
				d = -d;
			}
		}
		/**\brief find the distance of point p to plane defined by ax+by+cz+d=0
		  *\where a^2+b^2+c^2=1
		  */
		template<class PointXYZ>
		double distanceToPlane(const PointXYZ &p, double a, double b, double c, double d)
		{
			return abs(a*p.x + b*p.y + c*p.z + d);
		}
		/**\brief find the distance of point (x,y,z) to plane defined by ax+by+cz+d=0
		  *\where a^2+b^2+c^2=1
		  */
		double distanceToPlane(double x, double y, double z, double a, double b, double c, double d);
		/**\brief find the average distance of given points to a normalized plane with given parameters
		  *\param[in] pList the point vector
		  *\param[in] a,b,c,d the plane coefficient
		  *\param[out] average_dis the average distance of points to plane
		  */
		template<class PointXYZ, class allocator = std::allocator<PointXYZ> >
		void findAverageDistanceToPlane(const std::vector<PointXYZ, allocator> &pList, const double &a,
			const double &b, const double &c, const double &d, double &average_dis)
		{
			average_dis = 0;
			for (int i = 0; i < pList.size(); ++i) {
				average_dis += distanceToPlane(pList[i].x, pList[i].y, pList[i].z, a, b, c, d);
			}
			average_dis = average_dis / pList.size();
		}
		/**\brief region growing algorithm,give a seed unit and grow out a whole region.mark each unit being growned
		  *\param[in] seed the given seed point(type is unit,define whichever type you think can represent the growing unit)
		  *\param[in] findNeighbors a function or functor that can find the neighbors of a given unit
		  *\param[in] growingCondition the growing condition for this algorithm
		  *\param[out] growingResult the growing result,a vector that contains all the growing resultant units
		  *\function or functor declarations:
		  *\std::vector<Unit>* findNeighbors(const Unit &center)
		  *\bool growingCondition(const Unit &center,const Unit &neighbor)
		  *\
		  *\note: this tenmplate function is ill-designed, UnaryOperation and BinaryOperation should respectively 
		  *\		be functions as:`std::vector<Unit>* ()(Unit &center)` and `bool ()(Unit &center,Unit &target)`
		  */
		template<class Unit, class UnaryOperation, class BinaryPredict>
		void regionGrowing(Unit &seed, UnaryOperation findNeighbors,
			BinaryPredict growingCondition, std::vector<Unit> &growingResult)
		{
			std::vector<Unit> *recursion = new std::vector<Unit>();
			growingResult.push_back(seed);
			recursion->push_back(seed);
			while (recursion->size() > 0)
			{
				Unit indexesTemp = (*recursion)[recursion->size() - 1];
				recursion->pop_back();
				std::vector<Unit>* neighbors = findNeighbors(indexesTemp);
				for (auto iterator = neighbors->begin(); iterator != neighbors->end(); iterator++) {
					if (growingCondition(indexesTemp, *iterator)) {
						growingResult.push_back(*iterator);
						recursion->push_back(*iterator);
					}
				}
				delete neighbors;
			}
			delete recursion;
		}
        /**\brief find the orientation of the ordered triplet (p1,p2,p3)
          *\function returns following values:
          *\0-->if p1,p2,p3 are colinear (or if at least two of them are indentical)
          *\1-->if clockwise
          *\2-->if counterclockwise
          */
        template<class PointXY>
        int orientation(const PointXY &p1, const PointXY &p2, const PointXY &p3)
        {
            double val = (p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y);
            if (val == 0) return 0;  // colinear
            return (val > 0) ? 1 : 2; // clock or counterclock wise
        }
        /**\brief compare the two polar angles x-p1-p2 and x-p1-p3(counterclockwise)
          *\return true if x-p1-p2 is smaller than x-p1-p3(or if p1p2<p1p3 when p1,p2,p3 are colinear)
          */
        template<class PointXY>
        bool compareAngle(const PointXY &p1, const PointXY &p2, const PointXY &p3)
        {
            int o = orientation(p1, p2, p3);
            if (o == 0) {
                //return distance2DSquared(p1, p2) < distance2DSquared(p1, p3);
                //return distance2D(p1, p2) < distance2D(p1, p3);
                //assume p1 is the left most one
                double v12_x = p2.x - p1.x;
                double v12_y = p2.y - p1.y;
                double v23_x = p3.x - p2.x;
                double v23_y = p3.y - p2.y;
                double innerProduct = v12_x*v23_x + v12_y*v23_y;
                return innerProduct > 0;
            }
            return (o == 2) ? true : false;
        }
		/**\brief use Graham Scan to get the convex hull(anti-clockwise starting from min x) border of a given set of points
		  *\param[in] pList the given set of points
		  *\param[out] border the output set of border points
		  */
		template<class PointXY, class allocator=std::allocator<PointXY> >
		void findConvexHullBorder(const std::vector<PointXY,allocator> &pList, std::vector<PointXY,allocator> &border)
		{
			if (pList.size() < 3)
				throw std::runtime_error("number of points too small(<3) to define a convex hull");
			//copy to border,operate on border
			double x_min = pList[0].x;
			int x_min_index = 0;
			border.push_back(pList[0]);
			for (int i = 1; i < pList.size(); ++i) {
				border.push_back(pList[i]);
				if (pList[i].x < x_min || pList[i].x == x_min&&pList[i].y < pList[x_min_index].y) {
					x_min = pList[i].x;
					x_min_index = i;
				}
			}
			//move smallest x to head
			PointXY pointTemp = border[x_min_index];
			border[x_min_index] = border[0];
			//order anti-clock-wise
			//sort(border.begin() + 1, border.end(), [&pointTemp](PointType *p1, PointType *p2)
			//->bool {return compareAngle(*pointTemp, *p1, *p2); });
			border.erase(border.begin());
			merge_sort(border, [&pointTemp](PointXY p1, PointXY p2)
				->bool {return compareAngle(pointTemp, p1, p2); });
			border.insert(border.begin(), pointTemp);
			//remove collinear point, maintain only the furthest
			for (auto it = border.begin() + 1; it != border.end();) {
				if ((it + 1) != border.end() && orientation(pointTemp, *it, *(it + 1)) == 0) {
					it = border.erase(it);
				}
				else {
					it++;
				}
			}
			if (border.size() < 3)
				throw std::runtime_error("number of border points too small(<3)");
			if (border.size() == 3)
			{
				return;
			}
			//remove concave vertex
			auto first = border.begin() + 1;
			auto second = first + 1;
			auto third = second + 1;
			while (third != border.end()) {
				int orien = orientation(*first, *second, *third);
				if (orien != 2) {
					border.erase(second);
					second = first;
					first--;
					third = second + 1;
				}
				else {
					first++;
					second++;
					third++;
				}
			}
		}
		/**\brief find convex hull using javis march algorithm, the border point is in anti-clockwise order starting from mkinimum y value point
		  *\param[in] pList input point vector
		  *\param[out] border output vertices
		  *\note: the point class should have default constructor that requires no argument,should have clone & copy constructor
		  *\	  the temporal complexity is O(nH), where H is the number of vertices;
		  *\	  should be faster than graham-O(nlogn) when number of vertices are limited
		  */
		template<class PointXY, class allocator = std::allocator<PointXY> >
		void javisMarch(const std::vector<PointXY, allocator> &pList, std::vector<PointXY, allocator> &border)
		{
			if (pList.size() < 3)
				throw std::runtime_error("number of points too small(<3) to define a convex hull");
			//get the vector comprised by two points
			auto vector2P = [](const PointXY &p1, const PointXY &p2)->PointXY {
				PointXY result;
				result.x = p2.x - p1.x;
				result.y = p2.y - p1.y;
				return result;
			};
			auto turningAngle = [&vector2P](const PointXY *p1, const PointXY *p2, const PointXY *p3)->float
			{
				PointXY v1;
				PointXY v2;
				if (p1 == NULL) {
					v1.x = 1;
					v1.y = 0;
				}
				else {
					v1 = vector2P(*p1, *p2);
				}
				v2 = vector2P(*p2, *p3);
				return (v1.x*v2.x + v1.y*v2.y) / (sqrt(v1.x*v1.x + v1.y*v1.y)*sqrt(v2.x*v2.x + v2.y*v2.y));
			};
			double yMin = pList[0].y;
			int min_index = 0;
			std::vector<const PointXY*> remainPts;
			for (int i = 0; i < pList.size(); ++i) {
				remainPts.push_back(&pList[i]);
				if (pList[i].y < yMin) {
					yMin = pList[i].y;
					min_index = i;
				}
			}
			const PointXY *pTemp = remainPts[0];
			remainPts[0] = remainPts[min_index];
			remainPts[min_index] = pTemp;
			const PointXY *p1, *p2;
			p1 = NULL;
			p2 = remainPts[0];
			border.push_back(*p2);
			const PointXY *starter = p2;
			while (true) {
				if (remainPts.size() < 2)break;
				float maxCos = -2.0;
				auto p3_pntr = remainPts.begin();
				for (auto it = remainPts.begin(); it != remainPts.end(); ) {
					if ((*it)->x == p2->x && (*it)->y == p2->y) {
						if (it != remainPts.begin()) {
							it = remainPts.erase(it);
						}
						else {
							++it;
						}
					}
					else {
						float cosVal = turningAngle(p1, p2, *it);
						if (cosVal > maxCos) {
							maxCos = cosVal;
							p3_pntr = it;
						}
						++it;
					}
				}
				p1 = p2;
				p2 = *p3_pntr;
				border.push_back(*p2);
				remainPts.erase(p3_pntr);
				if (p2==starter) {
					border.pop_back();
					break;
				}
			}
		}
		/**\brief creationg mask image that marks the interior points inside a convex hull border
		  *\param[in] border the convex hull border (the point coordinate x y should be integer
		  *\param[in] rows,cols the rows and cols of the mask image
		  *\param[out] access the function that sets image value at x,y position(like access(x,y)=data),
		  *\	eg. auto accessMat=[&mask](int x,int y)->uchar&{return mask.at<uchar>(x,y);};
		  *\		//above lambda function can used to access a cv::Mat of type CV_8UC1
		  *\		convexhullIntereriorMask(border,mask.rows,mask.cols,accessMat,255);
		  *\param[in] value the foreground value(the value of points inside convex hull)
		  */
		template<class PointXY,class allocator=std::allocator<PointXY>,class DataType,class Accessor>
		void convexhullIntereriorMask(const std::vector<PointXY, allocator> &border, int rows, int cols, Accessor access, DataType value)
		{
			std::vector<int> left(rows, cols), right(rows, -1);//record left border y value and right border y value
			int xMinIdx = 0;
			int xMin = border[0].x;
			for (int i = 1; i < border.size();++i) {
				if (border[i].x < xMin) {
					xMin = border[i].x;
					xMinIdx = i;
				}
			}
			for (int i = xMinIdx; i < xMinIdx + border.size(); ++i) {//start from min x to max x then back to min x
				int x1 = border[i%border.size()].x;
				int y1 = border[i%border.size()].y;
				int x2 = border[(i + 1) % border.size()].x;
				int y2 = border[(i + 1) % border.size()].y;
				if (x1 < x2) {//left side
					for (int k = x1; k < x2; ++k) {
						left[k] = float(y1*(x2 - k) + y2*(k - x1)) / (x2 - x1);
					}
				}
				else if (x1 == x2) {
					left[x1] = y1;
					right[x2] = y2;
				}
				else {
					for (int k = x1; k > x2; --k) {
						right[k] = float(y1*(k - x2) + y2*(x1 - k)) / (x1 - x2);
					}
				}
			}
			for (int i = 0; i < rows; ++i) {//fill the mask image
				for (int j = left[i]; j < right[i]; ++j) {
					access(i, j) = value;
				}
			}
		}
		/**\brief calculate the area of a convex hull
		  *\param[in] border the border p[oint list��can be derived from lyxutils::algorithm::findConvexHull
		  *\note:the point list is in the order of clockwise or anti-clockwise
		  */
		template<class PointXY,class allocator=std::allocator<PointXY> >
		double getConvexHullArea(const std::vector<PointXY,allocator> &border)
		{
			if (border.size() < 3)
				throw std::runtime_error("number of points too small(<3) to enclose");
			auto assign = [](const PointXY &from, PointXY &to) {to.x = from.x; to.y = from.y; };//�����û�п�¡���캯��
			PointXY p1, p2, p3;
			double area = 0;
			int i = 1;
			assign(border[0], p1);
			while (i + 1 < border.size()) {
				assign(border[i], p2);
				assign(border[i + 1], p3);
				area += triangleArea(p1, p2, p3);
				++i;
			}
			return area;
		}
		/**\brief calculate area of a triangle(2D)
		*/
		template<class PointXY>
		double triangleArea(const PointXY &vertice1, const PointXY &vertice2, const PointXY &vertice3)
		{
			double m_x = vertice2.x - vertice1.x;
			double m_y = vertice2.y - vertice1.y;
			double n_x = vertice3.x - vertice1.x;
			double n_y = vertice3.y - vertice1.y;
			return abs(m_x*n_y - m_y*n_x) / 2;
		}
		/**\brief find the squared Euclidean distance between two points
		*/
		template<class PointXY>
		double distance2DSquared(const PointXY &p1, const PointXY &p2)
		{
			return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
		}
		/**\brief calculate distance between two 2D planes
		  */
		template<class PointXY>
		double distance2D(const PointXY &p1, const PointXY &p2)
		{
			return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
		}
		/**\brief insert sort
		  */
		template<class Iterator,class BinaryComp>
		void insert_sort(Iterator begin, Iterator end, BinaryComp pred) 
		{
			for (auto it = begin; it + 1 != end; ++it) {
				for (auto it2 = it + 1; it2 != begin; --it2) {
					if (pred(*it2, *(it2 - 1))) {
						auto temp = *(it2 - 1);
						*(it2 - 1) = *it2;
						*it2 = temp;
					}
					else {
						break;
					}
				}
			}
		}
		/**\brief merge sort, with nlogn temporal complexity
		*/
		template<class DataType,class BinaryComp>
		void merge_sort(std::vector<DataType> &list, BinaryComp pred)
		{
			if (list.size() <= 1)
				return;
			int half = list.size() / 2;
			std::vector<DataType> pre;
			std::vector<DataType> post;
			for (int i = 0; i < half; ++i) {
				pre.push_back(list[i]);
			}
			for (int i = half; i < list.size(); ++i) {
				post.push_back(list[i]);
			}
			merge_sort(pre, pred);
			merge_sort(post, pred);
			//�鲢
			auto it1 = pre.begin();
			auto it2 = post.begin();
			for (int i = 0; i < list.size(); ++i) {
				if (it1 == pre.end()) {
					list[i] = *it2++;
				}
				else if (it2 == post.end()) {
					list[i] = *it1++;
				}
				else {
					list[i] = pred(*it1, *it2) ? *it1++ : *it2++;
				}
			}
		}
	}
	/**\brief encapsulate color related functions
	*/
	namespace colors 
	{
		/**\brief convert rgb to hsv
		  *\param[in] rgb:0~255
		  *\param[out] hsv:h[0,359);s[0,1];v:0~255
		  */
		void rgb2hsv(const std::vector<int> &rgb, std::vector<float> &hsv);
		/**\brief convert hsv to rgb
		*\param[in] hsv:h[0,359);s[0,1];v:0~255
		*\param[out] rgb:0~255
		*/
		void hsv2rgb(const std::vector<float> &hsv, std::vector<int> &rgb);
	}
	/**\brief encapsulate numpy equivalent functionalities
	 * \this is poorly designed. don't use it unless you have nothing to use
	  */
	namespace numpy
	{
		/**\brief limited version of broadcast like that in numpy-broadcast one number to a vector
		  *\param[in] val the number
		  *\param[in] num number of elements in the output vector
		  *\param[out] vec the resulting vector
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void broadcast(const NumType &val, int num, std::vector<NumType, allocator> &vec)
		{
			vec.clear();
			vec.resize(num);
			for (auto &v : vec) { v = val; }
		}
		/**\brief limited version of broadcast like that in numpy-broadcast one vector to a 2d vector(matrix)
		  *\param[in] vec the input vector
		  *\param[in] repeats number of repeats in the output 2d vector
		  *\param[out] mat the resulting matrix
		  *\param[in] prepend true if broadcast along first dimension(namely repeat rows),false if broadcast alonmg second dimesion(namely repeat columns)
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void broadcast(const std::vector<NumType, allocator> &vec, int repeats, std::vector<std::vector<NumType, allocator> > &mat, bool prepend = true)
		{
			mat.clear();
			mat.resize(prepend ? repeats : vec.size());
			if (prepend) {
				for (int i = 0; i < repeats; ++i) {
					for (auto it = vec.begin(); it != vec.end(); ++it) {
						mat[i].push_back(*it);
					}
				}
			}
			else {
				for (int i = 0; i < mat.size(); ++i) {
					for (int j = 0; j < repeats; ++j) {
						mat[i].push_back(vec[i]);
					}
				}
			}
		}
		/**\brief limited version of broadcast like that in numpy-broadcast one number to a matrix
		  *\param[in] val the number
		  *\param[in] row,col the number of rows and cols in the resulting matrix
		  *\param[out] mat the resulting matrix
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void broadcast(const NumType &val, int row, int col, std::vector<std::vector<NumType, allocator> > &mat)
		{
			mat.clear();
			mat.resize(row);
			for (int i = 0; i < row; ++i) {
				for (int j = 0; j < col; ++j) {
					mat[i].push_back(val);
				}
			}
		}
		/**\brief transpose the matrix
		  *\param[in] mat the input matrix
		  *\param[out] out_mat the transposed matrix
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void transpose(const std::vector<std::vector<NumType, allocator> > &mat, std::vector<std::vector<NumType, allocator> > &out_mat)
		{
			out_mat.clear();
			if (mat.size() == 0) {
				return;
			}
			out_mat.resize(mat[0].size());
			for (int i = 0; i < mat[0].size(); ++i) {
				for (int j = 0; j < mat.size(); ++j) {
					out_mat[i].push_back(mat[j][i]);
				}
			}
		}
		template<class NumType,class allocator=std::allocator<NumType> >
		void range(std::vector<NumType, allocator> &result,NumType start, NumType end, NumType step = 1)
		{
			if (end == start || step == 0) {
				return;
			}
			if ((end - start)*step < 0) {//step is opposite to start-to-end direction
				step = -step;
			}
			NumType r2 = (end - start)*(end - start);
			result.clear();
			for (NumType i = start; (i - start)*(i - start) < r2; i += step) {
				result.push_back(i);
			}
		}
		template<class NumType, class allocator = std::allocator<NumType> >
		void range(std::vector<NumType, allocator> &result,NumType end)
		{
			range<NumType,allocator>(result, 0, end);
		}
		template<class NumType, class allocator = std::allocator<NumType> >
		void linspace(std::vector<NumType, allocator> &result,NumType start, NumType end, int num)
		{
			if (num < 1) {
				throw std::runtime_error("illegal parameter num,num should be greater than 1(or equal to 1,only when start==end)");
			}
			if (num == 1) {
				if (start == end) {
					result.push_back(start);
					return;
				}
				else {
					throw std::runtime_error("illegal parameter num,num should be greater than 1(or equal to 1,only when start==end)");
				}
			}
			NumType step = (end - start) / (num - 1);
			range(result,start, end, step);
			result.push_back(end);
		}
		/**\brief reshape the one dimensional vector into two dimensional
		  *\note: row times col should be equal to size of vec;if not the col will be automatically adjusted to satisfy the demand when row!=0 and size of vec is divisible by row
		  *\param[in] vec the 1d vector
		  *\param[in] row,col the output vector number of rows and cols
		  *\param[out] out_vec the resulting 2d vector
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void reshape(const std::vector<NumType, allocator> &vec, int row, int col, std::vector<std::vector<NumType, allocator> > &out_vec)
		{
			if (row*col != vec.size()) {
				if (row != 0 && vec.size() % row == 0) {
					col = vec.size() / row;
					std::cout << "row times col is not equal to size of vec, the col is changed to " << col << ", the row remains " << row << std::endl;
				}
				else {
					throw std::runtime_error("row times col is not equal to size of vec!");
				}
			}
			out_vec.clear();
			out_vec.resize(row);
			for (int i = 0; i < row; ++i) {
				for (int j = 0; j < col; ++j) {
					out_vec[i].push_back(vec[i*col + j]);
				}
			}
		}
		/**\brief expand the two dimensional vector into 1d one
		  *\param[in] vec the input two d vector
		  *\param[out] out_vec the expanded vector
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void flat(const std::vector<std::vector<NumType, allocator> > &vec, std::vector<NumType, allocator> &out_vec)
		{
			if (vec.size() == 0 || vec[0].size() == 0) {
				return;
			}
			out_vec.clear();
			for (int i = 0; i < vec.size(); ++i) {
				for (int j = 0; j < vec[i].size(); ++j) {
					out_vec.push_back(vec[i][j]);
				}
			}
		}
		template<class NumType, class allocator = std::allocator<NumType> >
		void print(const std::vector<NumType, allocator> &vec, bool column = false, bool keepSameLine = false, bool verbose = false)
		{
			std::cout << "[";
			if (vec.size() > 0)
			{
				std::cout << vec[0];
				if (!column) {
					for (auto it = vec.begin() + 1; it != vec.end(); ++it) {
						std::cout << ",\t" << *it;
					}
				}
				else {
					for (auto it = vec.begin() + 1; it != vec.end(); ++it) {
						std::cout << "," << std::endl << " " << *it;
					}
				}
			}
			std::cout << "]";
			if (verbose) {
				std::cout << std::endl << "number of elements:" << vec.size() << ",type name of element:" << std::string(typeid(NumType).name());
			}
			if (!keepSameLine) {
				std::cout << std::endl;
			}
		}
		template<class NumType, class allocator = std::allocator<NumType> >
		void print(const std::vector<std::vector<NumType, allocator> > &mat, bool keepSameLine = false, bool verbose = false)
		{
			std::cout << "[";
			for (auto it = mat.begin(); it != mat.end(); ++it) {
				if (it != mat.begin()) {
					std::cout << " ";
				}
				print(*it, false, it == (mat.end() - 1));
			}
			std::cout << "]";
			if (verbose) {
				std::cout << std::endl << "shape:" << mat.size() << "x" << mat[0].size() << ",type name of element:" << std::string(typeid(NumType).name());
			}
			if (!keepSameLine) {
				std::cout << std::endl;
			}
		}
		/**\brief get the summation
		  */
		template<class NumType,class allocator=std::allocator<NumType> >
		NumType sum(const std::vector<NumType,allocator> &vec)
		{
			NumType result = 0;
			for (auto it = vec.begin(); it != vec.end(); ++it)result += (*it);
			return result;
		}
		/**\brief get the summation along one axis
		  *\param[in] vec input 2D qrray(first dimension is row)
		  *\param[out] out_vec output 1D array-the summation
		  *\param[in] axis the direction of summation��default -1��namely add all elements up��
		  *\note:axis=0->sum along row direction; axis=1->sum along column direction; axis=-1->sum all the elements
		  *\example:
		  *\std::vector<float> sum;
		  *\lyxutils::stats::sum(input_vec,sum,axis=1);//sum all the rows
		*/
		template<class NumType, class allocator = std::allocator<NumType> >
		void sum(const std::vector<std::vector<NumType, allocator> > &mat, std::vector<NumType, allocator> &out_vec, int axis = -1)
		{
			if (mat.size() == 0 || mat[0].size() == 0) {
				out_vec.clear();
				out_vec.push_back(0);
				return;
			}
			if (axis == 0 || axis == -1)
			{
				out_vec.resize(mat[0].size());
				for (NumType &s : out_vec)s = 0;
				for (int i = 0; i < out_vec.size(); ++i)
				{
					for (auto &row : mat) {
						out_vec[i] += row[i];
					}
				}
				if (axis == -1) {
					NumType result = sum(out_vec);
					out_vec.clear();
					out_vec.push_back(result);
				}
			}
			else if (axis == 1)
			{
				out_vec.resize(mat.size());
				for (NumType &s : out_vec)s = 0;
				for (int i = 0; i < out_vec.size(); ++i)
				{
					for (auto &s : mat[i]) {
						out_vec[i] += s;
					}
				}
			}
			else {
				throw std::invalid_argument("invalid axis parameter!");
			}
		}
		/**\brief get average of 1D vector
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		NumType mean(const std::vector<NumType, allocator> &vec)
		{
			return sum(vec) / vec.size();
		}
		/**\brief get averages of 2D array, refer to lyxutils::stats::sum
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void mean(const std::vector<std::vector<NumType, allocator> > &mat, std::vector<NumType, allocator> &out_vec, int axis = -1)
		{
			sum(mat, out_vec, axis);
			if (axis == -1) {
				out_vec[0] /= (mat.size()*mat[0].size());
			}
			else if (axis == 0) {
				for (int i = 0; i < out_vec.size(); ++i) {
					out_vec[i] /= mat.size();
				}
			}
			else if (axis == 1) {
				for (int i = 0; i < out_vec.size(); ++i) {
					out_vec[i] /= mat[i].size();
				}
			}
			else {
				throw std::invalid_argument("invalid axis parameter!");
			}
		}
		/**\brief get the standard deviation of 1D vector
		  *\param[in] vec input vector
		  *\param[in] ddof degree of freedom
		  *\std=sqrt(sum((xi-xbar)^2)/(n-ddof)),where xi represents each ith element, xbar is the average x value
		  *\n is number of elements and ddof is degree of freedom
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		NumType std(const std::vector<NumType, allocator> &vec, int ddof = 1)
		{
			NumType xbar = mean(vec);
			int n = vec.size();
			std::vector<NumType, allocator> xi2;
			apply(vec, xi2, [ddof, xbar, n](NumType x)->float {return (x - xbar)*(x - xbar) / (n - ddof); });
			return sqrt(sum(xi2));
		}
		/**\brief get the standard deviation of 2D matrix
		  *\param[in] mat input matrix
		  *\param[out] out_vec output
		  *\param[in] ddof degree of freedom
		  *\std=sqrt(sum((xi-xbar)^2)/(n-ddof)),where xi represents each ith element, xbar is the average x value
		  *\n is number of elements and ddof is degree of freedom
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void std(const std::vector<std::vector<NumType, allocator> > &mat, std::vector<NumType, allocator> &out_vec, int ddof = 1, int axis = -1)
		{
			if (mat.size() == 0 || mat[0].size() == 0) {
				out_vec.clear();
				out_vec.push_back(0);
				return;
			}
			if (axis == -1) {
				std::vector<NumType, allocator> xbar;
				mean(mat, xbar);
				NumType m = xbar[0];
				int n = mat.size()*mat[0].size();
				NumType sqrSum = 0;
				for (int i = 0; i < mat.size(); ++i) {
					for (int j = 0; j < mat[i].size(); ++j) {
						sqrSum += (mat[i][j] - m)*(mat[i][j] - m) / (n - ddof);
					}
				}
				out_vec.clear();
				out_vec.push_back(sqrt(sqrSum));
			}
			else if (axis == 0) {
				out_vec.resize(mat[0].size());
				std::vector<NumType, allocator> rowMean;
				mean(mat, rowMean, 0);
				int n = mat.size();
				for (auto &s : out_vec) { s = 0; }
				for (int i = 0; i < out_vec.size(); ++i) {
					for (int j = 0; j < mat.size();++j) {
						out_vec[i] += (mat[j][i] - rowMean[i])*(mat[j][i] - rowMean[i]) / (n - ddof);
					}
				}
				for (auto &s : out_vec) { s = sqrt(s); }
			}
			else if (axis == 1) {
				out_vec.resize(mat.size());
				std::vector<NumType, allocator> colMean;
				mean(mat, colMean, 1);
				int n = mat[0].size();
				for (auto &s : out_vec) { s = 0; }
				for (int i = 0; i < out_vec.size(); ++i) {
					for (int j = 0; j < mat[0].size(); ++j) {
						out_vec[i] += (mat[i][j] - colMean[i])*(mat[i][j] - colMean[i]) / (n - ddof);
					}
				}
				for (auto &s : out_vec) { s = sqrt(s); }
			}
			else {
				throw std::invalid_argument("invalid axis parameter!");
			}
		}
		/**\brief add two vectors
		  *\param[in] vec1,vec2 the input vector of the same size
		  *\param[out] out_vec sum of input vectors
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void add(const std::vector<NumType, allocator> &vec1, const std::vector<NumType, allocator> &vec2, std::vector<NumType, allocator> &out_vec)
		{
			out_vec.clear();
			if (vec1.size() == 0 && vec2.size() == 0) {
				//do nothing if sizes of both vectors are zero
			}
			else if (vec1.size() == vec2.size()) {
				for (auto it1 = vec1.begin(), it2 = vec2.begin(); it1 != vec1.end(); ++it1, ++it2) {
					out_vec.push_back((*it1)+(*it2));
				}
			}
			else {
				throw std::runtime_error("size of vec1 is not compatible with that of vec2!");
			}
		}
		/**\brief add two matrices
		  *\param[in] matrix1,matrix2 the input matrices of the same size
		  *\param[out] out_mat sum of input matrices
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void add(const std::vector<std::vector<NumType, allocator> > &matrix1, const std::vector<std::vector<NumType, allocator> > &matrix2, std::vector<std::vector<NumType, allocator> > &out_mat)
		{
			out_mat.clear();
			if (matrix1.size() == 0 && matrix2.size() == 0) {

			}
			else if (matrix1.size() == matrix2.size()) {
				for (auto it1 = matrix1.begin(), it2 = matrix2.begin(); it1 != matrix1.end(); ++it1, ++it2) {
					out_mat.push_back(std::vector<NumType, allocator>());
					add(*it1, *it2, *(out_mat.end() - 1));
				}
			}
			else {
				throw std::runtime_error("size of matrix1 is not compatible with that of matrix2!");
			}
		}
		/**\brief minus two vectors
		  *\param[in] vec1,vec2 the input vector of the same size
		  *\param[out] out_vec result of input vectors
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void minus(const std::vector<NumType, allocator> &vec1, const std::vector<NumType, allocator> &vec2, std::vector<NumType, allocator> &out_vec)
		{
			out_vec.clear();
			if (vec1.size() == 0 && vec2.size() == 0) {
				//do nothing if sizes of both vectors are zero
			}
			else if (vec1.size() == vec2.size()) {
				for (auto it1 = vec1.begin(), it2 = vec2.begin(); it1 != vec1.end(); ++it1, ++it2) {
					out_vec.push_back((*it1) - (*it2));
				}
			}
			else {
				throw std::runtime_error("size of vec1 is not compatible with that of vec2!");
			}
		}
		/**\brief minus two matrices
		  *\param[in] matrix1,matrix2 the input matrices of the same size
		  *\param[out] out_mat result of input matrices
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void minus(const std::vector<std::vector<NumType, allocator> > &matrix1, const std::vector<std::vector<NumType, allocator> > &matrix2, std::vector<std::vector<NumType, allocator> > &out_mat)
		{
			out_mat.clear();
			if (matrix1.size() == 0 && matrix2.size() == 0) {

			}
			else if (matrix1.size() == matrix2.size()) {
				for (auto it1 = matrix1.begin(), it2 = matrix2.begin(); it1 != matrix1.end(); ++it1, ++it2) {
					out_mat.push_back(std::vector<NumType, allocator>());
					minus(*it1, *it2, *(out_mat.end() - 1));
				}
			}
			else {
				throw std::runtime_error("size of matrix1 is not compatible with that of matrix2!");
			}
		}
		/**\brief compute the dot product of two vectors
		  *\param[in] vec1,vec2 the input vector of the same size
		  *\param[out] out_vec the output vector, which is the dot product of input vectors
		  */
		template<class NumType,class allocator=std::allocator<NumType> >
		void dot(const std::vector<NumType, allocator> &vec1, const std::vector<NumType, allocator> &vec2, std::vector<NumType, allocator> &out_vec)
		{
			out_vec.clear();
			if (vec1.size() == 0 && vec2.size() == 0) {
				//do nothing if sizes of both vectors are zero
			}
			else if (vec1.size() == vec2.size()) {
				for (auto it1 = vec1.begin(), it2 = vec2.begin(); it1 != vec1.end(); ++it1, ++it2) {
					out_vec.push_back((*it1)*(*it2));
				}
			}
			else {
				throw std::runtime_error("size of vec1 is not compatible with that of vec2!");
			}
		}
		/**\brief compute the dot product of two matrices
		*\param[in] matrix1,matrix2 the input matrices of the same size
		*\param[out] out_mat the output matrix, which is the dot product of input matrices
		*/
		template<class NumType, class allocator = std::allocator<NumType> >
		void dot(const std::vector<std::vector<NumType, allocator> > &matrix1, const std::vector<std::vector<NumType, allocator> > &matrix2, std::vector<std::vector<NumType, allocator> > &out_mat)
		{
			out_mat.clear();
			if (matrix1.size() == 0 && matrix2.size() == 0) {

			}
			else if (matrix1.size() == matrix2.size()) {
				for (auto it1 = matrix1.begin(), it2 = matrix2.begin(); it1 != matrix1.end(); ++it1, ++it2) {
					out_mat.push_back(std::vector<NumType, allocator>());
					dot(*it1, *it2, *(out_mat.end() - 1));
				}
			}
			else {
				throw std::runtime_error("size of matrix1 is not compatible with that of matrix2!");
			}
		}
		/**\brief compute the cross product of two vectors
		  *\param[in] vec1,vec2 the input vector of the same size
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		NumType multiply(const std::vector<NumType, allocator> &vec1, const std::vector<NumType, allocator> &vec2)
		{
			NumType sum = 0;
			if (vec1.size() == 0 || vec2.size() == 0) {
				//do nothing if sizes of both vectors are zero
			}
			else if (vec1.size() == vec2.size()) {
				for (auto it1 = vec1.begin(), it2 = vec2.begin(); it1 != vec1.end(); ++it1, ++it2) {
					sum += (*it1)*(*it2);
				}
			}
			else {
				throw std::runtime_error("size of vec1 is not compatible with that of vec2!");
			}
			return sum;
		}
		/**\brief compute the cross product of two vectors, consider first vector as column vector and second as row vector
		  *\param[in] vec1,vec2 the input vector of the same size
		  *\param[out] result the resulting matrix of the cross product
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void multiply(const std::vector<NumType, allocator> &vec1, const std::vector<NumType, allocator> &vec2, std::vector<std::vector<NumType, allocator> > &result)
		{
			result.clear();
			result.resize(vec1.size());
			for (int i = 0; i < result.size(); ++i) {
				result[i].resize(vec2.size());
				for (int j = 0; j < vec2.size(); ++j) {
					result[i][j] = vec1[i] * vec2[j];
				}
			}
		}
		/**\brief compute the cross product of two matrices(underlying assumption:size of all second dimensions are the same)
		  *\param[in] matrix1,matrix2 the input matrices of the same size
		  *param[out] out_mat the output matrix of the cross product of the two input matrices
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void multiply(const std::vector<std::vector<NumType, allocator> > &matrix1, const std::vector<std::vector<NumType, allocator> > &matrix2, std::vector<std::vector<NumType, allocator> > &out_mat)
		{
			out_mat.clear();
			for (int row = 0; row < matrix1.size(); ++row) {
				out_mat.push_back(std::vector<NumType,allocator>());
				for (int col = 0; col < matrix2[0].size(); ++col) {
					out_mat[row].push_back(0);
					NumType &sum = out_mat[row][col];
					auto it1 = matrix1[row].begin();
					auto it2 = matrix2.begin();
					for (; it1 != matrix1[row].end(); ++it1, ++it2) {
						sum += *it1*(*it2)[col];
					}
				}
			}
		}
		/**\brief apply a function to each element of a vector,eg.dst=f(src)
		  *\param[in] src the input vector
		  *\param[out] dst the output vector
		  *\param,[in] f the applied function,f should take in a NumType1 argument and return a NumType2 value
		  */
		template<class NumType1, class NumType2, class allocator1 = std::allocator<NumType1>,class allocator2=std::allocator<NumType2>, class Func>
		void apply(const std::vector<NumType1, allocator1> &src, std::vector<NumType2, allocator2> &dst,Func f)
		{
			dst.clear();
			dst.resize(src.size());
			auto it1 = src.begin();
			auto it2 = dst.begin();
			for (; it1 != src.end(); ++it1, ++it2) {
				*it2 = f(*it1);
			}
		}
		/**\brief apply a function to each element of a vector,eg.f(src)
		  *\param[in] src the input vector
		  *\param,[in] f the applied function,f should take in a NumType argument
		  */
		template<class NumType, class allocator = std::allocator<NumType>, class Func>
		void apply(const std::vector<NumType, allocator> &src, Func f)
		{
			auto it = src.begin();
			for (; it != src.end(); ++it) {
				f(*it);
			}
		}
		/**\brief apply a function to each element of a matrix,eg.dst=f(src)
		  *\param[in] src the input matrix
		  *\param[out] dst the output matrix
		  *\param,[in] f the applied function,f should take in a NumType1 argument and return a NumType2 value
		  */
		template<class NumType1, class NumType2, class allocator1 = std::allocator<NumType1>, class allocator2 = std::allocator<NumType2>, class Func>
		void apply(const std::vector<std::vector<NumType1, allocator1> > &src, std::vector<std::vector<NumType2, allocator2> > &dst, Func f)
		{
			dst.clear();
			dst.resize(src.size());
			auto it1 = src.begin();
			auto it2 = dst.begin();
			for (; it1 != src.end(); ++it1, ++it2) {
				apply(*it1, *it2, f);
			}
		}
		/**\brief apply a function to each element of a matrix,eg.f(src)
		  *\param[in] src the input matrix
		  *\param,[in] f the applied function,f should take in a NumType argument
		  */
		template<class NumType, class allocator = std::allocator<NumType>, class Func>
		void apply(const std::vector<std::vector<NumType, allocator> > &src, Func f)
		{
			auto it = src.begin();
			for (; it != src.end(); ++it) {
				apply(*it, f);
			}
		}
		/**\brief find the elements in a vector that meet a certain requirement
		  *\param[in] src the given vector
		  *\param[out] indices the indices of elements that satisfy the requirement(s)
		  *\param[in] f the condition function, this function takes in a NumType argument and returns a boolean value that decides whether the NumType element satisfy the requirement(s)
		  */
		template<class NumType, class allocator = std::allocator<NumType>, class Func>
		void where(const std::vector<NumType, allocator> &src, std::vector<int> &indices, Func f)
		{
			indices.clear();
			for (int i = 0; i < src.size(); ++i) {
				if (f(src[i])) {
					indices.push_back(i);
				}
			}
		}
		/**\brief find the elements in a matrix that meet a certain requirement
		  *\param[in] src the given matrix
		  *\param[out] indix_x,index_y the indices of elements that satisfy the requirement(s)
		  *\param[in] f the condition function, this function takes in a NumType argument and returns a boolean value that decides whether the NumType element satisfy the requirement(s)
		  */
		template<class NumType, class allocator = std::allocator<NumType>, class Func>
		void where(const std::vector<std::vector<NumType, allocator> > &src, std::vector<int> &index_x, std::vector<int> &index_y, Func f)
		{
			index_x.clear();
			index_y.clear();
			for (int i = 0; i < src.size(); ++i) {
				for (int j = 0; j < src[i].size(); ++j) {
					if (f(src[i][j])) {
						index_x.push_back(i);
						index_y.push_back(j);
					}
				}
			}
		}
		/**\brief find the elements in a vector that meet a certain requirement and apply function to it
		  *\param[in] src the given vector
		  *\param[in] f the condition function, this function takes in a NumType argument and returns a boolean value that decides whether the NumType element satisfy the requirement(s)
		  *\param[in] pro the applied function
		  */
		template<class NumType, class allocator = std::allocator<NumType>, class Pred,class Process>
		void apply_where(std::vector<NumType, allocator> &src, Pred f, Process f2)
		{
			for (int i = 0; i < src.size(); ++i) {
				if (f(src[i])) {
					f2(src[i]);
				}
			}
		}
		/**\brief find the elements in a matrix that meet a certain requirement and apply function to it
		  *\param[in] src the given matrix
		  *\param[in] f the condition function, this function takes in a NumType argument and returns a boolean value that decides whether the NumType element satisfy the requirement(s)
		  *\param[in] pro the applied function
		  */
		template<class NumType, class allocator = std::allocator<NumType>, class Pred, class Process>
		void apply_where(std::vector<std::vector<NumType, allocator> > &src, Pred f, Process f2)
		{
			for (int i = 0; i < src.size(); ++i) {
				for (int j = 0; j < src[i].size(); ++j) {
					if (f(src[i][j])) {
						f2(src[i][j]);
					}
				}
			}
		}
	}
	/**\brief wrap up some Eigen functionalities
	 * \this is poorly designed. don't use it unless you have nothing else to use
	  */
	namespace eigen_wrapper
	{
		/**\brief find the covariance matrix of samples
		  *\param[in] samples the sample data(a n*m 2d vector,n represents number of features and m represents number of samples)
		  *\param[out] cov the covariance matrix of the given sample data
		  *\the covariance matrix is calculated as cov=(X-mu)@(X-mu)',where X represent the samples(n*m matrix) and mu is the mean values of X(mu=1/m*X);@ represents matrix multiplication;X' represents the transpose of X
		  *\note:above is the covariance matrix of the samples. the covariance matrix of the population is estimated as Var(x)=1/n*cov(x): as maximum likelihood estimation would tell you
		  */
		template<class NumType,class allocator=std::allocator<NumType> >
		void cov_samples(const std::vector<std::vector<NumType, allocator> > &samples, std::vector<std::vector<NumType, allocator> > &cov)
		{
			std::vector<NumType, allocator> *m = new std::vector<NumType, allocator>;
			std::vector<std::vector<NumType, allocator> > *m_broad = new std::vector<std::vector<NumType, allocator> >;
			std::vector<std::vector<NumType, allocator> > *Xi = new std::vector<std::vector<NumType, allocator> >;
			std::vector<std::vector<NumType, allocator> > *Xi_trans = new std::vector<std::vector<NumType, allocator> >;
			int num_feat = samples.size();
			int num_samp = samples[0].size();
			numpy::mean(samples, *m, 1);
            numpy::broadcast(*m, num_samp, *m_broad, false);
            numpy::minus(samples, *m_broad, *Xi);
            numpy::transpose(*Xi, *Xi_trans);
            numpy::multiply(*Xi, *Xi_trans, cov);
			delete m;
			delete m_broad;
			delete Xi;
			delete Xi_trans;
		}
        		/**\brief convert between 3x3 2d vectors and 3x3 Matrix3d type defined in Eigen(this is pretty much a useless wrapper up)
		  *\param[in] vec2d the input/output 2d vector(if convert the vector to Matrix3d,then the dimension of the vector should be 3x3)
		  *\parma[in] mat the input/output Matrix3d object
		  *\param[in] vec2mat if true,convert the vector to Matrix3d;else convert the Matrix3d to vector
		  */
		template<class NumType, class allocator = std::allocator<NumType> >
		void convert3d(std::vector<std::vector<NumType, allocator> > &vec2d, Eigen::Matrix3d &mat, bool vec2mat = true)
		{
			if (vec2mat) {
				if (vec2d.size() != 3 || vec2d[0].size() != 3 || vec2d[1].size() != 3 || vec2d[2].size() != 3) {
					throw std::runtime_error("this function \"conversion\" can only convert 3x3 2d vectors to Matrix3d or vice versa");
				}
			}
			if (vec2mat) {
				mat << vec2d[0][0], vec2d[0][1], vec2d[0][2],
					vec2d[1][0], vec2d[1][1], vec2d[1][2],
					vec2d[2][0], vec2d[2][1], vec2d[2][2];
			}
			else {
				vec2d.clear();
				vec2d.resize(3);
				for (int i = 0; i < 3; ++i) {
					for (int j = 0; j < 3; ++j) {
						vec2d[i].push_back(mat(i, j));
					}
				}
			}
		}
        /**\brief find the covariance matrix of a collection of points
          *\param[in] pList the point list
          *\param[out] cov the covariance matrix
          */
        template<class PointXYZ, class allocator = std::allocator<PointXYZ> >
        void cov_points(const std::vector<PointXYZ, allocator> &pList, Eigen::Matrix3d &cov)
        {
            std::vector<std::vector<double> > *samples = new std::vector<std::vector<double> >;
            for (int i = 0; i < 3; ++i) { samples->push_back(std::vector<double>()); }
            numpy::apply(pList, [samples](PointXYZ p) {
                (*samples)[0].push_back(p.x);
                (*samples)[1].push_back(p.y);
                (*samples)[2].push_back(p.z);
            });
            std::vector<std::vector<double> > cov_vec;
            cov_samples(*samples, cov_vec);
            convert3d(cov_vec, cov);
            delete samples;
        }
		/**\brief get the eigen values and eigen vectors of a given 3x3 matrix(this is only a wrapper of corrsponding Eigen function)
		  *\definition: v'@m@v=eigs,where v,v' are the eigen vector matrix and its transpose, m is the original matrix and eigs is the eigen value matrix
		  *\param[in] mat the given matrix
		  *\param[out] eigenValues the eigen values diagonal matrix, the eigen values are in descend order
		  *\param[out] eigenVectors the eigen vector matrix(each column is a eigen vector corresponding to eigen value at the same column,eigen vector is normalized)
		  */
		void eig3d(const Eigen::Matrix3d &mat, Eigen::Matrix3d &eigenValues, Eigen::Matrix3d &eigenVectors);
	}
}
namespace np = lyxutils::numpy;
namespace alg = lyxutils::algorithms;
namespace dbg = lyxutils::debug;
namespace lin = lyxutils::eigen_wrapper;//lin for linear algebra
#endif // !UTILS_H