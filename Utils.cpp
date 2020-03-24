//#include "stdafx.h"
#include "Utils.h"
#ifdef _WIN32
#include <io.h>
#include <direct.h>
#endif
#ifdef linux
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif
#include <fstream>
#include <sstream>
#include <exception>

std::vector<std::string> lyxutils::str_utils::split(const std::string & str, const std::string & delimiter)
{
	std::vector<std::string> result = std::vector<std::string>();
	if (str == "") return result;
	char *chs = new char[str.length() + 1];
	std::strcpy(chs, str.c_str());
	char *d = new char[delimiter.length() + 1];
	std::strcpy(d, delimiter.c_str());
	char *p = std::strtok(chs, d);
	while (p) {
		std::string s = p;
		result.push_back(s);
		p = std::strtok(NULL, d);
	}
	delete chs, d;
	return result;
}

std::vector<std::string> lyxutils::str_utils::split_with_string(const std::string & str, const std::string & delString)
{
	std::vector<std::string> result = std::vector<std::string>();
	if (str == "") return result;
	size_t start = 0;
	size_t end = str.find(delString, start);
	while (end < str.length()) {
		if(end>start)result.push_back(str.substr(start, end - start));
		start = end + delString.length();
		end = str.find(delString, start);
	}
	if (start < str.length()) {
		result.push_back(str.substr(start, str.length() - start));
	}
	return result;
}

std::string lyxutils::str_utils::lower(const std::string & str)
{
	std::string res = str;
	std::transform(res.begin(), res.end(), res.begin(), ::tolower);
	return res;
}

std::string lyxutils::str_utils::upper(const std::string & str)
{
	std::string res = str;
	std::transform(res.begin(), res.end(), res.begin(), ::toupper);
	return res;
}

std::string lyxutils::str_utils::replace(const std::string &str, const std::string &oldStr, const std::string &newStr, int count)
{
	std::string res = str;
	size_t pos = res.find(oldStr, 0);
	while (count--!=0 && pos < res.length()) {
        res.replace(pos, oldStr.length(), newStr);
        pos = res.find(oldStr, pos + newStr.length());
    }
	return res;
}

std::string lyxutils::str_utils::replace(const std::string &str, const std::string &oldStr, const std::string &newStr, int start, int end) {
    if(start<0)start+=str.length();
    if(end<0)end+=str.length();
    if(start<0||start>str.length()||end<0||end>str.length())throw std::out_of_range("string index out of range");
    if(end<start)throw std::logic_error("end position at left side of start position is not allowed");
    std::string res=str.substr(start,end-start);
    size_t pos = res.find(oldStr, 0);
    while (pos < res.length()) {
        res.replace(pos, oldStr.length(), newStr);
        pos = res.find(oldStr, pos + newStr.length());
    }
    res=str.substr(0,start)+res+str.substr(end,str.length()-end);
    return res;
}

int lyxutils::str_utils::count(const std::string &str, const std::string &sub, int start, int end) {
    if(start<0)start+=str.length();
    if(end<0)end+=str.length();
    if(start<0||start>str.length()||end<0||end>str.length())throw std::out_of_range("string index out of range");
    if(end<start)throw std::logic_error("end position at left side of start position is not allowed");
    std::string temp=str.substr(start,end-start);
    size_t pos=temp.find(sub,0);
    int count=0;
    while(pos<temp.length()){
        count++;
        pos=temp.find(sub,pos+sub.length());
    }
    return count;
}

std::string lyxutils::str_utils::multiply(const std::string &str, int count) {
    std::string result=str;
    while(--count>0){
        result+=str;
    }
    return result;
}

std::string lyxutils::str_utils::center(const std::string &str, int width, char c) {
    std::string result=str;
    if(width<=result.length())return result;
    int num=(width-str.length())/2;
    std::string pad=multiply(type2str(c),num);
    result=pad+result+pad;
    if(result.length()<width)result+=type2str(c);
    return result;
}

std::string lyxutils::str_utils::strip(const std::string &str, const std::string &chars) {
    auto it=str.begin();
    for(;it!=str.end();++it){
        if(lyxutils::str_utils::count(chars, type2str(*it), 0, chars.length())==0)break;
    }
    auto it2=str.end()-1;
    for(;it2!=it-1;--it2){
        if(lyxutils::str_utils::count(chars, type2str(*it2), 0, chars.length())==0)break;
    }
    std::string result=str.substr(it-str.begin(),it2+1-it);
    return result;
}

std::string lyxutils::str_utils::frame(const std::string &title, const std::string &text, int width, int pad, char border, int shift) {
    return std::__cxx11::string();
}

bool lyxutils::io::read_csv(const std::string &fileName, std::vector<std::vector<float> > &table, const std::string &sep, std::string &report,
                            std::vector<std::string> &headers, int fields, bool read_header, int numberOfRedirect, void(*directFunc)(const std::vector<float>&))
{
	if (!lyxutils::io::fileExists(fileName))
		throw std::runtime_error("文件不存在!");
	std::fstream inFile(fileName);
	std::stringstream ssReport;
	ssReport << lyxutils::debug::getTimeString() << ":打开文件 " << fileName << std::endl;
	ssReport << "文件报告:" << std::endl;
	std::string line;
	bool first_time = true;
	int count = 0;
	int error_count = 0;
	while (std::getline(inFile, line)) {
		if (line == "")continue;
		count++;
		std::vector<std::string> field_values = lyxutils::str_utils::split(line, sep);
		std::for_each(field_values.begin(), field_values.end(), [](std::string &str) {str = lyxutils::str_utils::replace(
                str, " ", "", 0); });
		if (first_time) {
			first_time = false;
			fields = fields > field_values.size() ? (int)(field_values.size()) : fields;
			if (fields == -1)
				fields = field_values.size();
			ssReport << "打开的文件共有" << field_values.size() << "个字段,读取了其中前" << fields << "个字段" << std::endl;
			table.clear();
			for (int i = 0; i < fields - numberOfRedirect; ++i) {
				table.push_back(std::vector<float>());
			}
			if (read_header) {
				headers = field_values;
				count--;
				continue;
			}
		}
		try {
			std::vector<float> line_data = std::vector < float>(fields);
			for (int i = 0; i < fields; ++i) {
				if (i < field_values.size())
					line_data[i] = std::stof(field_values[i]);//将字符串转为float
				else {
					line_data[i] = 0;//如果实际读取到的字段数不足要求的字段数，则不足的填为0
					if (i == field_values.size()) {
						ssReport << "第" << count << "行数据列数不足" << fields << ",只有" 
							<< field_values.size() << "列，不足的列已补为0" << std::endl;
					}
				}	
			}
			std::vector<float> redirectFields = std::vector<float>();
			for (int i = 0; i < fields; ++i) {
				if (i < numberOfRedirect) {
					redirectFields.push_back(line_data[i]);
				}
				else {
					table[i - numberOfRedirect].push_back(line_data[i]);
				}
			}
			if (numberOfRedirect > 0 && directFunc != NULL) {
				directFunc(redirectFields);
			}
		}
		catch (std::exception) {
			error_count++;
			ssReport << "第" << count << "行数据出错,已忽略;该行数据为:";
			std::for_each(field_values.begin(), field_values.end(), [&ssReport](std::string str) {ssReport << str << " "; });
			ssReport << std::endl;
		}
	}
	if (inFile.is_open()) {
		inFile.close();
	}
	ssReport << "读取文件共有" << count << "行非空数据行,其中错误" << error_count << "行" << std::endl;
	report = ssReport.str();
	return true;
}

bool lyxutils::io::fileExists(std::string fileName)
{
	return access(fileName.c_str(), 0) == 0;
}

void lyxutils::colors::rgb2hsv(const std::vector<int>& rgb, std::vector<float>& hsv)
{
	int r = rgb[0];
	int g = rgb[1];
	int b = rgb[2];
	int cmax = r > g ? r : g;
	cmax = b > cmax ? b : cmax;
	int cmin = r < g ? r : g;
	cmin = b < cmin ? b : cmin;
	hsv.clear();
	for (auto i : { 0,1,2 })hsv.push_back(0);
	if (cmax != cmin)
	{
		if (cmax == r) {
			hsv[0] = 60 * (float)(g - b) / (cmax - cmin);
		}
		else if (cmax == g) {
			hsv[0] = 120 + 60 * (float)(b - r) / (cmax - cmin);
		}
		else {
			hsv[0] = 240 + 60 * (float)(r - g) / (cmax - cmin);
		}
		if (hsv[0] < 0)hsv[0] += 360;
	}
	if (cmax != 0)hsv[1] = float(cmax - cmin) / cmax;
	hsv[2] = cmax;
}

void lyxutils::colors::hsv2rgb(const std::vector<float>& hsv, std::vector<int>& rgb)
{
	float h = hsv[0];
	float s = hsv[1];
	float v = hsv[2];
	int hi = int(h / 60) % 6;
	float f = h / 60 - hi;
	float p = v*(1 - s);
	float q = v*(1 - f*s);
	float t = v*(1 - (1 - f)*s);
	rgb.clear();
	for (auto i : { 1,1,1 })rgb.push_back(0);
	switch (hi)
	{
	case 0:rgb[0] = (int)(v); rgb[1] = (int)(t); rgb[2] = (int)(p); break;
	case 1:rgb[0] = (int)(q); rgb[1] = (int)(v); rgb[2] = (int)(p); break;
	case 2:rgb[0] = (int)(p); rgb[1] = (int)(v); rgb[2] = (int)(t); break;
	case 3:rgb[0] = (int)(p); rgb[1] = (int)(q); rgb[2] = (int)(v); break;
	case 4:rgb[0] = (int)(t); rgb[1] = (int)(p); rgb[2] = (int)(v); break;
	case 5:rgb[0] = (int)(v); rgb[1] = (int)(p); rgb[2] = (int)(q); break;
	}
}

double lyxutils::algorithms::distanceToPlane(double x, double y, double z, double a, double b, double c, double d)
{
	return abs(a*x + b*y + c*z + d);
}

std::string lyxutils::debug::getTimeString()
{
	time_t timep;
	time(&timep);
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S", localtime(&timep));
	return tmp;
}

bool lyxutils::debug::createFolder(const std::string folderPath)
{
	if (access(folderPath.c_str(), 0) == -1) {
#ifdef _WIN32
	    bool created=_mkdir(folderPath.c_str());
#endif
#ifdef linux
	    bool created=mkdir(folderPath.c_str(),775);
#endif
		if (created == 0) {
			return true;
		}
		else {
			return false;
		}
	}
	return true;
}

std::ofstream * lyxutils::debug::openLogFile(std::string fileName)
{
	std::ofstream *log = new std::ofstream();
	log->open(fileName, std::ios::app);
	return log;
}

void lyxutils::debug::writeLog(std::ofstream & log, std::string message, std::string end)
{
	if (log.is_open()) {
		log << lyxutils::debug::getTimeString() << "\t" << message << end;
	}
}

void lyxutils::debug::consoleProgressBar(int value, std::string delimiterL, std::string increment, std::string delimiterR, std::string status)
{
	if (value == 0) {
		std::cout << std::endl;
	}
	std::cout << "\r";
	std::cout << delimiterL;
	for (int i = 0; i < value; ++i) {
		std::cout << increment;
	}
	for (int i = value; i < 100; ++i) {
		std::cout << " ";
	}
	std::cout << delimiterR;
	std::cout << value << "% ";
	if (status == "" || value == 100) {
		std::cout << "completed";
	}
	else {
		std::cout << status;
	}
	if (value == 100) {
		std::cout << std::endl;
	}
}

void lyxutils::eigen_wrapper::eig3d(const Eigen::Matrix3d & mat, Eigen::Matrix3d & eigenValues, Eigen::Matrix3d & eigenVectors)
{
	Eigen::EigenSolver<Eigen::Matrix3d> es(mat);
	std::vector<double> idx;
	np::range<double>(idx, 3);
	eigenValues = es.pseudoEigenvalueMatrix();
	eigenVectors = es.pseudoEigenvectors();
	std::sort(idx.begin(), idx.end(), [&eigenValues](int pre, int post)->bool {return eigenValues(pre, pre) > eigenValues(post, post); });
	std::vector<std::vector<double> > elementary, eig, vec, eigOrdered, eigTemp, vecOrdered;
	np::broadcast<double>(0, 3, 3, elementary);
	for (int i = 0; i < idx.size(); ++i) {
		elementary[i][idx[i]] = 1;
	}
	//the stupid Eigen can not multiply two matrix. reported some C2338 error:YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX
	//and thus I have to use my self-defined multiplication for 2d vectors, thank god it's only 3x3 matrix
	lin::convert3d(eig, eigenValues, false);
	lin::convert3d(vec, eigenVectors, false);
	np::multiply(eig, elementary, eigTemp);
	np::multiply(elementary, eigTemp, eigOrdered);
	np::multiply(vec, elementary, vecOrdered);
	lin::convert3d(eigOrdered, eigenValues);
	lin::convert3d(vecOrdered, eigenVectors);
}
