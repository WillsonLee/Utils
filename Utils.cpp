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
    if(count<0)throw std::invalid_argument("times of repetition should not be negative!");
    std::string result;
    while(count-->0){
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
    if(pad<0)throw std::invalid_argument("pad should not be negative!");
    if(width<MIN_WORDWRAP_WIDTH+2+2*pad)throw std::invalid_argument("width should be at least 1+pad+MIN_WORDWRAP_WIDTH+pad+1!");
    if(title.length()>width)throw std::invalid_argument("width of title should be less than width!");
    std::vector<std::string> lines=word_wrap(text,width-2*pad-2);
    shift=shift<0?0:shift;
    std::string prefix=multiply(" ",shift);
    std::string top=prefix+center(title,width,border);
    std::string gap=prefix+border+center("",width-2,' ')+border;
    std::string bottom=prefix+center("",width,border);
    lines.insert(lines.begin(),gap);
    lines.insert(lines.begin(),top);
    lines.push_back(gap);
    lines.push_back(bottom);
    for(int i=2;i<lines.size()-2;++i){
        lines[i]=prefix+border+center(lines[i],width-2,' ')+border;
    }
    return join("\n",lines.begin(),lines.end());
}

std::string lyxutils::str_utils::left(const std::string &str, int width, char c) {
    std::string result=str;
    if(result.length()<width){
        result+=multiply(std::string(1,c),width-result.length());
    }
    return result;
}

std::string lyxutils::str_utils::right(const std::string &str, int width, char c) {
    std::string result=str;
    if(result.length()<width){
        result=multiply(std::string(1,c),width-result.length())+result;
    }
    return result;
}

std::vector<std::string> lyxutils::str_utils::word_wrap(const std::string &str, int width) {
    if(width<MIN_WORDWRAP_WIDTH)throw std::invalid_argument("width should be greater than MIN_WORDWRAP_WIDTH");
    std::vector<std::string> res;
    //get word wrap result of a single paragraph
    auto wrap_single_para=[&](const std::string &str, int width)->std::vector<std::string>{
        std::vector<std::string> result;
        if(str.length()<=width)result.push_back(lyxutils::str_utils::left(str,width,' '));
        else{
            auto isLetter=[](char c)->int{if((c>=65&&c<91)||(c>=97&&c<123))return 1;else return 0;};
            std::string line;
            for(int i=0;i<str.length();){
                if(line.length()==width-1&&isLetter(str[i])){
                    if(isLetter(*(line.end()-1))) {
                        line.push_back('-');
                    }
                    else{
                        line.push_back(' ');
                    }
                }
                else{
                    line.push_back(str[i++]);
                }
                if(line.length()==width||i==str.length()) {
                    result.push_back(lyxutils::str_utils::left(line,width,' '));
                    line.clear();
                }
            }
        }
        return result;
    };
    std::vector<std::string> paras=split(str,"\n");
    for(int i=0;i<paras.size();++i){
        std::vector<std::string> wrap_p=wrap_single_para(paras[i],width);
        for_each(wrap_p.begin(),wrap_p.end(),[&](const std::string &line){res.push_back(line);});
    }
    return res;
}

bool lyxutils::io::read_csv(const std::string &fileName, std::vector<std::vector<float> > &table, const std::string &sep, std::string &report,
                            std::vector<std::string> &headers, int fields, bool read_header, int numberOfRedirect, void(*directFunc)(const std::vector<float>&))
{
	if (!lyxutils::io::fileExists(fileName))
		throw std::runtime_error("file not exist!");
	std::fstream inFile(fileName);
	std::stringstream ssReport;
	ssReport << lyxutils::debug::getTimeString() << ":opened file " << fileName << std::endl;
	ssReport << "file report:" << std::endl;
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
	ssReport << "total non-empty line of read file:" << count << ",error line:" << error_count << std::endl;
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
	numpy::range<double>(idx, 3);
	eigenValues = es.pseudoEigenvalueMatrix();
	eigenVectors = es.pseudoEigenvectors();
	std::sort(idx.begin(), idx.end(), [&eigenValues](int pre, int post)->bool {return eigenValues(pre, pre) > eigenValues(post, post); });
	std::vector<std::vector<double> > elementary, eig, vec, eigOrdered, eigTemp, vecOrdered;
	numpy::broadcast<double>(0, 3, 3, elementary);
	for (int i = 0; i < idx.size(); ++i) {
		elementary[i][idx[i]] = 1;
	}
	//the stupid Eigen can not multiply two matrix. reported some C2338 error:YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX
	//and thus I have to use my self-defined multiplication for 2d vectors, thank god it's only 3x3 matrix
	convert3d(eig, eigenValues, false);
	convert3d(vec, eigenVectors, false);
	numpy::multiply(eig, elementary, eigTemp);
	numpy::multiply(elementary, eigTemp, eigOrdered);
	numpy::multiply(vec, elementary, vecOrdered);
	convert3d(eigOrdered, eigenValues);
	convert3d(vecOrdered, eigenVectors);
}

bool lyxutils::io::createFolder(const std::string folderPath)
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

std::ofstream * lyxutils::io::openLogFile(std::string fileName)
{
    std::ofstream *log = new std::ofstream();
    log->open(fileName, std::ios::app);
    return log;
}

void lyxutils::io::writeLog(std::ofstream & log, std::string message, std::string end)
{
    if (log.is_open()) {
        log << lyxutils::debug::getTimeString() << "\t" << message << end;
        log.flush();
    }
}

lyxutils::io::CLParser::CLParser() {

}

void lyxutils::io::CLParser::parse(std::string command) {
    std::vector<std::string> argv=str_utils::split(command," \t");
    parse(argv);
}

void lyxutils::io::CLParser::parse(int argc, char **argv) {
    std::vector<std::string> cmd;
    for(int i=0;i<argc;++i){
        cmd.push_back(argv[i]);
    }
    parse(cmd);
}

void lyxutils::io::CLParser::parse(std::string commands, std::string filename) {
    if (!lyxutils::io::fileExists(filename))throw std::runtime_error("file not exist!");
    std::fstream inFile(filename);
    std::string cl=commands;
    if(inFile.is_open()){
        std::string line;
        while(std::getline(inFile,line)){
            cl+=" "+line;//make sure there is a space between two lines
        }
        inFile.close();
    }
    else{
        throw std::runtime_error("open file " + filename + " failed");
    }
    cl=str_utils::replace(cl,"\n"," ");
    parse(cl);
}

void lyxutils::io::CLParser::setOptionPattern(int num_opts, const lyxutils::io::CLOption *options) {
    _opt_patterns.clear();
    _short2verbose.clear();
    for(int i=0;i<num_opts;++i){
        _opt_patterns[options[i].optName]=options[i];
        if(options[i].abbr) {
            _short2verbose[str_utils::type2str(*options[i].abbr)]=options[i].optName;
        }
    }
    _pattern_set=true;
}

void lyxutils::io::CLParser::clearOptionPattern() {
    _pattern_set=false;
    _opt_patterns.clear();
    _short2verbose.clear();
}

bool lyxutils::io::CLParser::hasOption(std::string optName) {
    return _opt_val.find(optName)!=_opt_val.end();
}

std::string lyxutils::io::CLParser::getOptionArg(std::string optName) {
    return _opt_val[optName];
}

std::vector<std::string> lyxutils::io::CLParser::getParameters() {
    return _parameters;
}

std::map<std::string, std::string> lyxutils::io::CLParser::getAllOptions() {
    return _opt_val;
}

void lyxutils::io::CLParser::parse(const std::vector<std::string> &argv) {
    _opt_val.clear();
    _parameters.clear();
    auto isShort=[](std::string chuck)->bool{return chuck.size() >= 2 && chuck[0] == '-' && chuck[1]!='-';};
    auto isLong=[](std::string chuck)->bool{return chuck.size() > 2 && chuck.substr(0, 2) == "--";};
    auto isParameter=[](std::string chuck)->bool{return chuck[0]!= '-';};
    int argc=argv.size();
    for(int i=1;i<argc;++i) {
        std::string chuck(argv[i]);
        std::vector<std::string> sv=str_utils::split(chuck,"=");
        std::string opt=sv[0];
        if (isShort(opt)) {//short option
            if(_pattern_set){
                if(_short2verbose.find(opt.substr(1,opt.length()-1))!=_short2verbose.end()){
                    opt=sv[0]="--"+_short2verbose[opt.substr(1,opt.length()-1)];
                    //handle the converted opt to if(isLong(opt)){...} to process, so no continue here
                }
                else{
                    throw std::invalid_argument("unrecognized option:"+opt);
                }
            }
            else{
                if(sv.size()>1){
                    _opt_val[opt.substr(1,opt.length()-1)]=sv[1];
                }
                else if(i+1<argc&&isParameter(argv[i+1])){
                    _opt_val[opt.substr(1,opt.length()-1)]=std::string(argv[i+1]);
                    ++i;//next element has been used,need to skip it in next loop
                }
                else{
                    _opt_val[opt.substr(1,opt.length()-1)]="";//no argument if not found
                }
                continue;
            }
        }
        if (isLong(opt)) {//long option
            if(_pattern_set){
                if(_opt_patterns.find(opt.substr(2,opt.length()-2))!=_opt_patterns.end()){
                    if(_opt_patterns[opt.substr(2,opt.length()-2)].hasArg==2){
                        if(sv.size()>1){
                            _opt_val[opt.substr(2,opt.length()-2)]=sv[1];
                        }
                        else if(i+1<argc&&isParameter(argv[i+1])){
                            _opt_val[opt.substr(2,opt.length()-2)]=std::string(argv[i+1]);
                            ++i;//next element has been used,need to skip it in next loop
                        }
                        else{
                            std::string full_name=opt.substr(2,opt.length()-2);
                            std::string type_hint="--"+full_name;
                            if(_opt_patterns[full_name].abbr)type_hint+="("+str_utils::type2str(_opt_patterns[full_name].abbr)+")";
                            throw std::invalid_argument("option "+type_hint+" has no argument!");
                        }
                        continue;
                    }
                    else if(_opt_patterns[opt.substr(2,opt.length()-2)].hasArg==1){
                        if(sv.size()>1){
                            _opt_val[opt.substr(2,opt.length()-2)]=sv[1];
                        }
                        else{
                            _opt_val[opt.substr(2,opt.length()-2)]="";
                        }
                    }
                    else{
                        _opt_val[opt.substr(2,opt.length()-2)]="";//hasArg=0,no argument required
                    }
                }
                else{
                    throw std::invalid_argument("unrecognized option:"+opt);
                }
            }
            else{
                if(sv.size()>1){
                    _opt_val[opt.substr(2,opt.length()-2)]=sv[1];
                }
                else if(i+1<argc&&isParameter(argv[i+1])){
                    _opt_val[opt.substr(2,opt.length()-2)]=std::string(argv[i+1]);
                    ++i;//next element has been used,need to skip it in next loop
                }
                else{
                    _opt_val[opt.substr(2,opt.length()-2)]="";
                }
                continue;
            }
        }
        else {//parameters
            _parameters.push_back(chuck);
        }
    }
}

std::string lyxutils::io::CLParser::generateHelp(std::string command, bool requireObj, std::string description,
                                                 const std::vector<std::string> &optNames,
                                                 const std::vector<std::string> &optDescription,
                                                 const std::vector<std::string> &optArgName) {
    const int leftCol=30;
    const int rightCol=50;
    if(_pattern_set&&optNames.size()!=_opt_patterns.size()){
        throw std::logic_error("optNames size("+str_utils::type2str(optNames.size())+") != number of defined options("+str_utils::type2str(_opt_patterns.size())+")");
    }
    std::string result="usage: "+command+" [options]... ";
    if(requireObj){
        result+="<object>...\n";
    }else{
        result+="[object]...\n";
    }
    std::vector<std::string> descrips=str_utils::word_wrap(description,80);
    for(std::string line:descrips)result+=line+"\n";
    result+="\noptions:\n";
    for(int i=0;i<optNames.size();++i){
        if(_pattern_set&&_opt_patterns.find(optNames[i])==_opt_patterns.end())throw std::logic_error("optName not contained in option pattern!");
        std::string optField;
        std::string o="    ";
        if(_pattern_set&&_opt_patterns[optNames[i]].abbr){
            o="-"+str_utils::type2str(*_opt_patterns[optNames[i]].abbr)+", ";
        }
        optField+="  "+o+"--"+optNames[i];
        std::string argName="value";
        if(optArgName.size()>i&&optArgName[i]!="")argName=optArgName[i];
        int arg_state=0;
        if(_pattern_set){
            arg_state=_opt_patterns[optNames[i]].hasArg;
        }
        else{
            if(optArgName.size()<=i||optArgName[i]=="")arg_state=0;
            else arg_state=1;
        }
        if(arg_state==2){
            optField+="="+argName;
        }
        else if(arg_state==1){
            optField+="[="+argName+"]";
        }
        if(optField.length()<leftCol-2){
            optField=str_utils::left(optField,leftCol,' ');
        }
        else{
            optField+="\n"+str_utils::multiply(" ",leftCol);
        }
        result+=optField;
        std::vector<std::string> optDes=str_utils::word_wrap(optDescription[i],rightCol);
        for(int i=0;i<optDes.size();++i){
            if(i!=0)result+=str_utils::multiply(" ",leftCol);
            result+=optDes[i]+"\n";
        }
    }
    return result;
}

lyxutils::io::CLOption::~CLOption() {
    if(abbr)delete abbr;
}
