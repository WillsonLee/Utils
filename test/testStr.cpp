//
// Created by lyx on 2020/3/24.
//

#include <Utils.h>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
namespace str=lyxutils::str_utils;
namespace dbg=lyxutils::debug;

int main(int argc,char **argv){
    string s="that day before ththese thieves got caught this town looks like hell";
    cout<<"original string:"<<s<<endl;
    vector<string> ss=str::split(s,"bt");
    cout<<"1.split with 'b' and 't':"<<endl<<"\t";
    dbg::printVector(ss);
    cout<<"2.split with 'th':"<<endl<<"\t";
    vector<string> ss2=str::split_with_string(s,"th");
    dbg::printVector(ss2);
    cout<<"3.split with non existed character 'z':"<<endl<<"\t";
    vector<string> ss3=str::split(s,"z");
    dbg::printVector(ss3);
    cout<<"3+.split considering multiple space as one:"<<endl<<"\t";
    string spstr="cmake  -v --find-package\t-build\t\trelease\t   ..";
    dbg::printVector(str::split(spstr," \t"));
    string lu="HelLo,mY naME Is aLIsA";
    cout<<"4.case conversion test:"<<endl<<"\t"<<lu<<"=====>lower case:"<<str::lower(lu)<<"========>upper case:"<<str::upper(lu)<<endl;
    string rep="I has a book and I has never been to Shanghai";
    cout<<"5.replace test:"<<endl<<"\t";
    cout<<"origin:"<<rep<<"======>"<<str::replace(rep,"has","have")<<endl<<"\t";
    cout<<"replace only first one:"<<str::replace(rep,"has","have",1)<<endl<<"\t";
    cout<<rep.substr(0,rep.length()/2)<<"<====half====>"<<rep.substr(rep.length()/2,rep.length()-1)<<endl<<"\t";
    cout<<"replace first half:"<<str::replace(rep,"has","have",0,rep.length()/2)<<endl<<"\t";
    cout<<"replace second half:"<<str::replace(rep,"has","have",rep.length()/2,-1)<<endl<<"\t";
    cout<<"replace in range (3,-2):"<<str::replace(rep,"has","have",3,-2)<<endl;
    int a=4;
    float b=3.6;
    double c=0.314;
    long d=444947132932109;
    char e='e';
    string st="10";
    cout<<"6.conversion test:"<<endl<<"\t";
    cout<<"int2str:"<<str::type2str(a)<<"\t\t\tstr2int:"<<str::str2type<int>(st)<<endl<<"\t";
    st="9.4";
    cout<<"float2str:"<<str::type2str(b)<<"\t\t\tstr2float:"<<str::str2type<float>(st)<<endl<<"\t";
    st="0.717";
    cout<<"double2str:"<<str::type2str(c)<<"\t\t\tstr2double:"<<str::str2type<double>(st)<<endl<<"\t";
    st="17801016201";
    cout<<"long2str:"<<str::type2str(d)<<"\t\t\tstr2long:"<<str::str2type<long>(st)<<endl<<"\t";
    st="j";
    cout<<"char2str:"<<str::type2str(e)<<"\t\t\tstr2char:"<<str::str2type<char>(st)<<endl;
    vector<int> v;
    for(int i=0;i<10;++i)v.push_back(2*i+1);
    cout<<"7.join test:"<<endl<<"\t";
    cout<<"original vector:";
    dbg::printVector(v);
    cout<<"\tconcatenate with 'aha':"<<str::join("aha",v.begin(),v.end())<<endl<<"\t";
    cout<<"concatenate second to last but one with white space:"<<str::join(" ",v.begin()+1,v.end()-1)<<endl;
    cout<<"8.count test:"<<endl<<"\t";
    string all="Roses are red, violets are blue. Sugar is sweet, and so are you.";
    cout<<"sentence is:"<<all<<endl<<"\t";
    cout<<"number of 'are' is:"<<str::count(all,"are",0,all.length())<<endl<<"\t";
    cout<<"number of 'are' within range(15,41):"<<str::count(all,"are",15,41)<<",the sub sentence is:"<<all.substr(15,41-15)<<endl;
    cout<<"9.alignment test:"<<endl<<"\t";
    std::string grt="I'm groot!!";
    cout<<"repeat x 4:"<<grt<<"=========>"<<str::multiply(grt,4)<<endl<<"\t";
    cout<<"**center**(width 21):"<<str::center(grt,21,'*')<<endl<<"\t";
    cout<<"**center**(width 20):"<<str::center(grt,20,'*')<<endl<<"\t";
    cout<<"**center**(width 10):"<<str::center(grt,10,'*')<<endl<<"\t";
    cout<<"left******(width 20):"<<str::left(grt,20,'*')<<endl<<"\t";
    cout<<"*****right(width 20):"<<str::right(grt,20,'*')<<endl<<endl;
    cout<<"10.strip test:"<<endl<<"\t";
    string st1="   hello thank you  ";
    string res1=str::strip(st1);
    cout<<"original string|"<<st1<<"|strip space=======>|"<<res1<<"|strip 'heou'|"<<str::strip(res1,"heou")<<"|"<<endl;
    cout<<"11.frame and word wrap test:"<<endl<<endl<<"word wrap:"<<endl<<endl;
    string txt="  This file is part of the GNU ISO C++ Library. This library is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.either version 3, or (at your option) any later version.";
    vector<string> wrap=str::word_wrap(txt, 60);
    dbg::printVector(wrap, "\n");
    cout<<endl<<"framed text:"<<endl<<endl;
    txt+="\n  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.";
    txt+="\n  Under Section 7 of GPL version 3, you are granted additional permissions described in the GCC Runtime Library Exception, version 3.1, as published by the Free Software Foundation.";
    cout<<str::frame("license",txt,80,5,'*',10);
    cout<<endl;
}
