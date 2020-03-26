#include <iostream>
#include <Utils.h>
using namespace std;
namespace xio=lyxutils::io;
namespace dbg=lyxutils::debug;
namespace str=lyxutils::str_utils;

int main(int argc,char **argv) {
    xio::CLParser parser;
    auto display_parse=[](xio::CLParser &parser){
        map<string,string> opts=parser.getAllOptions();
        vector<string> objs=parser.getParameters();
        cout<<"=>all options:"<<endl;
        for(auto it=opts.begin();it!=opts.end();++it){
            cout<<it->first<<":"<<it->second<<endl;
        }
        cout<<"=>all parameters:"<<endl;
        for_each(objs.begin(),objs.end(),[](string obj){cout<<obj<<"\t";});
        cout<<endl;
    };
    cout<<str::center("test without pattern",60,'=')<<endl;
    try{
        string cmd="ls --all -l --color=always test/ test2/";
        parser.parse(cmd);
        cout<<"the command is:"<<cmd<<endl;
        cout<<"has option `all`?:"<<parser.hasOption("all")<<endl;
        cout<<"color:"<<parser.getOptionArg("color")<<endl;
        display_parse(parser);
    }catch(invalid_argument ia){
        cout<<endl<<ia.what()<<endl;
    }
    cout<<str::center("test with pattern",60,'=')<<endl;
    xio::CLOption options[]={
            {"help",0,new char('h')},
            {"version",0,new char('v')},
            {"verbose",0,NULL},
            {"target",2,new char('t')},
            {"when",1,new char('w')},
            {"object",2,new char('o')},
            {"alignment",1,NULL}
    };
    parser.setOptionPattern(7,options);
    try{
        string cmd="g++ -h --verbose --target a.cpp -o liba.so --alignment=right --when folder1/ folder2/ folder3/";
        parser.parse(cmd);
        cout<<"correct test:"<<cmd<<endl;
        display_parse(parser);
        cmd="g++ --build-only folder1/";
        cout<<"unrecognized option test:"<<cmd<<endl;
        parser.parse(cmd);
        display_parse(parser);
    }catch(invalid_argument ia){
        cout<<endl<<ia.what()<<endl;
    }
    cout<<str::center("-",60,'-')<<endl;
    try{
        string cmd="g++ --target --verbose --alignment=right folder1/ folder2/";
        cout<<"lack required arguments test:"<<cmd<<endl;
        parser.parse(cmd);
        display_parse(parser);
    }catch(invalid_argument ia){
        cout<<endl<<ia.what()<<endl;
    }
    cout<<str::center("-",60,'-')<<endl;
    try{
        string cmd="g++ --verbose --alignment=right --target folder1/ folder2/";
        cout<<"required argument robs parameter test:"<<cmd<<endl;
        parser.parse(cmd);
        display_parse(parser);
    }catch(invalid_argument ia){
        cout<<endl<<ia.what()<<endl;
    }
    cout<<str::center("-",60,'-')<<endl;
    try{
        string cmd="g++ --verbose args -t a.cpp folder1/";
        cout<<"insert parameter after no-argument-option:"<<cmd<<endl;
        parser.parse(cmd);
        display_parse(parser);
    }catch(invalid_argument ia){
        cout<<endl<<ia.what()<<endl;
    }
    cout<<str::center("test options from file",60,'=')<<endl;
    try{
        parser.parse("work", "config.txt");
        display_parse(parser);
    }catch(invalid_argument ia){
        cout<<endl<<ia.what()<<endl;
    }catch(runtime_error re){
        cout<<endl<<re.what()<<endl;
    }
    cout<<str::center("free test",60,'=')<<endl;
    string line;
    cout<<"lyx@ubuntu:~$";
    while(getline(cin,line)){
        try{
            if(line=="exit")break;
            parser.parse(line);
            display_parse(parser);
        }
        catch(invalid_argument ia){
            cout<<endl<<ia.what()<<endl;
        }
        cout<<"lyx@ubuntu:~$";
    }
    return 0;
}