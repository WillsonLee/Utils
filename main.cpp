#include <iostream>
#include <sstream>
#include "Utils.h"
#include <algorithm>
#include <mutex>
//#include <Windows.h>
using namespace std;
using namespace lyxutils::colors;
#define ual lyxutils::algorithms
class Pnt{
public:
    double x;
    double y;
    double z;
    Pnt() :x(0), y(0), z(0) {}
    Pnt(float x_, float y_, float z_) :x(x_), y(y_), z(z_) {}
};

int main()
{
    std::vector<std::vector<double> > a,b;
    for(int i=0;i<3;++i)a.push_back(std::vector<double>());
    a[0].push_back(2);
    a[0].push_back(3);
    a[0].push_back(5);
    a[0].push_back(7);
    a[0].push_back(9);
    a[1].push_back(1);
    a[1].push_back(4);
    a[1].push_back(5);
    a[1].push_back(8);
    a[1].push_back(9);
    a[2].push_back(1);
    a[2].push_back(1);
    a[2].push_back(3);
    a[2].push_back(6);
    a[2].push_back(7);
    lin::cov_samples(a, b);
    cout << "cov=" << endl;
    np::print(b);
    std::vector<Pnt> pList;
    pList.push_back(Pnt(2, 1, 1));
    pList.push_back(Pnt(3, 4, 1));
    pList.push_back(Pnt(5, 5, 3));
    pList.push_back(Pnt(7, 8, 6));
    pList.push_back(Pnt(9, 9, 7));
    Eigen::Matrix3d mat;
    lin::cov_points(pList, mat);
    cout << "mat=" << endl << mat << endl;
    std::vector<int> x_idx, y_idx;
    np::where(b, x_idx, y_idx, [](double a)->bool {return a > 35; });
    cout << "indices that marks elements greater than 35:" << endl;
    np::print(x_idx);
    np::print(y_idx);
    std::vector<float> t1;
    np::range<float>(t1, 1, 15, 2);
    np::where(t1, x_idx, [](float a)->bool {return a > 10; });
    cout << "test one d vector search..." << endl << "t1=" << endl;
    np::print(t1);
    cout << "where:" << endl;
    np::print(x_idx);
    np::apply_where(b, [](double a)->bool {return a > 35; }, [](double &a) {a = 35; });
    np::apply_where(t1, [](float a)->bool {return a > 10; }, [](float &a) {a = 10; });
    cout << "apply_where test:" << endl;
    np::print(b);
    np::print(t1);
//    system("pause");
    return 0;
}