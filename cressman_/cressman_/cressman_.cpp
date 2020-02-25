// cressman_.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime> 
#include <cmath>
#include "windows.h"
using namespace std;

class Point {
public:
    Point(double lat, double lon, double value) { _lat = lat; _lon = lon; _value = value; }
    double _lat;
    double _lon;
    double _value;
};

vector<double> cressman_interpolate(double min_lat, double max_lat, double min_lon, double max_lon, int num,
    double min_value, double max_value, int grid_hres, double radius, double approximation) {//网格分辨率单位是米
    //生成随机点
    vector<Point> points_set;
    srand((unsigned)time(0));
    for (int i = 0; i < num; ++i) {
        points_set.push_back(Point(((max_lat - min_lat) * ((double)(rand() / (double)RAND_MAX)) + min_lat),
            ((max_lon - min_lon) * ((double)(rand() / (double)RAND_MAX)) + min_lon),
            ((max_value - min_value) * ((double)(rand() / (double)RAND_MAX)) + min_value)));
        //cout << points_set[i]._lat << " " << points_set[i]._lon << " " << points_set[i]._value << endl;
    }

    //采用WGS84坐标下经纬度长度，y是纵向有多少格网，x是横向有多少格网
    int y = (int)(((max_lat - min_lat) * 111000) / grid_hres), x = (int)(((max_lon - min_lon) * 85390) / grid_hres);
    vector<double> out_interpolate;

    //设置默认值
    for (int i = 0; i < y; ++i)
        for (int j = 0; j < x; ++j)
            out_interpolate.push_back(-999.9);

    //计算每一个栅格插值，默认栅格从左上开始计算，若半径内无观察点则保持默认值
    for (int i = 0; i < y; ++i) {
        for (int j = 0; j < x; ++j) {
            //判断有多少点在当前栅格点半径内
            vector<Point> in_points;
            vector<double> distance;
            vector<double> weight;
            for (auto k = points_set.begin(); k < points_set.end(); ++k) {
                double d_x = abs(k->_lon - (min_lon + (j - 0.5) * ((max_lon - min_lon) / x))) * 85390;
                double d_y = abs(k->_lat - (min_lat + (i - 0.5) * ((max_lat - min_lat) / y))) * 111000;
                double dist = sqrt(pow(d_x, 2) + pow(d_y, 2));
                //cout << dist << endl;
                if (dist < radius) {
                    in_points.push_back(*k);
                    distance.push_back(dist);
                }
            }
            if (in_points.size() == 0)
                continue;
            else {
                //用有效观察点的平均数作为第一猜测值
                double a =-999.9;//a为最终值
                double a0 = 0;//a0为猜测值
                double a_d = 0;//a_d是订正值
                for (auto p = in_points.begin(); p < in_points.end(); ++p)
                    a0 += p->_value;
                a0 /= in_points.size();
                
                while (1) {                  
                    double a_denominator = 0, a_numerator = 0;//代表计算a_d时的分子分母累计值
                   
                    for (int k = 0; k < in_points.size(); ++k) {
                        //计算权重
                        weight.push_back(((radius * radius) - (distance[k] * distance[k])) / ((radius * radius) + (distance[k] * distance[k])));
                        //cout << weight[k] << endl;
                        a_denominator += (weight[k] * weight[k] * (in_points[k]._value - a0));
                        a_numerator += weight[k];
                    }
                    //计算a_d
                    a_d = a_denominator / a_numerator;
                    a = a0 + a_d;
                    
                    //cout << a0 << endl;
                    if (abs(a_d) <= approximation)
                        break;
                    a0 = a;
                    
                }
                out_interpolate[i * x + j] = a;
            }
        }
    }


    return out_interpolate;
}



int main()
{
    vector<double> p=cressman_interpolate(114, 114.5, 30.3, 30.7, 1900, 20.0, 37.0, 1000, 2000, 0.1);//逼近值设为0.1
    for (int i = 0; i < p.size(); ++i)
        cout << p[i] << endl;
    cout << p.size();
}


