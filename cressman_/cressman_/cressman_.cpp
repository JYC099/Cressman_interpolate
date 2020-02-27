// cressman_.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime> 
#include <cmath>
#include "windows.h"
#include <string>
#include <fstream>
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
using namespace std;

class Point {
public:
    Point(double lat, double lon, double value) { _lat = lat; _lon = lon; _value = value; }
    double _lat;
    double _lon;
    double _value;
};

vector<double> cressman_interpolate(char *file_path, int grid_hres, double radius, double approximation,int &row) {//grid_hres网格分辨率单位是米，radius是搜索半径，approximation为近似值，row是输出结果每行栅格个数
    vector<Point> points_set;
    vector<double> out_interpolate;
    double min_lat = 999.0, max_lat = 0.0, min_lon = 999.0, max_lon = 0.0;

    GDALAllRegister();
    //读取shp文件
    GDALDataset* dataset = (GDALDataset*)GDALOpenEx(file_path, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (dataset == NULL)
    {
        cout<<"Open failed"<<endl;
        return out_interpolate;
    }
    OGRLayer* poLayer = dataset->GetLayer(0); //读取层
    OGRFeature* poFeature;
    poLayer->ResetReading();
    int i = 0;
    while ((poFeature = poLayer->GetNextFeature()) != NULL)
    {
        i = i++;
        OGRFeatureDefn* poFDefn = poFeature->GetDefnRef();
        OGRPoint* p = poFeature->GetGeometryRef()->toPoint();

        int iField;
        int n = poFDefn->GetFieldCount(); //获得字段的数目，不包括前两个字段（FID,Shape);
        for (iField = 0; iField < n; iField++) 
            points_set.push_back(Point(p->getY(), p->getX(), atof(poFeature->GetFieldAsString(iField))));//输入每个字段的值,Y是纬度
        
        
           
        OGRFeature::DestroyFeature(poFeature);
    }
    GDALClose(dataset);

  
    //获得点集合的最小外包矩形      
    for (int i = 0; i < points_set.size(); ++i) {
        if (points_set[i]._lat < min_lat)
            min_lat = points_set[i]._lat;
        if (points_set[i]._lon < min_lon)
            min_lon = points_set[i]._lon;
        if (points_set[i]._lat > max_lat)
            max_lat = points_set[i]._lat;
        if (points_set[i]._lon > max_lon)
            max_lon = points_set[i]._lon;
    }

    //采用WGS84坐标下经纬度长度，y是纵向有多少格网，x是横向有多少格网
    int y = (int)(((max_lat - min_lat) * 111000) / grid_hres), x = (int)(((max_lon - min_lon) * 85390) / grid_hres);
    row = x;

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
    char c[] = "D:/test_data/observe_points.shp";
    int r_num;
    vector<double> p = cressman_interpolate(c, 100, 8000, 0.1, r_num);//逼近值设为0.1

    int ScreenSizeX = r_num;
    int ScreenSizeY = p.size()/ r_num;

    ofstream fout("D:/OutputImage.ppm");
    fout << "P3\n" << ScreenSizeX << " " << ScreenSizeY << "\n255\n";
    cout << "开始渲染" << endl;
    for (int i = 0; i < p.size(); ++i) {

        if (p[i] < -300) {
            int ir = 0;
            int ig = 0;
            int ib = 0;
            fout << ir << " " << ig << " " << ib << "\n";
            continue;
        }
        if (p[i] > -300 && p[i] < 10) {
            int ir = 0;
            int ig = 0;
            int ib = 255;
            fout << ir << " " << ig << " " << ib << "\n";
            continue;
        }
        if (p[i] >= 10 && p[i] < 20) {
            int ir = 0;
            int ig = 250;
            int ib = 154;
            fout << ir << " " << ig << " " << ib << "\n";
            continue;
        }
        if (p[i] >= 20 && p[i] < 30) {
            int ir = 255;
            int ig = 193;
            int ib = 37;
            fout << ir << " " << ig << " " << ib << "\n";
            continue;
        }
        if (p[i] >= 30) {
            int ir = 255;
            int ig = 0;
            int ib = 0;
            fout << ir << " " << ig << " " << ib << "\n";
            continue;
        }      
    }

    cout << "渲染完毕" << endl;
    fout.close();
}


