#include <iostream>
#include <opencv2/opencv.hpp>
#include <Eigen/Eigen>

using namespace std;
using namespace cv;
using namespace Eigen;

class cameramodel
{
public:
    cameramodel(){}

    ~cameramodel(){}

    void polyProjection();

    void equaldistProjection();

    void planarProjection();

private:
    double ma[4];
    double mf;
    double mx0,my0;

    double mk[4];
    double mp[2];
    double mc[2];

    double mdownsample;


};