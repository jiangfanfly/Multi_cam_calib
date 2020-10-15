#include <iostream>
#include <opencv2/opencv.hpp>
#include <Eigen/Eigen>
#include <dirent.h>

#include "exincalib.h"

using namespace std;
using namespace cv;
using namespace Eigen;

bool cmp(std::string const &arg_a, std::string const &arg_b)
{
    return arg_a.size() < arg_b.size() || (arg_a.size() == arg_b.size() && arg_a < arg_b);
}

vector<string> getFiles(string cate_dir)
{
    vector<string> files;//存放文件名

    DIR *dir;
    struct dirent *ptr;
    char base[1000];

    if ((dir=opendir(cate_dir.c_str())) == NULL)
        {
        perror("Open dir error...");
                exit(1);
        }

    while ((ptr=readdir(dir)) != NULL)
    {
        if(strcmp(ptr->d_name,".")==0 || strcmp(ptr->d_name,"..")==0)    ///current dir OR parrent dir
                continue;
        else if(ptr->d_type == 8)    ///file
            //printf("d_name:%s/%s\n",basePath,ptr->d_name);
            files.push_back(ptr->d_name);
        else if(ptr->d_type == 10)    ///link file
            //printf("d_name:%s/%s\n",basePath,ptr->d_name);
            continue;
        else if(ptr->d_type == 4)    ///dir
        {
            files.push_back(ptr->d_name);
            /*
                memset(base,'\0',sizeof(base));
                strcpy(base,basePath);
                strcat(base,"/");
                strcat(base,ptr->d_nSame);
                readFileList(base);
            */
        }
    }
    closedir(dir);


    //排序，按从小到大排序
    sort(files.begin(), files.end(),cmp);
    return files;
}

// test useless
void readmatches(string &path)
{
            int mnum;                               //匹配点数
        vector<cv::Point2d> mvP2D;             // 2d point
        vector<cv::Point3d> mvP3D;
        vector<string> mvname;
    ifstream of1;
    of1.open(path.data());
    int i=1;
    cv::Point2d p2d;
    cv::Point3d p3d;
    string name;

    while(!of1.eof())
    {
        if(i == 1)
        {
            of1>>mnum>>mnum;
            // mvP2D.resize(mnum);
            // mvP3D.resize(mnum);
            i++;
        }
        else
        {
            of1>>name>>p2d.x>>p2d.y>>p3d.x>>p3d.y>>p3d.z;
            mvP2D.push_back(p2d);
            mvP3D.push_back(p3d);
            mvname.push_back(name);
        }
    }
    mnum=mvP2D.size();
    cout<<"读入 "<<mvname.size()<<"组数据！"<<endl;
    
    vector<cv::Point2d> vp2d;
    cv::Point2d p;
    for(int i=0;i<mvP2D.size();i++)
    {
        // fisheye -> rectify planar
        double x=mvP2D[i].x-1728;
        double y=mvP2D[i].y-1728;

        double r=sqrt(x*x+y*y);
        double xr=1300*x/r*tan(r/1300);
        double yr=1300*y/r*tan(r/1300);
        
        double r2=xr*xr+yr*yr;
        double x_=xr+1728;
        double y_=yr+1728;
        p.x=x_;
        p.y=y_;
        vp2d.push_back(p);

    }
    cv::Mat K;
    //K = cv::Mat::zeros(3, 3, CV_32FC1);
    K = cv::Mat::zeros(3, 3, CV_64FC1);
	K.at<double>(0, 0) = 1300;
	K.at<double>(0, 2) = 1728;
	K.at<double>(1, 1) = 1300;
	K.at<double>(1, 2) = 1728;
	K.at<double>(2, 2) = 1.0;
    cv::Mat r = cv::Mat::zeros(3, 1, CV_64FC1);
	cv::Mat t = cv::Mat::zeros(3, 1, CV_64FC1);
    cv::Mat R,in;
    //cv::solvePnP(mvP3D,mvP2D,K,Mat(),r,t,false,SOLVEPNP_EPNP);
    cv::solvePnPRansac(mvP3D,vp2d,K,Mat(),r,t,false,100,10,0.9, in,SOLVEPNP_EPNP);
    cv::Rodrigues(r,R);
    cout<<"R:"<<R<<endl<<"t:"<<t<<endl;

}

int main()
{
    // string pathimageA="/home/jiangfan/data/20190326/20190326/ONEPKRMNS/IMG_20190326_103739_260/A.tif";
    // //string pathimageB="/home/jiangfan/data/20190326/20190326/ONEPKRMNS/IMG_20190326_103739_260/B.tif";

    // string pathtxtfileA="/home/jiangfan/data/20190326/20190326/fisheye/ONE2/IMG_20190326_103739_260/A.ctrl";
    // //string pathtxtfileB="/home/jiangfan/data/20190326/20190326/fisheye/ONE2/IMG_20190326_103739_260/B.ctrl";

    // string path_totxt="/home/jiangfan/data/20190326/20190326/fisheye/ONE/";
    
    // vector<string> ctrlpath;
    // ctrlpath=getFiles(path_totxt);

    // cv::Mat img=cv::imread(pathimageA);
    // double x0=img.cols/2;
    // double y0=img.rows/2;
    // double f=1300;

    // string path1=path_totxt+ctrlpath[0]+"/B.ctrl";
    // readmatches(path1);

    // exincalib calib(x0,y0,f,img);

    // // model=0: poly   model =1 :equal distance
    // int model=0; 

    // // for(int i=0;i<ctrlpath.size();i++)
    // for(int i=0;i<1;i++)
    // {

    //     string path=path_totxt+ctrlpath[i]+"/B.ctrl";
    //     calib.readmatchestxtfile(path);

    //     //calib.caliboptimizeini();

    //     calib.ransacinitial(10,500);
    //     calib.ransac_solver(model);
    // }

    string path="/home/jiangfan/桌面/5_14/G1_reasults/51/p3d2d1_5_d.txt";
    string imgpath="/home/jiangfan/桌面/5_14/G1/1-2/Teche0001_rect.jpg";
    string exinpath="/home/jiangfan/桌面/5_14/exin_optimize/51/5.txt";
    string reasultspath="/home/jiangfan/桌面/5_14/G1_reasults/inex_op_1.txt";
    cv::Mat img=imread(imgpath);
    double x0=img.cols/2;
    double y0=img.rows/2;
    double f=1300;
    //double f=396.006564204808;
    Eigen::Matrix2d F,F1;
    F<<1.004451,0.000401,-0.000676,1;
    F1=F.inverse();
    cout<<"F1:"<<F1<<endl;
    cout<<"F1:"<<F1(0,1)<<" "<<F1(1,0)<<endl;

    // read exin para
    double inia[4];
    double inib[3];
    double inip0[2];
    double iniq[4];
    double init[3];
    double rr[9];
    double a;
    double h,w;

    ifstream of;
    of.open(exinpath.data());
    of>>a>>inia[0]>>a>>inia[1]>>inia[2]>>inia[3];
    of>>inip0[1]>>inip0[0];
    of>>inib[0]>>inib[1]>>inib[2];
    of>>h>>w;
    of>>rr[0]>>rr[1]>>rr[2]>>rr[3]>>rr[4]>>rr[5]>>rr[6]>>rr[7]>>rr[8];
    of>>init[0]>>init[1]>>init[2];
    of.close();
    Eigen:Matrix3d R;
    R<<rr[0],rr[1],rr[2],
        rr[3],rr[4],rr[5],
        rr[6],rr[7],rr[8];
    cout<<"R"<<R<<endl;
    Eigen::Quaterniond q=Eigen::Quaterniond(R);
    q.normalized();
    iniq[0]=q.w();iniq[1]= q.x();iniq[2]= q.y();iniq[3]= q.z();
    cout<<"q"<<iniq[0]<<" "<<iniq[1]<<" "<<iniq[2]<<" "<<iniq[3]<<" "<<endl;

    exincalib calib(x0,y0,f,img);
    calib.readmatchestxtfile(path);
    calib.readexininital(inia,inip0,inib,iniq,init);
    calib.ransacinitial(10,80);
    int model=0; 
    calib.ransac_solver(model);
    calib.savereasults(reasultspath);


    return 0;
}