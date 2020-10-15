#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <Eigen/Eigen>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <ceres/ceres.h>
#include <ceres/rotation.h>


using namespace std;
using namespace cv;
using namespace Eigen;
using namespace ceres;

#define PI (3.1415926535897932346f)

struct polyProjecterror_in
{
    polyProjecterror_in(cv::Point2d p2D,cv::Point3d p3D,double f):mP2D(p2D),mP3D(p3D),mf(f){};
    
    //a 多项式系数a0,a2,a3,a4;     p0 像主点 x0,y0;        b 仿射变换参数b0,1,2 cde
    template<typename T>
    bool operator()(const T* const a,const T* const p0,const T* b,T * residual) const                 
    {
        // 2d project planar f
        T x=T(mP2D.x)-p0[0];
        T y=T(mP2D.y)-p0[1];
        
        Eigen::Matrix<T,2,2> F1;
        F1<<b[0],b[1],b[2],1;
        Eigen::Matrix<T,2,2> F;
        F=F1.inverse();
        T xx=x*F(0,0)+y*F(0,1);
        T yy=x*F(1,0)+y*F(1,1);
        T rho=sqrt(xx*xx+yy*yy);
        T rho2=rho*rho;
        T rho3=rho2*rho;
        T rho4=rho2*rho2;
        
        T Z=a[0]+a[1]*rho2+a[2]*rho3+a[3]*rho4;
        Z=-Z;
        T x_=mf*xx/Z;
        T y_=mf*yy/Z;

        // 3d project planar f
        T xp=mf*T(mP3D.x)/T(mP3D.z);
        T yp=mf*T(mP3D.y)/T(mP3D.z);
        T xxp=xp*b[0]+yp*b[1];
        T yyp=xp*b[2]+yp;

        residual[0]=x_-xxp;
        residual[1]=y_-yyp;
        return true;
    }

    static ceres::CostFunction* Create(const cv::Point2d p2d,const cv::Point3d p3d,double f)
    {
        // 2 残差维度    ---- parameters block  4 a参数维度 2 p0维度 3 b维度
        return (new ceres::AutoDiffCostFunction<polyProjecterror_in, 2, 4, 2,3>(
                 new polyProjecterror_in(p2d, p3d,f)));                             
    }

private:
    const cv::Point2d mP2D;             // 2d point
    const cv::Point3d mP3D;
    const double mf;

};


//定义鱼眼成像模型 多项式模型  内参、外参  quaternion
struct polyProjecterror_exin
{
    polyProjecterror_exin(cv::Point2d p2D,cv::Point3d p3D,double f):mP2D(p2D),mP3D(p3D),mf(f){};
    
    //内参：a 多项式系数a0,a2,a3,a4;     p0 像主点 x0,y0;        b 仿射变换参数b0,1,2 cde  
    //外参：quaternion 四元数               t 平移矩阵
    template<typename T>
    bool operator()(const T* const a,const T* const p0,const T* b, const T* quaternion, const T* t,T * residual) const                 
    {
        // transform 3d point from world coordinate to camera coordinate 
        T p1[3];
        T p2[3];
        p1[0]=T(mP3D.x);
        p1[1]=T(mP3D.y);
        p1[2]=T(mP3D.z);
        ceres::QuaternionRotatePoint(quaternion,p1,p2);
        p2[0] += t[0]; p2[1] += t[1]; p2[2] += t[2];

        // 2d project planar f
        T x=T(mP2D.x)-p0[0];
        T y=T(mP2D.y)-p0[1];

        T temp=x;
        x=y;y=temp;

        Eigen::Matrix<T,2,2> F;
        F<<b[0],b[1],b[2],T(1);
        F.inverse();
        T xx=x*F(0,0)+y*F(0,1);
        T yy=x*F(1,0)+y*F(1,1);
        T rho=sqrt(xx*xx+yy*yy);
        T rho2=rho*rho;
        T rho3=rho2*rho;
        T rho4=rho2*rho2;
        
        T Z=a[0]+a[1]*rho2+a[2]*rho3+a[3]*rho4;
        Z=-Z;

        T x_=T(mf)*xx/Z;
        T y_=T(mf)*yy/Z;

        temp=x_;
        x_=y_;y_=temp;
        x_=x_+p0[0];
        y_=y_+p0[1];

        // 3d project planar f
        
        T xp=T(mf)*p2[0]/p2[2]+p0[0];
        T yp=T(mf)*p2[1]/p2[2]+p0[1];

        residual[0]=x_-xp;
        residual[1]=y_-yp;
        return true;
    }

    static ceres::CostFunction* Create(const cv::Point2d p2d,const cv::Point3d p3d,double f)
    {
        // 2 残差维度    ---- parameters block  4 a参数维度 ,2 p0维度 ,3 b维度 , 4 四元数 ,3 平移

        return (new ceres::AutoDiffCostFunction<polyProjecterror_exin, 2, 4, 2, 3, 4, 3>(
                 new polyProjecterror_exin(p2d, p3d,f)));   
        // return (new ceres::NumericDiffCostFunction<polyProjecterror_exin, ceres::RIDDERS,2, 4, 2,3,4,3>(
        //          new polyProjecterror_exin(p2d, p3d,f)));                          
    }

private:
    const cv::Point2d mP2D;             // 2d point
    const cv::Point3d mP3D;
    const double mf;

};

// // Rodrigues
// struct polyProjecterror_exin 
// {
//     polyProjecterror_exin(cv::Point2d p2D,cv::Point3d p3D,double f):mP2D(p2D),mP3D(p3D),mf(f){};
    
//     //内参：a 多项式系数a0,a2,a3,a4;     p0 像主点 x0,y0;        b 仿射变换参数b0,1,2 cde  
//     //外参：rvec Rodrigues              t 平移矩阵
//     template<typename T>
//     bool operator()(const T* const a,const T* const p0,const T* b, const T* rvec, const T* t,T * residual) const                 
//     {
//         // 2d project planar f
//         T x=T(mP2D.x)-p0[0];
//         T y=T(mP2D.y)-p0[1];
//         T xx=x*b[0]+y*b[1];
//         T yy=x*b[2]+y;
//         T rho=sqrt(xx*xx+yy*yy);
//         T rho2=rho*rho;
//         T rho3=rho2*rho;
//         T rho4=rho2*rho2;
        
//         T Z=a[0]+a[1]*rho2+a[2]*rho3+a[3]*rho4;
//         Z=-Z;
//         T x_=mf*xx/Z+p0[0];
//         T y_=mf*yy/Z+p0[1];
        

//         // 3d project planar f

//         // transform 3d point from world coordinate to camera coordinate 
//         T p1[3];
//         T p2[3];
//         p1[0]=T(mP3D.x);
//         p1[1]=T(mP3D.y);
//         p1[2]=T(mP3D.z);
//         ceres::AngleAxisRotatePoint(rvec, p1, p2);
//         p2[0] =p2[0]+t[0]; p2[1] =p2[1]+ t[1]; p2[2] =p2[2]+ t[2];

//         T xp=mf*p2[0]/p2[2]+p0[0];
//         T yp=mf*p2[1]/p2[2]+p0[1];

//         residual[0]=x_-xp;
//         residual[1]=y_-yp;
//         return true;
//     }

//     static ceres::CostFunction* Create(const cv::Point2d p2d,const cv::Point3d p3d,double f)
//     {
//         // 2 残差维度    ---- parameters block  4 a参数维度 ,2 p0维度 ,3 b维度 , 3 Rodrigues ,3 平移

//         return (new ceres::AutoDiffCostFunction<polyProjecterror_exin, 2, 4, 2,3,3,3>(
//                  new polyProjecterror_exin(p2d, p3d,f)));   
//         // return (new ceres::NumericDiffCostFunction<polyProjecterror_exin, ceres::RIDDERS,2, 4, 2,3,4,3>(
//         //          new polyProjecterror_exin(p2d, p3d,f)));                          
//     }

// private:
//     const cv::Point2d mP2D;             // 2d point
//     const cv::Point3d mP3D;
//     const double mf;

// };

// struct polyProjecterror_exin
// {
//     polyProjecterror_exin(cv::Point2d p2D,cv::Point3d p3D,double f):mP2D(p2D),mP3D(p3D),mf(f){};
    
//     //内参：a 多项式系数a0,a2,a3,a4;     p0 像主点 x0,y0;        b 仿射变换参数b0,1,2 cde  
//     //外参：rvec Rodrigues              t 平移矩阵
//     template<typename T>
//     bool operator()(const T* const a,const T* const p0,const T* b, const T* rvec, const T* t,T * residual) const                 
//     {
//         // 2d project planar f
//         T x=T(mP2D.x)-p0[0];
//         T y=T(mP2D.y)-p0[1];
//         T xx=x*b[0]+y*b[1];
//         T yy=x*b[2]+y;
//         T rho=sqrt(xx*xx+yy*yy);
//         T rho2=rho*rho;
//         T rho3=rho2*rho;
//         T rho4=rho2*rho2;
        
//         T Z=a[0]+a[1]*rho2+a[2]*rho3+a[3]*rho4;
//         Z=-Z;
//         // T x_=mf*xx/Z+p0[0];
//         // T y_=mf*yy/Z+p0[1];
//         T r1=sqrt(xx*xx+yy*yy+Z*Z);
//         T x1_n=xx/r1;
//         T y1_n=yy/r1;
//         T z1_n=Z/r1;
        

//         // 3d project planar f

//         // transform 3d point from world coordinate to camera coordinate 
//         T p1[3];
//         T p2[3];
//         p1[0]=T(mP3D.x);
//         p1[1]=T(mP3D.y);
//         p1[2]=T(mP3D.z);
//         ceres::AngleAxisRotatePoint(rvec, p1, p2);
//         p2[0] =p2[0]+t[0]; p2[1] =p2[1]+ t[1]; p2[2] =p2[2]+ t[2];

//         // T xp=mf*p2[0]/p2[2]+p0[0];
//         // T yp=mf*p2[1]/p2[2]+p0[1];
//         T r2=sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[1]*p2[1]);
//         T x2_n=p2[0]/r2;
//         T y2_n=p2[1]/r2;
//         T z2_n=p2[2]/r2;

//         residual[0]=x1_n-x2_n;
//         residual[1]=y1_n-y2_n;
//         residual[2]=z1_n-z2_n;
//         return true;
//     }

//     static ceres::CostFunction* Create(const cv::Point2d p2d,const cv::Point3d p3d,double f)
//     {
//         // 3 残差维度    ---- parameters block  4 a参数维度 ,2 p0维度 ,3 b维度 , 3 Rodrigues ,3 平移

//         return (new ceres::AutoDiffCostFunction<polyProjecterror_exin, 3, 4, 2,3,3,3>(
//                  new polyProjecterror_exin(p2d, p3d,f)));   
//         // return (new ceres::NumericDiffCostFunction<polyProjecterror_exin, ceres::RIDDERS,2, 4, 2,3,4,3>(
//         //          new polyProjecterror_exin(p2d, p3d,f)));                          
//     }

// private:
//     const cv::Point2d mP2D;             // 2d point
//     const cv::Point3d mP3D;
//     const double mf;

// };


//定义鱼眼成像模型 等距投影模型  内参
struct equaldistprojecterror_in
{
    equaldistprojecterror_in(cv::Point2d p2D,cv::Point3d p3D,double f):mP2D(p2D),mP3D(p3D),mf(f){};

     template<typename T>
    bool operator()(const T* const k,const T* const p,const T* c,const T* const p0,T * residual) const                 
    {

        
        // fisheye -> rectify planar
        T x=T(mP2D.x)-p0[0];
        T y=T(mP2D.y)-p0[1];

        T r=sqrt(x*x+y*y);
        T xr=mf*x/r*tan(r/mf);
        T yr=mf*y/r*tan(r/mf);
        
        T r2=xr*xr+yr*yr;
        T dx=xr*(k[0]*r2+k[1]*r2*r2+k[2]*r2*r2*r2+k[3]*r2*r2*r2)+T(2)*p[0]*xr*yr+p[1]*(r2+T(2)*xr*xr)+c[0]*xr+c[1]*yr;
        T dy=yr*(k[0]*r2+k[1]*r2*r2+k[2]*r2*r2*r2+k[3]*r2*r2*r2)+p[0]*(r2+T(2)*yr*yr)+T(2)*p[0]*yr*xr+c[0]*xr+c[1]*yr;
        T x_=xr+dx;
        T y_=yr+dy;

        // 3D -> planar
        T xp=mf*T(mP3D.x)/T(mP3D.z);
        T yp=mf*T(mP3D.y)/T(mP3D.z);

        residual[0]=x_-xp;
        residual[1]=y_-yp;
        return true;
    }

    static ceres::CostFunction* Create(const cv::Point2d p2d,const cv::Point3d p3d,double f)
    {
        // 2 残差维度    ---- parameters block  4 k参数维度 2 p维度 2 c维度 2 p0像主点 
        return (new ceres::AutoDiffCostFunction<equaldistprojecterror_in, 2, 4, 2, 2, 2>(
                 new equaldistprojecterror_in(p2d, p3d,f)));                             
    }



    private:
    const cv::Point2d mP2D;             // 2d point
    const cv::Point3d mP3D;
    const double mf;
    
};

//定义鱼眼成像模型 等距投影模型  内参外参
struct equaldistprojecterror_exin  
{
    equaldistprojecterror_exin(cv::Point2d p2D,cv::Point3d p3D,double f):mP2D(p2D),mP3D(p3D),mf(f){};

     template<typename T>
    bool operator()(const T* const k,const T* const p,const T* c,const T* const p0,const T* quaternion, const T* t,T * residual) const                 
    {
        // transform 3d point from world coordinate to camera coordinate 
        T p1[3];
        T p2[3];
        p1[0]=T(mP3D.x);
        p1[1]=T(mP3D.y);
        p1[2]=T(mP3D.z);
        ceres::QuaternionRotatePoint(quaternion,p1,p2);
        p2[0] += t[0]; p2[1] += t[1]; p2[2] += t[2];

        // fisheye -> rectify planar

        T x=T(mP2D.x)-p0[0];
        T y=T(mP2D.y)-p0[1];

        T r=sqrt(x*x+y*y);
        T xr=mf*x/r*tan(r/mf);
        T yr=mf*y/r*tan(r/mf);
        
        T r2=xr*xr+yr*yr;
        T dx=xr*(k[0]*r2+k[1]*r2*r2+k[2]*r2*r2*r2+k[3]*r2*r2*r2)+T(2)*p[0]*xr*yr+p[1]*(r2+T(2)*xr*xr)+c[0]*xr+c[1]*yr;
        T dy=yr*(k[0]*r2+k[1]*r2*r2+k[2]*r2*r2*r2+k[3]*r2*r2*r2)+p[0]*(r2+T(2)*yr*yr)+T(2)*p[0]*yr*xr+c[0]*xr+c[1]*yr;
        T x_=xr+dx+p0[0];
        T y_=yr+dy+p0[1];

        // 3D -> planar
        T xp=mf*p2[0]/p2[2]+p0[0];
        T yp=mf*p2[1]/p2[2]+p0[1];

        residual[0]=x_-xp;
        residual[1]=y_-yp;
        return true;
    }

    static ceres::CostFunction* Create(const cv::Point2d p2d,const cv::Point3d p3d,double f)
    {
        // 2 残差维度    ---- parameters block  4 k参数维度 2 p维度 2 c维度 2 p0像主点 4 四元数 ,3 平移
        return (new ceres::AutoDiffCostFunction<equaldistprojecterror_exin, 2, 4, 2, 2, 2,4,3>(
                 new equaldistprojecterror_exin(p2d, p3d,f)));                             
    }



    private:
    const cv::Point2d mP2D;             // 2d point
    const cv::Point3d mP3D;
    const double mf;
    
};


class exincalib
{
    public:
        exincalib(double &x0,double &y0,double &f,cv::Mat image)
        {
            minip0[0]=x0;
            minip0[1]=y0;
            mf=f;
            mimage=image;
        }

        ~exincalib(){}

        //读取txt 对应点
        void readmatchestxtfile(string &path);

        // ceres 优化 内参、外参初始化  poly model
        void caliboptimizeiniforpoly(vector<cv::Point2d> &vp2d,vector<cv::Point3d> &vp3d);
        // ceres 优化 内参、外参初始化  equal distance model
        void caliboptimizeiniforequaldist(vector<cv::Point2d> &vp2d,vector<cv::Point3d> &vp3d);

        // ceres 优化 内参、外参初始化  poly model 赋值
        void readexininital(double a[4],double p0[2],double b[3],double quaternion[4],double t[3]);
        void caliboptimizeiniforpoly();

        // ransac 初始化
        void ransacinitial(int iterationnum,int setnum);
        
        // poly model 优化求解
        void optimizesolverforpoly(vector<cv::Point2d> &vp2d,vector<cv::Point3d> &vp3d);
        // equal distance model 优化求解
        void optimizesolverforequaldist(vector<cv::Point2d> &vp2d,vector<cv::Point3d> &vp3d);
        //ransac 优化
        void ransac_solver(int &model);
       
        // poly model  检测内点，用于最后优化以及定权
        int checkinlinersforpoly( vector<bool> &vpinliners,double &error,int th=1);
        // equal distance model  检测内点，用于最后优化以及定权
        int checkinlinersforequaldist(vector<bool> &vpinliners,double &error,int th=1);


        // //Projection: poly / equal-distance / planar 
        // // 相机坐标系中的3D点
        // void Fisheye_C_PolyProjetion(cv::Point2d &p2d,cv::Point3d &p3d);             // poly model  fisheye->3D

        // void Fisheye_C_EqualdistProjection(cv::Point2d &p2d,cv::Point3d &p3d);       // Equal distance model  fisheye->3D

        // void C_2D_PlanarProjection(cv::Point2d &p2d,cv::Point3d &p3d);               // 3D -> 2D

        void CalculateRotationMatFromEulerAngle(double Rx, double Ry, double Rz, double * R); //Rz*Ry*Rx

        void CalculateEulerAngleFromRotationMat(double &Rx, double &Ry, double &Rz, double * R);


        void savereasults(string path);

    private:

        cv::Mat mimage;

        int mnum;                               //匹配点数
        vector<cv::Point2d> mvP2D;             // 2d point
        vector<cv::Point3d> mvP3D;
        vector<string> mvname;
        ceres::Problem exincalibproblem;

        int mmodel;                            // model=0: poly   model =1 :equal distance

        // 内外参初始值
        double inia[4];
        double inib[3];
        double inip0[2];

        double iniq[4];
        double init[3];

        // 内参
        double mf;
        double minip0[2]; 
        double mp0[2];                         //像主点

        //ploy model
        double ma[4];                          //多项式系数    初始化值
        double mb[3];                          //仿射变换系数

        //equal distance model
        double mk[4];
        double mp[2];
        double mc[2];

        // 外参
        cv::Mat mMR;
        cv::Mat mMt;
        cv::Mat mMT;

        double mquaternion[4];
        double mt[3];
        double mrvec[3];

        double mR[9];


        //Ransac
        vector<vector<int>> rvvnum;                // ransac num of every set 
        int miterationnum;
        int msetnum;
        vector<bool> mvbBestInliers;

        int minliner;                              // number of inliners
        double error;
};