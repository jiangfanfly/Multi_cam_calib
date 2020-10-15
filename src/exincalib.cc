#include "exincalib.h"
#define PI (3.1415926535897932346f)
void exincalib::readmatchestxtfile(string &path)
{
    ifstream of1;
    of1.open(path.data());
    int i=1;
    cv::Point2d p2d;
    cv::Point3d p3d;
    string name;

    while(!of1.eof())
    {
        // if(i == 1)
        // {
        //     of1>>mnum>>mnum;
        //     // mvP2D.resize(mnum);
        //     // mvP3D.resize(mnum);
        //     i++;
        // }
        // else
        {
            //of1>>name>>p2d.x>>p2d.y>>p3d.x>>p3d.y>>p3d.z;
            of1>>name>>p3d.x>>p3d.y>>p3d.z>>p2d.x>>p2d.y;
            //p3d.z=-p3d.z;
            mvP2D.push_back(p2d);
            mvP3D.push_back(p3d);
            mvname.push_back(name);
        }
    }
    mnum=mvP2D.size();
    cout<<"读入 "<<mvname.size()<<"组数据！"<<endl;
}

void exincalib::caliboptimizeiniforpoly(vector<cv::Point2d> &vp2d,vector<cv::Point3d> &vp3d)
{

    // 外参初始化

    //R=[0,1,0,        t=[0,-41.7024893356828,0]         近似初始化值   (-pi/2,-pi/2,0,'xyz')
    //   0, 0,-1, 
    //   -1 ,0,0]

    
    cv::Mat cameraMatrix=cv::Mat::zeros(3, 3, CV_64FC1);
    cameraMatrix.at<double>(0, 0) = mf;
	cameraMatrix.at<double>(0, 2) = minip0[0];
	cameraMatrix.at<double>(1, 1) = mf;
	cameraMatrix.at<double>(1, 2) = minip0[1];
	cameraMatrix.at<double>(2, 2) = 1.0;
    cv::Mat r = cv::Mat::zeros(3, 1, CV_64FC1);
	cv::Mat t = cv::Mat::zeros(3, 1, CV_64FC1);
    cv::Mat in;
    cv::solvePnPRansac(vp3d,vp2d,cameraMatrix,Mat(),r,t,true,50,10,0.9, in,SOLVEPNP_EPNP);
    cv::Rodrigues(r,mMR);
    mMt=t.clone();
    

    Eigen::Matrix3d eR;
    cv::cv2eigen(mMR,eR);
    
    // mR[0]=eR(0,0);mR[1]=eR(0,1);mR[2]=eR(0,2);
    // mR[3]=eR(1,0);mR[4]=eR(1,1);mR[5]=eR(1,2);
    // mR[6]=eR(2,0);mR[7]=eR(2,1);mR[8]=eR(2,2);
    //cout<<"eR:"<<eR<<endl;

    eR<<0.809330,0.587295,0.008355,
    0.024946,-0.020158,-0.999486,
    -0.586824,0.809122,-0.030965;
    Eigen::AngleAxisd V2;
    V2.fromRotationMatrix(eR);
    Eigen::Vector3d axis;
    axis=V2.axis();
    double angle;
    angle=V2.angle();
    axis=angle*axis;
    
    mrvec[0]=axis[0];mrvec[1]=axis[1];mrvec[2]=axis[2];

    Eigen::Quaterniond q=Eigen::Quaterniond(eR);
    q.normalize();
    // mrvec[0]=r.at<double>(0,0);
    // mrvec[1]=r.at<double>(1,0);
    // mrvec[2]=r.at<double>(2,0);
    mquaternion[0]=q.w();mquaternion[1]= q.x();mquaternion[2]= q.y();mquaternion[3]= q.z();
    mt[0]=mMt.at<double>(0,0);mt[1]=mMt.at<double>(0,1);mt[2]=mMt.at<double>(0,2);

    mt[0]=-4014.163818 ;mt[1]=-140.713547;mt[2]=-895.277405;

    cout<<"R:"<<eR<<endl<<"t:"<<mt<<endl;


    // 内参初始化
    Eigen::MatrixXd A(2*vp2d.size(),4);
    Eigen::MatrixXd L(2*vp2d.size(),1);
    for(size_t i=0;i<vp2d.size();i++)
    {
        double xi=vp2d[i].x-minip0[0];
        double yi=vp2d[i].y-minip0[1];              // 像素坐标转换到图像坐标（鱼眼）

        double p1[3];                               // 世界坐标系转换到相机坐标系
        double p2[3];
        p1[0]=vp3d[i].x;
        p1[1]=vp3d[i].y;
        p1[2]=vp3d[i].z;
        //ceres::QuaternionRotatePoint(mquaternion,p1,p2);
        ceres::AngleAxisRotatePoint(mrvec, p1, p2);

        p2[0] += mt[0]; p2[1] += mt[1]; p2[2] += mt[2];   

        double xri=p2[0]/p2[2]*mf;
        double yri=p2[1]/p2[2]*mf;     //相机坐标转换到图像坐标（纠正）

        double rho=sqrt(xi*xi+yi*yi);
        double rho2=rho*rho;
        double rho3=rho*rho2;
        double rho4=rho2*rho2;

        A(2*i,0)=xri;A(2*i,1)=xri*rho2;A(2*i,2)=xri*rho3;A(2*i,3)=xri*rho4;            // xr/f=xi/f(rhp) => xr*f(rhp)=xi*f
        A(2*i+1,0)=yri;A(2*i+1,1)=yri*rho2;A(2*i+1,2)=yri*rho3;A(2*i+1,3)=yri*rho4;

        L(2*i,0)=-mf*xi;
        L(2*i+1,0)=-mf*yi;
    }

    Eigen::MatrixXd ATA=A.transpose()*A;
    Eigen::MatrixXd ATL=A.transpose()*L;

    Eigen::MatrixXd a;
    a=ATA.inverse()*ATL;

    mb[0]=1.0;mb[1]=0.0;  mb[2]=0.0;
    ma[0]=a(0,0);ma[1]=a(1,0);ma[2]=a(2,0);ma[3]=a(3,0);
    mp0[0]=minip0[0];mp0[1]=minip0[1];
    vector<bool> b;
   

    ma[0]=-1.959116e+03;ma[1]=2.188031e-04;ma[2]=-4.770103e-08;ma[3]=2.324483e-11;
    mb[0]=1.003344;mb[1] =0.000988;mb[2]=-0.000726;
    mp0[0]=1524.465110;mp0[1]=2031.354394 ;

     //checkinlinersforpoly(b);
    cout<<"initialize successful!"<<endl;
    cout<<"before solver:"<<endl<<"a: "<<ma[0]<<" "<<ma[1]<<" "<<ma[2]<<" "<<ma[3]<<" "<<endl;
    cout<<"p0: "<<mp0[0]<<" "<<mp0[1]<<endl;
    cout<<"b: "<<mb[0]<<" "<<mb[1]<<" "<<mb[2]<<endl;

}

void exincalib::readexininital(double a[4],double p0[2],double b[3],double quaternion[4],double t[3])
{
    inia[0]=a[0];inia[1]=a[1];inia[2]=a[2];inia[3]=a[3];
    inip0[0]=p0[0];inip0[1]=p0[1];
    inib[0]=b[0];inib[1]=b[1];inib[2]=b[2];

    iniq[0]=quaternion[0];iniq[1]=quaternion[1];iniq[2]=quaternion[2];iniq[3]=quaternion[3];
    init[0]=t[0];init[1]=t[1];init[2]=t[2];
}

void exincalib::caliboptimizeiniforpoly()
{
    ma[0]=inia[0];ma[1]=inia[1];ma[2]=inia[2];ma[3]=inia[3];
    mp0[0]=inip0[0];mp0[1]=inip0[1];
    mb[0]=inib[0];mb[1]=inib[1];mb[2]=inib[2];

    mquaternion[0]=iniq[0];mquaternion[1]=iniq[1];mquaternion[2]=iniq[2];mquaternion[3]=iniq[3];
    mt[0]=init[0];mt[1]=init[1];mt[2]=init[2];

    cout<<"initialize successful!"<<endl;
    cout<<"before solver:"<<endl<<"a: "<<ma[0]<<" "<<ma[1]<<" "<<ma[2]<<" "<<ma[3]<<" "<<endl;
    cout<<"p0: "<<mp0[0]<<" "<<mp0[1]<<endl;
    cout<<"b: "<<mb[0]<<" "<<mb[1]<<" "<<mb[2]<<endl;
    //cout<<"R:"<<eR<<endl<<"t:"<<mt<<endl;
    // vector<bool> vpinliners;
    // int in;
    // double error;
    // in=checkinlinersforpoly(vpinliners,error,5);
}

void exincalib::caliboptimizeiniforequaldist(vector<cv::Point2d> &vp2d,vector<cv::Point3d> &vp3d)
{

    // 内参初始化
    mk[0]=0;mk[1]=0;mk[2]=0;mk[3]=0;
    mp[0]=0;mp[1]=0;mc[0]=0;mc[1]=0;
    mp0[0]=minip0[0];mp0[1]=minip0[1];

    // 外参初始化
    // 2d纠正后再pnp
    vector<cv::Point2d> p2d;
    cv::Point2d p;
    for(int i=0;i<vp2d.size();i++)
    {
        // fisheye -> rectify planar
        double x=vp2d[i].x-minip0[0];
        double y=vp2d[i].y-minip0[1];

        double r=sqrt(x*x+y*y);
        double xr=mf*x/r*tan(r/mf);
        double yr=mf*y/r*tan(r/mf);
        
        double r2=xr*xr+yr*yr;
        double dx=xr*(mk[0]*r2+mk[1]*r2*r2+mk[2]*r2*r2*r2+mk[3]*r2*r2*r2)+2*mp[0]*xr*yr+mp[1]*(r2+2*xr*xr)+mc[0]*xr+mc[1]*yr;
        double dy=yr*(mk[0]*r2+mk[1]*r2*r2+mk[2]*r2*r2*r2+mk[3]*r2*r2*r2)+mp[0]*(r2+2*yr*yr)+2*mp[1]*yr*xr+mc[0]*xr+mc[1]*yr;
        double x_=xr+dx+minip0[0];
        double y_=yr+dy+minip0[1];
        p.x=x_;
        p.y=y_;
        p2d.push_back(p);

    }

    cv::Mat cameraMatrix=cv::Mat::zeros(3, 3, CV_64FC1);
    cameraMatrix.at<double>(0, 0) = mf;
	cameraMatrix.at<double>(0, 2) = minip0[0];
	cameraMatrix.at<double>(1, 1) = mf;
	cameraMatrix.at<double>(1, 2) = minip0[1];
	cameraMatrix.at<double>(2, 2) = 1.0;
    cv::Mat r = cv::Mat::zeros(3, 1, CV_64FC1);
	cv::Mat t = cv::Mat::zeros(3, 1, CV_64FC1);
    cv::Mat in;
    cv::solvePnPRansac(vp3d,vp2d,cameraMatrix,Mat(),r,t,true,50,10,0.9, in,SOLVEPNP_EPNP);
    cv::Rodrigues(r,mMR);
    mMt=t.clone();

    Eigen::Matrix3d eR;
    cv::cv2eigen(mMR,eR);
    Eigen::Vector3d rr;
    Eigen::Quaterniond q=Eigen::Quaterniond(eR);
    mquaternion[0]=q.w();mquaternion[1]= q.x();mquaternion[2]= q.y();mquaternion[3]= q.z();
    mt[0]=mMt.at<double>(0,0);mt[1]=mMt.at<double>(0,1);mt[2]=mMt.at<double>(0,2);

    vector<bool> b;
    double error;
    int i=checkinlinersforequaldist(b,error,100);

    cout<<"initialize successful!"<<endl;
    cout<<"before solver:"<<endl<<"k: "<<mk[0]<<" "<<mk[1]<<" "<<mk[2]<<" "<<mk[3]<<" "<<endl;
    cout<<"p: "<<mp[0]<<" "<<mp[1]<<endl;
    cout<<"c: "<<mc[0]<<" "<<mc[1]<<endl;
    cout<<"p0: "<<mp0[0]<<" "<<mp0[1]<<endl;
    cout<<"R:"<<mMR<<endl<<"t:"<<mMt<<endl;
}

void exincalib::ransacinitial(int iterationnum=10,int setnum=500)     
{
    vector<int> temp;
    for(int i=0;i<mvP2D.size();i++)
    {
        temp.push_back(i);
    }
    //random_shuffle(temp.begin(), temp.end());
    rvvnum = vector< vector<int> >(iterationnum,vector<int>(setnum,0));
    int count=0;
    for(int i =0 ;i<iterationnum;i++)
    {
        // for(int j=0;j<setnum;j++)
        // {
        //     rvvnum[i][j]=temp[count];
        //     count++;
        // }
        random_shuffle(temp.begin(), temp.end());
        copy(temp.begin(),temp.begin()+setnum,rvvnum[i].begin());
    }

    miterationnum = iterationnum;
    msetnum = setnum;

}

void exincalib::optimizesolverforpoly(vector<cv::Point2d> &vp2d,vector<cv::Point3d> &vp3d)
{
    for(size_t i=0;i<vp2d.size();i++)
    {
        double p1[3];
        double p2[3];
        p1[0]=vp3d[i].x;
        p1[1]=vp3d[i].y;
        p1[2]=vp3d[i].z;
        // ceres::AngleAxisRotatePoint(mrvec, p1, p2);
        // //ceres::QuaternionRotatePoint(mquaternion,p1,p2);

        // p2[0] += mt[0]; p2[1] += mt[1]; p2[2] += mt[2];   // 世界坐标系转换到相机坐标系
        // if(p2[2]<0)
        //     continue;
        CostFunction* cost_function = polyProjecterror_exin::Create(vp2d[i],vp3d[i],mf);

        exincalibproblem.AddResidualBlock(cost_function,new CauchyLoss(0.5),ma,mp0,mb,mquaternion,mt);  //new CauchyLoss(0.5)
       // exincalibproblem.AddResidualBlock(cost_function,NULL,ma,mp0,mb);
       exincalibproblem.SetParameterBlockConstant(mb);

    }

    double cost;
    std::vector<double> residuals;
    std::vector<double> gradient;
    ceres::CRSMatrix jacobian;
    exincalibproblem.Evaluate(Problem::EvaluateOptions(), &cost, &residuals, &gradient, &jacobian);
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    //options.minimizer_progress_to_stdout = true;
    options.max_num_iterations=200;
    options.function_tolerance=1e-3;
    options.gradient_tolerance = 1e-16;
    options.minimizer_type=ceres::TRUST_REGION;
    options.trust_region_strategy_type=ceres::DOGLEG;
    //options.check_gradients=true;
    ceres::Solver::Summary summary;
    ceres::Solve(options,&exincalibproblem,&summary);
    std::cout << summary.BriefReport() << "\n";  

    cout<<"after solver:"<<endl<<"a: "<<ma[0]<<" "<<ma[1]<<" "<<ma[2]<<" "<<ma[3]<<" "<<endl;
    cout<<"p0: "<<mp0[0]<<" "<<mp0[1]<<endl;
    cout<<"b: "<<mb[0]<<" "<<mb[1]<<" "<<mb[2]<<endl;
    Eigen::Quaterniond q(mquaternion[0],mquaternion[1],mquaternion[2],mquaternion[3]);
    Eigen::Matrix3d R=q.matrix();
    cout<<"R: "<<R<<endl;
    //cout<<"qt: "<<mquaternion[0]<<" "<<mquaternion[1]<<" "<<mquaternion[2]<<" "<<mquaternion[3]<<endl;
    cout<<"t: "<<mt[0]<<" "<<mt[1]<<" "<<mt[2]<<endl;
}

void exincalib::optimizesolverforequaldist(vector<cv::Point2d> &vp2d,vector<cv::Point3d> &vp3d)
{
        for(size_t i=0;i<vp2d.size();i++)
    {
        CostFunction* cost_function = equaldistprojecterror_exin::Create(vp2d[i],vp3d[i],mf);

        exincalibproblem.AddResidualBlock(cost_function,NULL,mk,mp,mc,mp0,mquaternion,mt); //new CauchyLoss(0.5)
       // exincalibproblem.AddResidualBlock(cost_function,NULL,ma,mp0,mb);

    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.minimizer_progress_to_stdout = true;
    options.max_num_iterations=200;
    options.function_tolerance=1e-3;
    options.check_gradients=true;
    ceres::Solver::Summary summary;
    ceres::Solve(options,&exincalibproblem,&summary);
    std::cout << summary.BriefReport() << "\n";  

    cout<<"after solver:"<<endl<<"k: "<<mk[0]<<" "<<mk[1]<<" "<<mk[2]<<" "<<mk[3]<<" "<<endl;
    cout<<"p: "<<mp[0]<<" "<<mp[1]<<endl;
    cout<<"c: "<<mc[0]<<" "<<mc[1]<<endl;
    cout<<"p0: "<<mp0[0]<<" "<<mp0[1]<<endl;
    cout<<"qt: "<<mquaternion[0]<<" "<<mquaternion[1]<<" "<<mquaternion[2]<<" "<<mquaternion[3]<<endl;
    cout<<"t: "<<mt[0]<<" "<<mt[1]<<" "<<mt[2]<<endl;

}

void exincalib::ransac_solver(int &model)
{

    int count=0;
    int in;
    int bestinliners=0;
    int leasterror=0;
    
    if(model != 1)     //model ploy
    {
        cout<<"use poly model.........."<<endl;
        //caliboptimizeiniforpoly(mvP2D,mvP3D);
        //optimizesolverforpoly(mvP2D,mvP3D);
        for(int i=0;i<miterationnum;i++)
        {
            vector<cv::Point2d> rvp2d;
            vector<cv::Point3d> rvp3d;
            for(int j=0;j<msetnum;j++)
            {
                rvp2d.push_back(mvP2D[rvvnum[i][j]]);
                rvp3d.push_back(mvP3D[rvvnum[i][j]]);
                }
            cout<<endl<<"ransac optimizer solver "<<i+1<<" sets:------------------"<<endl;
        
            //caliboptimizeiniforpoly(rvp2d,rvp3d);
            caliboptimizeiniforpoly();
            vector<bool> vpinliners;
            in=checkinlinersforpoly(vpinliners,error,5);
            cout<<"the number of initial inliners is :"<<in<<endl;
            cout<<"the initial reproject error is :"<<error<<endl;
            optimizesolverforpoly(rvp2d,rvp3d);
            //vector<bool> vpinliners;

            in=checkinlinersforpoly(vpinliners,error,5);
            cout<<"the number of inliners is :"<<in<<endl;
            cout<<"the reproject error is :"<<error<<endl;

            if(error<leasterror)
            {
                leasterror=error;
                mvbBestInliers=vpinliners;
            }
        }


        // do ceres_calib using bestinliners 
        vector<cv::Point2d> vbestp2d;
        vector<cv::Point3d> vbestp3d;
        for(int i=0;i<mvP2D.size();i++)
        {
            if(!mvbBestInliers.empty() && mvbBestInliers[i])
            {
                vbestp2d.push_back(mvP2D[i]);
                vbestp3d.push_back(mvP3D[i]);
            }
        }

        cout<<endl<<"Last best inliers optimize solve!--------"<<endl;
        optimizesolverforpoly(vbestp2d,vbestp3d);
        vector<bool> vpinliners;
        minliner=checkinlinersforpoly(vpinliners,error,5);
        cout<<"the number of inliners is :"<<in<<endl;
        cout<<"the reproject error is :"<<error<<endl;

    }
    else               //model qual dist
    {
        cout<<"use equal distance model.........."<<endl;
        caliboptimizeiniforequaldist(mvP2D,mvP3D);
        double error;
        for(int i=0;i<miterationnum;i++)
        {
            vector<cv::Point2d> rvp2d;
            vector<cv::Point3d> rvp3d;
            for(int j=0;j<msetnum;j++)
            {
                rvp2d.push_back(mvP2D[rvvnum[i][j]]);
                rvp3d.push_back(mvP3D[rvvnum[i][j]]);
            }
            cout<<endl<<"ransac optimizer solver "<<i+1<<" sets:------------------"<<endl;

            caliboptimizeiniforequaldist(rvp2d,rvp3d);
            optimizesolverforequaldist(rvp2d,rvp3d);
            vector<bool> vpinliners;
            in=checkinlinersforequaldist(vpinliners,error);
            cout<<"the number of inliners is :"<<in<<endl;
        

            if(in>bestinliners)
            {
                bestinliners=in;
                mvbBestInliers=vpinliners;
            }
        }

        // do ceres_calib using bestinliners 
        vector<cv::Point2d> vbestp2d;
        vector<cv::Point3d> vbestp3d;
        for(int i=0;i<mvP2D.size();i++)
        {
            if(!mvbBestInliers.empty() && mvbBestInliers[i])
            {
                vbestp2d.push_back(mvP2D[i]);
                vbestp3d.push_back(mvP3D[i]);
            }
        }

        cout<<endl<<"Last best inliers optimize solve!--------"<<endl;
        optimizesolverforequaldist(vbestp2d,vbestp3d);

    }

}



int exincalib::checkinlinersforpoly(vector<bool> &vpinliners,double &error,int th)
{
    // Mat img = Mat::zeros(Size(mimage.cols,mimage.rows), CV_8UC3);
    // img.setTo(255);
    int inliners=0;
    int num=0;

    Eigen::Quaterniond q(mquaternion[0],mquaternion[1],mquaternion[2],mquaternion[3]);
    ///q.normalized();
    Eigen::Matrix3d R=q.toRotationMatrix();
    cout<<"R: "<<R<<endl;

    //Eigen::AngleAxisd r(mrvec[0],mrvec[1],mrvec[2]);

    error=0;
    
    cv::Mat img=mimage.clone();
    vpinliners.resize(mnum);
    for(int i=0;i<mvP2D.size();i++)
    {

        // 2d project planar f
        double x=mvP2D[i].x-mp0[0];
        double y=mvP2D[i].y-mp0[1];

        double temp=x;
        x=y;y=temp;

        Eigen::Matrix2d F1;
        F1<<mb[0],mb[1],mb[2],1;
        Eigen::Matrix2d F;
        F=F1.inverse();
        //cout<<"Fn:"<<F(0,0)<<" "<<F(0,1)<<" "<<F(1,0)<<" "<<F(0,1)<<endl;
        double xx=x*F(0,0)+y*F(0,1);
        double yy=x*F(1,0)+y*F(1,1);
        double rho=sqrt(xx*xx+yy*yy);
        double rho2=rho*rho;
        double rho3=rho*rho*rho;
        double rho4=rho*rho*rho*rho;
        
        double Z=ma[0]+ma[1]*rho2+ma[2]*rho3+ma[3]*rho4;
        Z=-Z;


        double x_=mf*xx/Z;
        double y_=mf*yy/Z;
        temp=x_;
        x_=y_;y_=temp;
        x_=x_+mp0[0];
        y_=y_+mp0[1];


        // 3d project planar f

        double p1[3];
        double p2[3];
        p1[0]=mvP3D[i].x;
        p1[1]=mvP3D[i].y;
        p1[2]=mvP3D[i].z;
        //ceres::AngleAxisRotatePoint(mrvec, p1, p2);
        double mq[4];
        mq[0]=mquaternion[0];mq[1]=mquaternion[1];mq[2]=mquaternion[2];mq[3]=mquaternion[3];
        
        double result[3];
         const double scale = double(1) / sqrt(mquaternion[0] * mquaternion[0] +
                              mquaternion[1] * mquaternion[1] +
                              mquaternion[2] * mquaternion[2] +
                              mquaternion[3] * mquaternion[3]);

        // Make unit-norm version of q.
        const double unit[4] = {
            scale * mquaternion[0],
            scale * mquaternion[1],
            scale * mquaternion[2],
            scale * mquaternion[3],
        };
        const double t2 =  mquaternion[0] * mquaternion[1];
        const double t3 =  mquaternion[0] * mquaternion[2];
        const double t4 =  mquaternion[0] * mquaternion[3];
        const double t5 = -mquaternion[1] * mquaternion[1];
        const double t6 =  mquaternion[1] * mquaternion[2];
        const double t7 =  mquaternion[1] * mquaternion[3];
        const double t8 = -mquaternion[2] * mquaternion[2];
        const double t9 =  mquaternion[2] * mquaternion[3];
        const double t1 = -mquaternion[3] * mquaternion[3];

        result[0] = 2 * ((t8 + t1) * p1[0] + (t6 - t4) * p1[1] + (t3 + t7) * p1[2]) + p1[0];  // NOLINT
        result[1] = 2 * ((t4 + t6) * p1[0] + (t5 + t1) * p1[1] + (t9 - t2) * p1[2]) + p1[1];  // NOLINT
        result[2] = 2 * ((t7 - t3) * p1[0] + (t2 + t9) * p1[1] + (t5 + t8) * p1[2]) + p1[2];


        ceres::QuaternionRotatePoint(mq,p1,p2);
        // Eigen::Vector3d pp1,pp2;
        // pp1<<p1[0],p1[1],p1[2];
        // pp2=R*pp1;
        // pp2[0] += mt[0]; pp2[1] += mt[1]; pp2[2] += mt[2];

        p2[0] += mt[0]; p2[1] += mt[1]; p2[2] += mt[2];   // 世界坐标系转换到相机坐标系

        // double xp=mf*pp2[0]/pp2[2]+mp0[0];
        // double yp=mf*pp2[1]/pp2[2]+mp0[1];

        double xp=mf*p2[0]/p2[2]+mp0[0];
        double yp=mf*p2[1]/p2[2]+mp0[1];;



        //plot
        cv::Point P1;
        P1.x=x_;P1.y=y_;
        cv::circle(img,P1,1,cv::Scalar(255, 0, 0),-1);     // 2d points
        cv::Point P2;
        P2.x=xp;P2.y=yp;
        cv::circle(img,P2,1,cv::Scalar(0, 0, 255),-1);   // 2d points

        double dist=sqrt((x_-xp)*(x_-xp)+(y_-yp)*(y_-yp));
        
            
        if(dist <th)
        {
            inliners++;
            vpinliners[i]=true;
            
        }
        else
        {
            vpinliners[i]=false;
        }
        if(dist<20)
        {
        error=error+dist;
        num++;
        }
    }

    error=error/num;
    cv::putText(img,to_string(inliners),Point(500,1800),FONT_HERSHEY_PLAIN,20,cv::Scalar(255, 255, 255),10,8);
    cv::namedWindow("camera", CV_WINDOW_NORMAL);
    imshow("camera",img);
    //imwrite("/home/jiangfan/桌面/5_14/G1/1-2/op1.jpg",img);
    waitKey(1000);

    return inliners;
}

void exincalib::savereasults(string path)
{
    ofstream of;
    of.open(path.data());
    of<<"pol: "<<ma[0]<<" "<<0 <<" "<<ma[1]<<" "<<ma[2]<<" "<<ma[3]<<endl;
    of<<"x0,y0: "<<mp0[0]<< " " <<mp0[1]<<endl;
    of<<"c,d,e:"<<mb[0]<<" "<<mb[1]<<" "<<mb[2]<<endl;
    //of<<"angle-aixs:"<<mrvec[0]<<" "<<mrvec[1]<<" "<<mrvec[1]<<endl;
    of<<"mquaternion:"<<mquaternion[0]<<" "<<mquaternion[1]<<" "<<mquaternion[2]<<" "<<mquaternion[3]<<endl;
    Eigen::Matrix3d R;
    of<<"t:"<<mt[0]<<" "<<mt[1]<<" "<<mt[2]<<endl;
    of<<"error:"<<error<<endl;
    of.close();
}

int exincalib::checkinlinersforequaldist(vector<bool> &vpinliners,double &error,int th)
{
    Mat img = Mat::zeros(Size(mimage.cols,mimage.rows), CV_8UC3);
    img.setTo(255);

    int inliners=0;
    vpinliners.resize(mnum);
    for(int i=0;i<mvP2D.size();i++)
    {
        // fisheye -> rectify planar
        double x=mvP2D[i].x-mp0[0];
        double y=mvP2D[i].y-mp0[1];

        double r=sqrt(x*x+y*y);
        double xr=mf*x/r*tan(r/mf);
        double yr=mf*y/r*tan(r/mf);
        
        double r2=xr*xr+yr*yr;
        double dx=xr*(mk[0]*r2+mk[1]*r2*r2+mk[2]*r2*r2*r2+mk[3]*r2*r2*r2)+2*mp[0]*xr*yr+mp[1]*(r2+2*xr*xr)+mc[0]*xr+mc[1]*yr;
        double dy=yr*(mk[0]*r2+mk[1]*r2*r2+mk[2]*r2*r2*r2+mk[3]*r2*r2*r2)+mp[0]*(r2+2*yr*yr)+2*mp[1]*yr*xr+mc[0]*xr+mc[1]*yr;
        double x_=xr+dx+mp0[0];
        double y_=yr+dy+mp0[1];

        // 3D -> planar
        double p1[3];
        double p2[3];
        p1[0]=mvP3D[i].x;
        p1[1]=mvP3D[i].y;
        p1[2]=mvP3D[i].z;
        ceres::QuaternionRotatePoint(mquaternion,p1,p2);
        p2[0] += mt[0]; p2[1] += mt[1]; p2[2] += mt[2];   // 世界坐标系转换到相机坐标系

        double xp=mf*p2[0]/p2[2]+mp0[0];
        double yp=mf*p2[1]/p2[2]+mp0[1];

        //plot
        cv::Point P1;
        P1.x=x_;P1.y=y_;
        cv::circle(img,P1,5,cv::Scalar(255, 0, 0),-1);
        cv::Point P2;
        P2.x=xp;P2.y=yp;
        cv::circle(img,P2,5,cv::Scalar(0, 0, 255),-1);
        
        double dist=sqrt((x_-xp)*(x_-xp)+(y_-yp)*(y_-yp));
            
        if(dist <th)
        {
            inliners++;
            vpinliners[i]=true;
        }
        else
        {
            vpinliners[i]=false;
        }
    }
    cv::putText(img,to_string(inliners),Point(500,1800),FONT_HERSHEY_PLAIN,20,cv::Scalar(0, 0, 0),10,8);
    cv::namedWindow("camera", CV_WINDOW_NORMAL);
    imshow("camera",img);
     waitKey(0);
    //getchar();
    return inliners;
}

// void exincalib::Fisheye_C_PolyProjetion(cv::Point2d &p2d,cv::Point3d &p3d)
// {
//         double x=p3d.x-mp0[0];
//         double y=p2d.y-mp0[1];
//         double xx=x*mb[0]+y*mb[1];
//         double yy=x*mb[2]+y;
//         double rho=sqrt(xx*xx+yy*yy);
//         double rho2=rho*rho;
//         double rho3=rho2*rho;
//         double rho4=rho2*rho2;
        
//         double Z=ma[0]+ma[1]*rho2+ma[2]*rho3+ma[3]*rho4;
//         Z=-Z;
//         p3d.x=xx;
//         p3d.y=yy;
//         p3d.z=Z;
// }

// void Fisheye_C_EqualdistProjection(cv::Point2d &p2d,cv::Point3d &p3d); 
// {

// }

void exincalib::CalculateRotationMatFromEulerAngle(double Rx, double Ry, double Rz, double * R)   //Rz*Ry*Rx
{
    double cRx, cRy, cRz, sRx, sRy, sRz;
    cRx = cos(Rx); cRy = cos(Ry); cRz = cos(Rz);
    sRx = sin(Rx); sRy = sin(Ry); sRz = sin(Rz);
    R[0] = cRz * cRy;
    R[1] = cRz * sRy * sRx - sRz * cRx;
    R[2] = cRz * sRy * cRx + sRz * sRx;
    R[3] = sRz * cRy;
    R[4] = sRz * sRy* sRx + cRz * cRx;
    R[5] = sRz * sRy* cRx - cRz * sRx;
    R[6] = -sRy;
    R[7] = cRy * sRx;
    R[8] = cRy * cRx;
}

void exincalib::CalculateEulerAngleFromRotationMat(double &Rx, double &Ry, double &Rz, double * R)
{
    Rx=atan2(R[7],R[8]);
    Ry=atan2(-R[6],sqrt(R[7]*R[7]+R[8]*R[8]));
    Rz=atan2(R[3],R[0]);
}