// DH_Model_Example.cpp : 定义控制台应用程序的入口点。
//

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "complex"
#include "Eigen/Dense"
using namespace std;
#define PI 3.1415926535897932384626433832795

// 平移
Eigen::Matrix4d translate(Eigen::Vector3d xyz)               //平移算子
{
	//Identify()函数为单位对角阵，由左上角到右下角对角线值为1，其余为0。 Matrix4d对应4行4列矩阵。
	Eigen::Matrix4d T = Eigen::Matrix4d::Identity(4, 4);

	T(0, 3) = xyz(0);
	T(1, 3) = xyz(1);
	T(2, 3) = xyz(2);

	return T;
}

// 旋转
Eigen::Matrix4d rotation(char axis, double radian)      //旋转算子
{
	Eigen::Matrix4d T = Eigen::Matrix4d::Identity(4, 4);

	double ct = cos(radian);
	double st = sin(radian);

	switch (axis)
	{
	case 'x':
	case 'X':
		T << 1, 0, 0, 0,
			0, ct, -st, 0,
			0, st, ct, 0,
			0, 0, 0, 1;
		break;
	case 'y':
	case 'Y':
		T << ct, 0, st, 0,
			0, 1, 0, 0,
			-st, 0, ct, 0,
			0, 0, 0, 1;
		break;
	case 'z':
	case 'Z':
		T << ct, -st, 0, 0,
			st, ct, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1;
		break;
	}

	return T;
}

// 微分
Eigen::VectorXd tr2diff(Eigen::MatrixXd T1, Eigen::MatrixXd T2)     //位姿之差
{
	Eigen::VectorXd diff(6);
	Eigen::VectorXd vCross = Eigen::VectorXd::Zero(3);

	for (int i = 0; i<3; i++)
	{
		Eigen::Vector3d vT1, vT2;
		vT1 << T1(0, i), T1(1, i), T1(2, i);
		vT2 << T2(0, i), T2(1, i), T2(2, i);

		vCross += vT1.cross(vT2);      //求取向量积
	}

	for (int i = 0; i<3; i++)
	{
		diff(i) = T2(i, 3) - T1(i, 3);
		diff(i + 3) = 0.5*vCross(i);
	}

	return diff;
}

// 正解
Eigen::MatrixXd  DH_Foward(Eigen::MatrixXd DH, Eigen::VectorXd angle)    //正运动学
{
	Eigen::Matrix4d T_60 = Eigen::Matrix4d::Identity(4, 4);

	for (int i = 0; i<DH.rows(); i++)
	{
		Eigen::Vector3d xyz;
		xyz << DH(i, 1), 0, DH(i, 3);

		DH(i, 2) = DH(i, 2) + angle(i);

		T_60 = T_60 * rotation('z', DH(i, 2)) * translate(xyz) * rotation('x', DH(i, 0));     //  DH模型

		//xyz << MDH(i, 1), 0, MDH(i, 3);

		//MDH(i, 2) = MDH(i, 2) + angle(i);

		//T_60 = T_60 * rotation('x', MDH(i, 0)) * translate(xyz) * rotation('z', MDH(i, 2));     //  M_DH模型
	}

	return T_60;
}

// 雅克比n
Eigen::MatrixXd DH_Jacobn(Eigen::MatrixXd DH, Eigen::VectorXd angle)      //机械臂末端雅克比矩阵
{
	//Eigen::MatrixXd Jacobn = Eigen::MatrixXd::Identity(6,6);  
	Eigen::MatrixXd Jacobn = Eigen::MatrixXd::Identity(6, 7);    //行数等于机械臂在笛卡尔空间的自由度数量，列数等于操作臂关节数量
	Eigen::Matrix4d T_i = Eigen::Matrix4d::Identity(4, 4);
	for (int i = DH.rows() - 1; i >= 0; i--)
	{
		Eigen::Vector3d xyz;
		xyz << DH(i, 1), 0, DH(i, 3);
		DH(i, 2) = DH(i, 2) + angle(i);
		T_i = rotation('z', DH(i, 2)) * translate(xyz) * rotation('x', DH(i, 0)) * T_i;

		Eigen::Vector3d n, o, a, p;
		a << T_i(0, 0), T_i(1, 0), T_i(2, 0);
		o << T_i(0, 1), T_i(1, 1), T_i(2, 1);
		n << T_i(0, 2), T_i(1, 2), T_i(2, 2);
		p << T_i(0, 3), T_i(1, 3), T_i(2, 3);

		Jacobn(0, i) = -a(0) * p(1) + a(1) * p(0);
		Jacobn(1, i) = -o(0) * p(1) + o(1) * p(0);
		Jacobn(2, i) = -n(0) * p(1) + n(1) * p(0);
		Jacobn(3, i) = a(2);
		Jacobn(4, i) = o(2);
		Jacobn(5, i) = n(2);

	}


	return Jacobn;
}

// 雅克比0
Eigen::MatrixXd DH_Jacob0(Eigen::MatrixXd DH, Eigen::VectorXd angle)     //机械臂基座雅克比矩阵
{
	Eigen::MatrixXd Jacob0 = DH_Jacobn(DH, angle);         //Jacob0是6*7矩阵
	Eigen::MatrixXd T = DH_Foward(DH, angle);

	Eigen::MatrixXd Transform(6, 6);
	Transform << T(0, 0), T(0, 1), T(0, 2), 0, 0, 0,
		T(1, 0), T(1, 1), T(1, 2), 0, 0, 0,
		T(2, 0), T(2, 1), T(2, 2), 0, 0, 0,
		0, 0, 0, T(0, 0), T(0, 1), T(0, 2),
		0, 0, 0, T(1, 0), T(1, 1), T(1, 2),
		0, 0, 0, T(2, 0), T(2, 1), T(2, 2);
	Jacob0 = Transform * Jacob0;           //Jacob0得到的还是6*7矩阵

	return Jacob0;
}

// 反解
Eigen::VectorXd DH_Inverse(Eigen::MatrixXd DH, Eigen::MatrixXd T, Eigen::VectorXd ref)         //运动学逆解
{
	//Eigen::VectorXd q = Eigen::VectorXd::Zero(6);
	//Eigen::VectorXd q = Eigen::VectorXd::Zero(7);               

	Eigen::VectorXd q = ref;

	double nm = 1.0;
	double mini = 1e-12f;
	int count = 0;
	int limitNum = 1000;

	double qlimit[7] = { 170, 120, 170, 120, 170, 120, 175 };  //新加入的关节角限制

	while (nm > mini)
	{
		Eigen::VectorXd e = tr2diff(DH_Foward(DH, q), T);
		//cout << "rows and cols number of e is:" << e.rows() << "      " << e.cols() << endl;

		Eigen::MatrixXd Jacob0 = DH_Jacob0(DH, q);
		//cout << "rows and cols number of jacob0.pinv().transpose() is:" << Jacob0.pinv().transpose().rows() << "     " << Jacob0.pinv().transpose().cols() << endl;

	
		Eigen::VectorXd dq = (Jacob0.pinv().transpose()) *e;      //  Jacob0.pinv().transpose()是7*6矩阵，e是6*1矩阵

		q = q + dq;

		/*for (int i = 0; i<q.rows(); i++)
		{
			if (abs(q(i)) > 2 * PI)
			{
				q(i) = q(i) - ((int)(q(i) / (2 * PI))) * 2 * PI;
			}
		}*/

		nm = dq.norm();    //norm,求dq的矢量范数,即矢量dq的长度

		count++;
		if (count > limitNum)
		{
			std::cout << "反解无法收敛\n" << endl;
			break;
		}
	}
		for (int i = 0; i<q.rows(); i++)
		{
			if (abs(q(i)) > qlimit[i]* PI / 180)
			{
				cout << "关节 " << i+1 << "  超限" << endl;
			}
		}
	
	return q;
}


int main()
{
	//DH is DH model of the manipulator. 
	// rows is DOF, cols is alpha/a/theta/d.
	// the angle uses radian.
	// the distance uses mm.
	Eigen::MatrixXd DH(7, 4);

	//DH模型
	DH << -PI / 2, 0, 0, 340,
		PI / 2, 0, 0, 0,
		PI / 2, 0, 0, 400,
		-PI / 2, 0, 0, 0,
		-PI / 2, 0, 0, 400,
		PI / 2, 0, 0, 0,
		0, 0, 0, 152;

	Eigen::VectorXd angle(7);
	angle << 30, 50 , 60, 42, 51, 85, 26;        //提供一个目标位姿
	angle = angle * PI / 180;

	Eigen::VectorXd ref(7);
	ref << 32, 48, 58, 45, 47, 83, 27;          //上一周期关节角作为逆解的参考角
	ref = ref * PI / 180;


	Eigen::Matrix4d T1 = DH_Foward(DH,angle);         //求取目标点的位姿
	std::cout<<"其正解末端位姿为：\n"<<T1<<"\n"<<endl;            

	Eigen::MatrixXd Jacobn = DH_Jacobn(DH,angle);
	std::cout<<"当前位姿下的末端雅克比为：\n"<<Jacobn<<"\n"<<endl;

	Eigen::MatrixXd Jacob0 = DH_Jacob0(DH,angle);
	std::cout<<"当前位姿下的基座雅克比为: \n"<<Jacob0<<"\n"<<endl;

	Eigen::VectorXd angle2 = DH_Inverse(DH,T1,ref);         //经逆解求得的关节角
	std::cout<<"当前位姿下反解关节角为：\n"<<angle2*180/PI<<"\n"<<endl;        

	Eigen::Matrix4d T2 = DH_Foward(DH,angle2);       //求取逆解出来的关节角对应的位姿，与目标点位姿进行比较，若一致，则逆解正确
	std::cout<<"反解对应的末端位姿为：\n"<<T2<<"\n"<<endl;

	printf("Press Enter to exit...\n");
	getchar();
}