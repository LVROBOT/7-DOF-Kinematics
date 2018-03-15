// DH_Model_Example.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "complex"
#include "Eigen/Dense"
using namespace std;
#define PI 3.1415926535897932384626433832795

// ƽ��
Eigen::Matrix4d translate(Eigen::Vector3d xyz)               //ƽ������
{
	//Identify()����Ϊ��λ�Խ��������Ͻǵ����½ǶԽ���ֵΪ1������Ϊ0�� Matrix4d��Ӧ4��4�о���
	Eigen::Matrix4d T = Eigen::Matrix4d::Identity(4, 4);

	T(0, 3) = xyz(0);
	T(1, 3) = xyz(1);
	T(2, 3) = xyz(2);

	return T;
}

// ��ת
Eigen::Matrix4d rotation(char axis, double radian)      //��ת����
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

// ΢��
Eigen::VectorXd tr2diff(Eigen::MatrixXd T1, Eigen::MatrixXd T2)     //λ��֮��
{
	Eigen::VectorXd diff(6);
	Eigen::VectorXd vCross = Eigen::VectorXd::Zero(3);

	for (int i = 0; i<3; i++)
	{
		Eigen::Vector3d vT1, vT2;
		vT1 << T1(0, i), T1(1, i), T1(2, i);
		vT2 << T2(0, i), T2(1, i), T2(2, i);

		vCross += vT1.cross(vT2);      //��ȡ������
	}

	for (int i = 0; i<3; i++)
	{
		diff(i) = T2(i, 3) - T1(i, 3);
		diff(i + 3) = 0.5*vCross(i);
	}

	return diff;
}

// ����
Eigen::MatrixXd  DH_Foward(Eigen::MatrixXd DH, Eigen::VectorXd angle)    //���˶�ѧ
{
	Eigen::Matrix4d T_60 = Eigen::Matrix4d::Identity(4, 4);

	for (int i = 0; i<DH.rows(); i++)
	{
		Eigen::Vector3d xyz;
		xyz << DH(i, 1), 0, DH(i, 3);

		DH(i, 2) = DH(i, 2) + angle(i);

		T_60 = T_60 * rotation('z', DH(i, 2)) * translate(xyz) * rotation('x', DH(i, 0));     //  DHģ��

		//xyz << MDH(i, 1), 0, MDH(i, 3);

		//MDH(i, 2) = MDH(i, 2) + angle(i);

		//T_60 = T_60 * rotation('x', MDH(i, 0)) * translate(xyz) * rotation('z', MDH(i, 2));     //  M_DHģ��
	}

	return T_60;
}

// �ſ˱�n
Eigen::MatrixXd DH_Jacobn(Eigen::MatrixXd DH, Eigen::VectorXd angle)      //��е��ĩ���ſ˱Ⱦ���
{
	//Eigen::MatrixXd Jacobn = Eigen::MatrixXd::Identity(6,6);  
	Eigen::MatrixXd Jacobn = Eigen::MatrixXd::Identity(6, 7);    //�������ڻ�е���ڵѿ����ռ�����ɶ��������������ڲ����۹ؽ�����
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

// �ſ˱�0
Eigen::MatrixXd DH_Jacob0(Eigen::MatrixXd DH, Eigen::VectorXd angle)     //��е�ۻ����ſ˱Ⱦ���
{
	Eigen::MatrixXd Jacob0 = DH_Jacobn(DH, angle);         //Jacob0��6*7����
	Eigen::MatrixXd T = DH_Foward(DH, angle);

	Eigen::MatrixXd Transform(6, 6);
	Transform << T(0, 0), T(0, 1), T(0, 2), 0, 0, 0,
		T(1, 0), T(1, 1), T(1, 2), 0, 0, 0,
		T(2, 0), T(2, 1), T(2, 2), 0, 0, 0,
		0, 0, 0, T(0, 0), T(0, 1), T(0, 2),
		0, 0, 0, T(1, 0), T(1, 1), T(1, 2),
		0, 0, 0, T(2, 0), T(2, 1), T(2, 2);
	Jacob0 = Transform * Jacob0;           //Jacob0�õ��Ļ���6*7����

	return Jacob0;
}

// ����
Eigen::VectorXd DH_Inverse(Eigen::MatrixXd DH, Eigen::MatrixXd T, Eigen::VectorXd ref)         //�˶�ѧ���
{
	//Eigen::VectorXd q = Eigen::VectorXd::Zero(6);
	//Eigen::VectorXd q = Eigen::VectorXd::Zero(7);               

	Eigen::VectorXd q = ref;

	double nm = 1.0;
	double mini = 1e-12f;
	int count = 0;
	int limitNum = 1000;

	double qlimit[7] = { 170, 120, 170, 120, 170, 120, 175 };  //�¼���Ĺؽڽ�����

	while (nm > mini)
	{
		Eigen::VectorXd e = tr2diff(DH_Foward(DH, q), T);
		//cout << "rows and cols number of e is:" << e.rows() << "      " << e.cols() << endl;

		Eigen::MatrixXd Jacob0 = DH_Jacob0(DH, q);
		//cout << "rows and cols number of jacob0.pinv().transpose() is:" << Jacob0.pinv().transpose().rows() << "     " << Jacob0.pinv().transpose().cols() << endl;

	
		Eigen::VectorXd dq = (Jacob0.pinv().transpose()) *e;      //  Jacob0.pinv().transpose()��7*6����e��6*1����

		q = q + dq;

		/*for (int i = 0; i<q.rows(); i++)
		{
			if (abs(q(i)) > 2 * PI)
			{
				q(i) = q(i) - ((int)(q(i) / (2 * PI))) * 2 * PI;
			}
		}*/

		nm = dq.norm();    //norm,��dq��ʸ������,��ʸ��dq�ĳ���

		count++;
		if (count > limitNum)
		{
			std::cout << "�����޷�����\n" << endl;
			break;
		}
	}
		for (int i = 0; i<q.rows(); i++)
		{
			if (abs(q(i)) > qlimit[i]* PI / 180)
			{
				cout << "�ؽ� " << i+1 << "  ����" << endl;
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

	//DHģ��
	DH << -PI / 2, 0, 0, 340,
		PI / 2, 0, 0, 0,
		PI / 2, 0, 0, 400,
		-PI / 2, 0, 0, 0,
		-PI / 2, 0, 0, 400,
		PI / 2, 0, 0, 0,
		0, 0, 0, 152;

	Eigen::VectorXd angle(7);
	angle << 30, 50 , 60, 42, 51, 85, 26;        //�ṩһ��Ŀ��λ��
	angle = angle * PI / 180;

	Eigen::VectorXd ref(7);
	ref << 32, 48, 58, 45, 47, 83, 27;          //��һ���ڹؽڽ���Ϊ���Ĳο���
	ref = ref * PI / 180;


	Eigen::Matrix4d T1 = DH_Foward(DH,angle);         //��ȡĿ����λ��
	std::cout<<"������ĩ��λ��Ϊ��\n"<<T1<<"\n"<<endl;            

	Eigen::MatrixXd Jacobn = DH_Jacobn(DH,angle);
	std::cout<<"��ǰλ���µ�ĩ���ſ˱�Ϊ��\n"<<Jacobn<<"\n"<<endl;

	Eigen::MatrixXd Jacob0 = DH_Jacob0(DH,angle);
	std::cout<<"��ǰλ���µĻ����ſ˱�Ϊ: \n"<<Jacob0<<"\n"<<endl;

	Eigen::VectorXd angle2 = DH_Inverse(DH,T1,ref);         //�������õĹؽڽ�
	std::cout<<"��ǰλ���·���ؽڽ�Ϊ��\n"<<angle2*180/PI<<"\n"<<endl;        

	Eigen::Matrix4d T2 = DH_Foward(DH,angle2);       //��ȡ�������ĹؽڽǶ�Ӧ��λ�ˣ���Ŀ���λ�˽��бȽϣ���һ�£��������ȷ
	std::cout<<"�����Ӧ��ĩ��λ��Ϊ��\n"<<T2<<"\n"<<endl;

	printf("Press Enter to exit...\n");
	getchar();
}