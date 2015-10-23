//Copyright (C) 2012-2015, guardiancrow
// 三原色＋白色から変換マトリックスを求める方法
// boostを使用

#include <iostream>
#include "boost/multiprecision/cpp_dec_float.hpp"
#include "boost/numeric/ublas/lu.hpp"

using namespace std;
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;

//素直にテンプレートにすべきですね
#ifdef USE_DEC_FLOAT
#define double cpp_dec_float_50
#endif

//逆行列を求めます
void InvertMatrix(matrix<double>&mat, matrix<double>&matinv)
{
	permutation_matrix<> pm(mat.size1());

	matrix<double> copiyMatrix = matrix<double>(mat);
	matrix<double> invMatrix = identity_matrix<double>(mat.size1());

	lu_factorize(copiyMatrix, pm);
	lu_substitute(copiyMatrix, pm, invMatrix);

	matinv.resize(invMatrix.size1(), invMatrix.size2(), false);
	matinv.assign(invMatrix);
}

//三原色と白色から変換マトリックスを求めます
int CalcXYZMatrix(double red_x, double red_y, double green_x, double green_y, double blue_x, double blue_y, double white_x, double white_y, matrix<double> &mat3x3)
{
	//マトリックスを求めます
	//rRx/Ry + gGx/Gy + bBx/By = Wx/Wy
	//r      + g      + b      = 1
	//rRz/Ry + gGz/Gy + bBz/By = Wz/Wy
	//の連立方程式からr,g,bを求めます
	matrix<double> in_mat(3,4);

	in_mat(0,0) = red_x / red_y;
	in_mat(0,1) = green_x / green_y;
	in_mat(0,2) = blue_x / blue_y;
	in_mat(0,3) = white_x / white_y;
	in_mat(1,0) = 1.;
	in_mat(1,1) = 1.;
	in_mat(1,2) = 1.;
	in_mat(1,3) = 1.;
	in_mat(2,0) = (1 - red_x - red_y) / red_y;
	in_mat(2,1) = (1 - green_x - green_y) / green_y;
	in_mat(2,2) = (1 - blue_x - blue_y) / blue_y;
	in_mat(2,3) = (1 - white_x - white_y) / white_y;

	//Gauss-Jordan
	for(int k = 0; k < in_mat.size1(); k++){
		double p = in_mat(k,k);
		int i;
		for (i = k; i <  in_mat.size2(); i++){
			in_mat(k,i) /= p;
		}

		for (i = 0; i < in_mat.size1(); i++){
			if(i != k){
				double d = in_mat(i,k);
				for(int j = k; j < in_mat.size2(); j++)
					in_mat(i,j) -= d * in_mat(k,j);
			}
		}
	}

	matrix<double> out_mat(3,3);
	out_mat(0,0) = in_mat(0,3) * red_x / red_y;
	out_mat(0,1) = in_mat(1,3) * green_x / green_y;
	out_mat(0,2) = in_mat(2,3) * blue_x / blue_y;
	out_mat(1,0) = in_mat(0,3);
	out_mat(1,1) = in_mat(1,3);
	out_mat(1,2) = in_mat(2,3);
	out_mat(2,0) = in_mat(0,3) * (1 - red_x - red_y) / red_y;
	out_mat(2,1) = in_mat(1,3) * (1 - green_x - green_y) / green_y;
	out_mat(2,2) = in_mat(2,3) * (1 - blue_x - blue_y) / blue_y;

	mat3x3 = out_mat;
}

//単一テスト
void transform_test(double R, double G, double B, matrix<double>&in_mat, matrix<double>&out_mat)
{
	double X, Y, Z;
	X = Y = Z = 0.;
	cout << setprecision(10) << "source :R = " << R << " : G = " << G << " : B = " << B << endl;
	X = in_mat(0,0) * R + in_mat(0,1) * G + in_mat(0,2) * B;
	Y = in_mat(1,0) * R + in_mat(1,1) * G + in_mat(1,2) * B;
	Z = in_mat(2,0) * R + in_mat(2,1) * G + in_mat(2,2) * B;
	R = out_mat(0,0) * X + out_mat(0,1) * Y + out_mat(0,2) * Z;
	G = out_mat(1,0) * X + out_mat(1,1) * Y + out_mat(1,2) * Z;
	B = out_mat(2,0) * X + out_mat(2,1) * Y + out_mat(2,2) * Z;

	cout << setprecision(10) << "trans  :X = " << X << " : Y = " << Y << " : Z = " << Z << endl;
	cout << setprecision(10) << "retrans:R = " << R << " : G = " << G << " : B = " << B << endl;
}

//階調テスト
void transform_test_steps(unsigned int steps, matrix<double>&in_mat, matrix<double>&out_mat)
{
	unsigned int i;
	double R, G, B, X, Y, Z;
	cout << "transform_test_steps : " << steps << " (+1) gradation" << endl;
	for(i=0;i<=steps;i++){
		R = G = B = X = Y = Z = 0.;
		R = G = B = (double)i / (double)steps * (double)1.0;
		cout << i << ":" << endl;
		cout << setprecision(10) << "source :(R,G,B)=(" << R << "," << G << "," << B << ")" << endl;
		X = in_mat(0,0) * R + in_mat(0,1) * G + in_mat(0,2) * B;
		Y = in_mat(1,0) * R + in_mat(1,1) * G + in_mat(1,2) * B;
		Z = in_mat(2,0) * R + in_mat(2,1) * G + in_mat(2,2) * B;
		R = out_mat(0,0) * X + out_mat(0,1) * Y + out_mat(0,2) * Z;
		G = out_mat(1,0) * X + out_mat(1,1) * Y + out_mat(1,2) * Z;
		B = out_mat(2,0) * X + out_mat(2,1) * Y + out_mat(2,2) * Z;
		cout << setprecision(10) << "trans  :(X,Y,Z)=(" << X << "," << Y << "," << Z << ")" << endl;
		cout << setprecision(10) << "retrans:(R,G,B)=(" << R << "," << G << "," << B << ")" << endl;
	}
}

//順行列・逆行列のダンプ
void dump(matrix<double>&in_mat, matrix<double>&out_mat)
{
	cout << setprecision(10) << "X = " << in_mat(0,0) << " * Red + " << in_mat(0,1) << " * Green + " << in_mat(0,2) << " * Blue" << endl;
	cout << setprecision(10) << "Y = " << in_mat(1,0) << " * Red + " << in_mat(1,1) << " * Green + " << in_mat(1,2) << " * Blue" << endl;
	cout << setprecision(10) << "Z = " << in_mat(2,0) << " * Red + " << in_mat(2,1) << " * Green + " << in_mat(2,2) << " * Blue" << endl;
	cout << setprecision(10) << "Red   = " << out_mat(0,0) << " * X + " << out_mat(0,1) << " * Y + " << out_mat(0,2) << " * Z" << endl;
	cout << setprecision(10) << "Green = " << out_mat(1,0) << " * X + " << out_mat(1,1) << " * Y + " << out_mat(1,2) << " * Z" << endl;
	cout << setprecision(10) << "Blue  = " << out_mat(2,0) << " * X + " << out_mat(2,1) << " * Y + " << out_mat(2,2) << " * Z" << endl;
}

int main(int argc, char **argv)
{
	double red_x;
	double red_y;
	double green_x;
	double green_y;
	double blue_x;
	double blue_y;
	double white_x;
	double white_y;

	matrix<double> matSRGB;
	matrix<double> matSRGBinv;
	matrix<double> matAdobeRGB;
	matrix<double> matAdobeRGBinv;
	matrix<double> matCustom;
	matrix<double> matCustominv;

	double R,G,B;

	matSRGB = matrix<double>(3,3);
	matSRGBinv = matrix<double>(3,3);
	matAdobeRGB = matrix<double>(3,3);
	matAdobeRGBinv = matrix<double>(3,3);

	matCustom = matrix<double>(3,3);
	matCustominv = matrix<double>(3,3);

	//sRGB (D65)
	red_x = 0.64;
	red_y = 0.33;
	green_x = 0.30;
	green_y = 0.60;
	blue_x = 0.15;
	blue_y = 0.06;
	white_x = 0.3127;
	white_y = 0.3290;

	CalcXYZMatrix(red_x, red_y, green_x, green_y, blue_x, blue_y, white_x, white_y, matSRGB);
	InvertMatrix(matSRGB, matSRGBinv);

	cout << "-----=====-----" << endl;
	cout << "sRGB(Illuminant D65)" << endl;
	cout << "Red(0.64, 0.33)" << endl;
	cout << "Green(0.30, 0.60)" << endl;
	cout << "Blue(0.15, 0.06)" << endl;
	cout << "White(0.3127, 0.3290)" << endl;
	dump(matSRGB, matSRGBinv);

	//AdobeRGB
	red_x = 0.64;
	red_y = 0.33;
	green_x = 0.21;
	green_y = 0.71;
	blue_x = 0.15;
	blue_y = 0.06;
	white_x = 0.3127;
	white_y = 0.3290;

	CalcXYZMatrix(red_x, red_y, green_x, green_y, blue_x, blue_y, white_x, white_y, matAdobeRGB);
	InvertMatrix(matAdobeRGB, matAdobeRGBinv);

	cout << "-----=====-----" << endl;
	cout << "AdobeRGB" << endl;
	cout << "Red(0.64, 0.33)" << endl;
	cout << "Green(0.21, 0.71)" << endl;
	cout << "Blue(0.15, 0.06)" << endl;
	cout << "White(0.3127, 0.3290)" << endl;
	dump(matAdobeRGB, matAdobeRGBinv);

	//sRGBの色相を入れ替えたもの(R<->G)
	red_x = 0.30;
	red_y = 0.60;
	green_x = 0.64;
	green_y = 0.33;
	blue_x = 0.15;
	blue_y = 0.06;
	white_x = 0.3127;
	white_y = 0.3290;

	CalcXYZMatrix(red_x, red_y, green_x, green_y, blue_x, blue_y, white_x, white_y, matCustom);
	InvertMatrix(matCustom, matCustominv);

	cout << "-----=====-----" << endl;
	cout << "sRGB(D65) Red <-> Green" << endl;
	cout << "Red(0.30, 0.60)" << endl;
	cout << "Green(0.64, 0.33)" << endl;
	cout << "Blue(0.15, 0.06)" << endl;
	cout << "White(0.3127, 0.3290)" << endl;
	dump(matCustom, matCustominv);
	cout << "-----=====-----" << endl;

	//変換テスト(1) sRGB -> Adobe RGB
	//(R,G,B) = (0.0,0.0,0.0)
	R = G = B = 0.;
	cout << "-----=====-----" << endl;
	cout << "sRGB -> AdobeRGB : (R,G,B) = (0.0,0.0,0.0)" << endl;
	transform_test(R, G, B, matSRGB, matAdobeRGBinv);
	cout << "-----=====-----" << endl;

	//(R,G,B) = (1.0,1.0,1.0)
	R = G = B = 1.;
	cout << "-----=====-----" << endl;
	cout << "sRGB -> AdobeRGB : (R,G,B) = (1.0,1.0,1.0)" << endl;
	transform_test(R, G, B, matSRGB, matAdobeRGBinv);
	cout << "-----=====-----" << endl;

	//変換テスト(2) sRGB -> 入れ替えたもの
	R = 0.25;
	G = 0.5;
	B = 0.75;
	cout << "-----=====-----" << endl;
	cout << "sRGB -> RGchange : (R,G,B) = (0.25,0,5,0.75)" << endl;
	transform_test(R, G, B, matSRGB, matCustominv);
	cout << "-----=====-----" << endl;

	//階調変換テスト 8段階 sRGB -> Adobe RGB
	cout << "-----=====-----" << endl;
	cout << "sRGB -> AdobeRGB" << endl;
	transform_test_steps(8, matSRGB, matAdobeRGBinv);
	cout << "-----=====-----" << endl;

	return 0;
}
