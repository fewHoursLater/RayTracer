#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <limits>
#include <random>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include "CImg.h"


#define Pi 3.1415926535f;
constexpr float BIAS = 1e-8f;
constexpr float epsi = 1e-3f;

using namespace cimg_library;
using namespace std;

class Vector_3_float
{
public:
	Vector_3_float();
	Vector_3_float(const float, const float, const float);
	Vector_3_float(const Vector_3_float&);
	Vector_3_float& operator=(const Vector_3_float&);
	~Vector_3_float();

	Vector_3_float operator+(const Vector_3_float&);
	Vector_3_float operator-(const Vector_3_float&);
	float operator*(const Vector_3_float&);
	Vector_3_float operator*(const float&);
	float operator^(const Vector_3_float&);

	Vector_3_float V_product(const Vector_3_float&);

	void normalize(void);
	float magnitude(void);

	float get_v1(void);
	float get_v2(void);
	float get_v3(void);

private:
	float vec1_;
	float vec2_;
	float vec3_;

};


Vector_3_float::Vector_3_float()
{
	vec1_ = 1.0;
	vec2_ = 1.0;
	vec3_ = 1.0;
}

Vector_3_float::Vector_3_float(const float v1, const float v2, const float v3)
{
	vec1_ = v1;
	vec2_ = v2;
	vec3_ = v3;
}

Vector_3_float::Vector_3_float(const Vector_3_float& other)
{
	this->vec1_ = other.vec1_;
	this->vec2_ = other.vec2_;
	this->vec3_ = other.vec3_;
}

Vector_3_float& Vector_3_float::operator=(const Vector_3_float& other)
{
	this->vec1_ = other.vec1_;
	this->vec2_ = other.vec2_;
	this->vec3_ = other.vec3_;

	return *this;
}

Vector_3_float Vector_3_float::operator+(const Vector_3_float& other)
{
	Vector_3_float rez(this->vec1_ + other.vec1_, this->vec2_ + other.vec2_, this->vec3_ + other.vec3_);
	return rez;
}

Vector_3_float Vector_3_float::operator-(const Vector_3_float& other)
{
	Vector_3_float rez(this->vec1_ - other.vec1_, this->vec2_ - other.vec2_, this->vec3_ - other.vec3_);
	return rez;
}

float Vector_3_float::operator*(const Vector_3_float& other)
{
	float dot_product = this->vec1_ * other.vec1_ + this->vec2_ * other.vec2_ + this->vec3_ * other.vec3_;
	return dot_product;
}

Vector_3_float Vector_3_float::operator*(const float& value)
{

	Vector_3_float rez(value * this->get_v1(), value * this->get_v2(), value * this->get_v3());

	return rez;
}

float Vector_3_float::operator^(const Vector_3_float& other)
{
	float cos_angle = (this->get_v1() * other.vec1_ + this->get_v2() * other.vec2_ + this->get_v3() * other.vec3_) / (this->magnitude() * (sqrtf(other.vec1_*other.vec1_+other.vec2_*other.vec2_+other.vec3_*other.vec3_)));

	float angle = 180.0f * acos(cos_angle) / Pi;

	return angle;

}


Vector_3_float Vector_3_float::V_product(const Vector_3_float& other)
{
	Vector_3_float rez(vec2_ * other.vec3_ - vec3_ * other.vec2_, vec3_ * other.vec1_ - vec1_ * other.vec3_, vec1_ * other.vec2_ - vec2_ * other.vec1_);

	float tmp = vec2_ * other.vec3_ - vec3_ * other.vec2_;


	return rez;
}


float Vector_3_float::get_v1(void)
{
	return vec1_;
}

float Vector_3_float::get_v2(void)
{
	return vec2_;
}

float Vector_3_float::get_v3(void)
{
	return vec3_;
}


void Vector_3_float::normalize(void)
{
	float tmp = this->magnitude();

	vec1_ = vec1_ / tmp;
	vec2_ = vec2_ / tmp;
	vec3_ = vec3_ / tmp;
}

float Vector_3_float::magnitude(void)
{
	return sqrtf(vec1_*vec1_ + vec2_*vec2_+vec3_*vec3_);
}

Vector_3_float::~Vector_3_float()
{
}



float Volume(Vector_3_float, Vector_3_float, Vector_3_float);


float Volume(Vector_3_float U, Vector_3_float V, Vector_3_float W)
{
	float volume = U.get_v1() * (V.get_v2() * W.get_v3() - V.get_v3() * W.get_v2()) - U.get_v2() * (V.get_v1() * W.get_v3() - V.get_v3() * W.get_v1()) + U.get_v3() * (V.get_v1() * W.get_v2() - V.get_v2() * W.get_v1());

	return volume;
}

class Figure
{
private:


public:
	Figure();
	~Figure();

	virtual bool ray_intersect(const float, const float, const float, const float, const float, const float) = 0;
	virtual Vector_3_float ret_point(const float, const float, const float, const float, const float, const float) = 0;
	virtual Vector_3_float ret_normal(const float, const float, const float) = 0;

};

Figure::Figure()
{
}

Figure::~Figure()
{
}


class Sphere : public Figure
{
private:

	float x_;
	float y_;
	float z_;

	float R_;

public:
	Sphere();
	Sphere(const float, const float, const float, const float);
	~Sphere();

	void set_x(const float);
	void set_y(const float);
	void set_z(const float);
	void set_R(const float);

	bool ray_intersect(const float, const float, const float, const float, const float, const float) override;
	Vector_3_float ret_point(const float, const float, const float, const float, const float, const float) override;
	Vector_3_float ret_normal(const float, const float, const float) override;

	float get_x(void);
	float get_y(void);
	float get_z(void);
	float get_R(void);

};



Sphere::Sphere()
{
	x_ = 0.;
	y_ = 0.;
	z_ = 0.;

	R_ = 1.;


}
Sphere::Sphere(const float x, const float y, const float z, const float R)
{
	if (R < 0)
	{
		throw std::runtime_error("sphere.\n");
	}

	x_ = x;
	y_ = y;
	z_ = z;
	R_ = R;
}

float Sphere::get_x(void)
{
	return x_;
}
float Sphere::get_y(void)
{
	return y_;
}
float Sphere::get_z(void)
{
	return z_;
}
float Sphere::get_R(void)
{
	return R_;
}

bool Sphere::ray_intersect(const float origin_x, const float origin_y, const float origin_z, const float dir_x, const float dir_y, const float dir_z)  // Луч, выходящий из точки в данном направлении пересекает ли сферу?
{

	float a = (dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
	float b = (2.0 * ((dir_x) * (origin_x - x_) + dir_y * (origin_y - y_) + dir_z * (origin_z - z_)));
	float c = origin_x * origin_x + origin_y * origin_y + origin_z * origin_z + x_ * x_ + y_ * y_ + z_ * z_ - 2.0 * ((origin_x * x_) + origin_y * y_ + origin_z * z_) - (R_ * R_);

	float D = b * b - 4.0 * a * c;



	

	if (D < 0)
	{
		return false;
	}

	if (a == 0)
	{
		cout << "OOps" << endl;
	}

	float t = (-b - sqrtf(D)) / (2.0f * a);


	if (D >= 0)
	{
		if (t > 0)
		{
			return true;
		}

		if (t < 0)
		{
			return false;
		}

	}

	return false;

}


Vector_3_float Sphere::ret_point(const float origin_x, const float origin_y, const float origin_z, const float dir_x, const float dir_y, const float dir_z)
{
	float a = dir_x * dir_x + dir_y * dir_y + dir_z * dir_z;
	float b = 2.0f * (dir_x * (origin_x - x_) + dir_y * (origin_y - y_) + dir_z * (origin_z - z_));
	float c = origin_x * origin_x + origin_y * origin_y + origin_z * origin_z + x_ * x_ + y_ * y_ + z_ * z_ - 2.0f * (origin_x * x_ + origin_y * y_ + origin_z * z_) - R_ * R_;

	float D = b * b - 4.0f * a * c;

	float t = (-b - sqrtf(D)) / (2.0f * a);

	Vector_3_float pixel(origin_x + dir_x * t, origin_y + dir_y * t, origin_z + dir_z * t);

	return pixel;
}

Vector_3_float Sphere::ret_normal(const float x, const float y, const float z)
{
	Vector_3_float normal_surface(x - x_, y - y_, z - z_);

	return normal_surface;
}

void Sphere::set_x(const float x)
{
	x_ = x;
}
void Sphere::set_y(const float y)
{
	y_ = y;
}
void Sphere::set_z(const float z)
{
	z_ = z;
}
void Sphere::set_R(const float R)
{
	R_ = R;
}

Sphere::~Sphere()
{
}

class Box : public Figure
{
private:
	Vector_3_float begin;
	Vector_3_float end;

public:
	Box();
	Box(const float, const float, const float, const float, const float, const float);
	~Box();

	bool ray_intersect(const float, const float, const float, const float, const float, const float) override;
	Vector_3_float ret_point(const float, const float, const float, const float, const float, const float) override;
	Vector_3_float ret_normal(const float, const float, const float) override;


};

Box::Box()
{
	Vector_3_float A((float)0.0, (float)0.0, (float)0.0);
	Vector_3_float B((float)1.0, (float)1.0, (float)1.0);

}

Box::Box(const float x1, const float y1, const float z1, const float x2, const float y2, const float z2)
{

	if (x2 - x1 <= 0 || y2 - y1 <= 0 || z2 - z1 <= 0)
	{
		throw std::runtime_error("Failed to create box.\n");
	}


	Vector_3_float A(x1, y1, z1);
	Vector_3_float B(x2, y2, z2);

	begin = A;
	end = B;

}

Box::~Box()
{
}

bool Box::ray_intersect(const float origin_x, const float origin_y, const float origin_z, const float dir_x, const float dir_y, const float dir_z)
{



	if (dir_x !=0 && dir_y !=0 && dir_z !=0)
	{
		float t1 = (begin.get_v1() - origin_x) / dir_x;
		float Y1 = origin_y + dir_y * t1;
		float Z1 = origin_z + dir_z * t1;

		if (t1 > 0 && begin.get_v2() <= Y1 && Y1 <= end.get_v2() && begin.get_v3() <= Z1 && Z1 <= end.get_v3())
		{

			return true;
		}

		float t2 = (begin.get_v2() - origin_y) / dir_y;
		float X2 = origin_x + dir_x * t2;
		float Z2 = origin_z + dir_z * t2;

		if (t2 > 0 && begin.get_v1() <= X2 && X2 <= end.get_v1() && begin.get_v3() <= Z2 && Z2 <= end.get_v3())
		{
			return true;
		}

		float t3 = (begin.get_v3() - origin_z) / dir_z;
		float X3 = origin_x + dir_x * t3;
		float Y3 = origin_y + dir_y * t3;

		if (t3 > 0 && begin.get_v1() <= X3 && X3 <= end.get_v1() && begin.get_v2() <= Y3 && Y3 <= end.get_v2())
		{

			return true;
		}

		float t4 = (end.get_v1() - origin_x) / dir_x;
		float Y4 = origin_y + dir_y * t4;
		float Z4 = origin_z + dir_z * t4;

		if (t4 > 0 && begin.get_v2() <= Y4 && Y4 <= end.get_v2() && begin.get_v3() <= Z4 && Z4 <= end.get_v3())
		{

			return true;
		}

		float t5 = (end.get_v2() - origin_y) / dir_y;
		float X5 = origin_x + dir_x * t5;
		float Z5 = origin_z + dir_z * t5;

		if (t5 > 0 && begin.get_v1() <= X5 && X5 <= end.get_v1() && begin.get_v3() <= Z5 && Z5 <= end.get_v3())
		{

			return true;
		}

		float t6 = (end.get_v3() - origin_z) / dir_z;
		float X6 = origin_x + dir_x * t6;
		float Y6 = origin_y + dir_y * t6;

		if (t6 > 0 && begin.get_v1() <= X6 && X6 <= end.get_v1() && begin.get_v2() <= Y6 && Y6 <= end.get_v2())
		{

			return true;
		}

	}

	if (dir_x == 0 && dir_y != 0 && dir_z != 0)
	{
		if (begin.get_v1() <= origin_x && origin_x <= end.get_v1())
		{


			float t7 = (begin.get_v2() - origin_y) / dir_y;
			float Z7 = origin_z + dir_z * t7;

			if (t7 > 0 && begin.get_v3() <= Z7 && Z7 <= end.get_v3())
			{

				return true;
			}

			float t8 = (end.get_v2() - origin_y) / dir_y;
			float Z8 = origin_z + dir_z * t8;

			if (t8 > 0 && begin.get_v3() <= Z8 && Z8 <= end.get_v3())
			{

				return true;
			}

			float t9 = (begin.get_v3() - origin_z) / dir_z;
			float Y9 = origin_y + dir_y * t9;

			if (t9 > 0 && begin.get_v2() <= Y9 && Y9 <= end.get_v2())
			{

				return true;

			}

			float t10 = (end.get_v3() - origin_z) / dir_z;
			float Y10 = origin_y + dir_y * t10;

			if (t10 > 0 && begin.get_v2() <= Y10 && Y10 <= end.get_v2())
			{

				return true;

			}

		}
	}

	if (dir_y == 0 && dir_x != 0 && dir_z != 0)
	{
		if (begin.get_v2() <= origin_y && origin_y <= end.get_v2())
		{

			float t11 = (begin.get_v1() - origin_x) / dir_x;
			float Z11 = origin_z + dir_z * t11;

			if (t11 > 0 && begin.get_v3() <= Z11 && Z11 <= end.get_v3())
			{

				return true;
			}


			float t12 = (end.get_v1() - origin_x) / dir_x;
			float Z12 = origin_z + dir_z * t12;

			if (t12 > 0 && begin.get_v3() <= Z12 && Z12 <= end.get_v3())
			{

				return true;
			}


			float t13 = (begin.get_v3() - origin_z) / dir_z;
			float X13 = origin_x + dir_x * t13;

			if (t13 > 0 && begin.get_v1() <= X13 && X13 <= end.get_v1())
			{

				return true;
			}


			float t14 = (end.get_v3() - origin_z) / dir_z;
			float X14 = origin_x + dir_x * t14;

			if (t14 > 0 && begin.get_v1() <= X14 && X14 <= end.get_v1())
			{

				return true;
			}
		}
	}

	if (dir_z == 0 && dir_x != 0 && dir_y != 0)
	{

		if (begin.get_v3() <= origin_z && origin_z <= end.get_v3())
		{

			float t15 = (begin.get_v1() - origin_x) / dir_x;
			float Y15 = origin_y + dir_y * t15;

			if (t15 > 0 && begin.get_v2() <= Y15 && Y15 <= end.get_v2())
			{
				return true;
			}


			float t16 = (end.get_v1() - origin_x) / dir_x;
			float Y16 = origin_y + dir_y * t16;

			if (t16 > 0 && begin.get_v2() <= Y16 && Y16 <= end.get_v2())
			{

				return true;
			}


			float t17 = (begin.get_v2() - origin_y) / dir_y;
			float X17 = origin_x + dir_x * t17;

			if (t17 > 0 && begin.get_v1() <= X17 && X17 <= end.get_v1())
			{

				return true;
			}

			float t18 = (end.get_v2() - origin_y) / dir_y;
			float X18 = origin_x + dir_x * t18;

			if (t18 > 0 && begin.get_v1() <= X18 && X18 <= end.get_v1())
			{

				return true;
			}

		}

	}

	if (dir_x == 0 && dir_y == 0 && dir_z != 0)
	{
		if (begin.get_v1() <= origin_x && origin_x <= end.get_v1() && begin.get_v2() <= origin_y && origin_y <= end.get_v2())
		{

			float t19 = (begin.get_v3() - origin_z) / dir_z;

			if (t19 > 0)
			{

				return true;
			}

			float t20 = (end.get_v3() - origin_z) / dir_z;

			if (t20 > 0)
			{

				return true;
			}


		}
	}

	if (dir_y == 0 && dir_z == 0 && dir_x != 0)
	{
		if (begin.get_v2() <= origin_y && origin_y <= end.get_v2() && begin.get_v3() <= origin_z && origin_z <= end.get_v3())
		{
			float t21 = (begin.get_v1() - origin_x) / dir_x;

			if (t21 > 0)
			{

				return true;
			}

			float t22 = (end.get_v1() - origin_x) / dir_x;

			if (t22 > 0)
			{

				return true;
			}

		}

	}

	if (dir_x == 0 && dir_z == 0 && dir_y != 0)
	{
		if (begin.get_v1() <= origin_x && origin_x <= end.get_v1() && begin.get_v3() <= origin_z && origin_z <= end.get_v3())
		{
			float t23 = (begin.get_v2() - origin_y) / dir_y;

			if (t23 > 0)
			{

				return true;
			}

			float t24 = (end.get_v2() - origin_y) / dir_y;

			if (t24 > 0)
			{

				return true;
			}
		}
	}

	return false;

}

Vector_3_float Box::ret_point(const float origin_x, const float origin_y, const float origin_z, const float dir_x, const float dir_y, const float dir_z)
{


	vector<float> parameters;

	if (dir_x != 0 && dir_y != 0 && dir_z != 0)
	{
		float t1 = (begin.get_v1() - origin_x) / dir_x;
		float Y1 = origin_y + dir_y * t1;
		float Z1 = origin_z + dir_z * t1;

		if (t1 > 0 && begin.get_v2() <= Y1 && Y1 <= end.get_v2() && begin.get_v3() <= Z1 && Z1 <= end.get_v3())
		{
			parameters.push_back(t1);
		}

		float t2 = (begin.get_v2() - origin_y) / dir_y;
		float X2 = origin_x + dir_x * t2;
		float Z2 = origin_z + dir_z * t2;

		if (t2 > 0 && begin.get_v1() <= X2 && X2 <= end.get_v1() && begin.get_v3() <= Z2 && Z2 <= end.get_v3())
		{
			parameters.push_back(t2);
		}

		float t3 = (begin.get_v3() - origin_z) / dir_z;
		float X3 = origin_x + dir_x * t3;
		float Y3 = origin_y + dir_y * t3;

		if (t3 > 0 && begin.get_v1() <= X3 && X3 <= end.get_v1() && begin.get_v2() <= Y3 && Y3 <= end.get_v2())
		{
			parameters.push_back(t3);
		}

		float t4 = (end.get_v1() - origin_x) / dir_x;
		float Y4 = origin_y + dir_y * t4;
		float Z4 = origin_z + dir_z * t4;

		if (t4 > 0 && begin.get_v2() <= Y4 && Y4 <= end.get_v2() && begin.get_v3() <= Z4 && Z4 <= end.get_v3())
		{
			parameters.push_back(t4);
		}

		float t5 = (end.get_v2() - origin_y) / dir_y;
		float X5 = origin_x + dir_x * t5;
		float Z5 = origin_z + dir_z * t5;

		if (t5 > 0 && begin.get_v1() <= X5 && X5 <= end.get_v1() && begin.get_v3() <= Z5 && Z5 <= end.get_v3())
		{
			parameters.push_back(t5);
		}

		float t6 = (end.get_v3() - origin_z) / dir_z;
		float X6 = origin_x + dir_x * t6;
		float Y6 = origin_y + dir_y * t6;

		if (t6 > 0 && begin.get_v1() <= X6 && X6 <= end.get_v1() && begin.get_v2() <= Y6 && Y6 <= end.get_v2())
		{
			parameters.push_back(t6);
		}

	}

	if (dir_x == 0 && dir_y != 0 && dir_z != 0)
	{
		if (begin.get_v1() <= origin_x && origin_x <= end.get_v1())
		{


			float t7 = (begin.get_v2() - origin_y) / dir_y;
			float Z7 = origin_z + dir_z * t7;

			if (t7 > 0 && begin.get_v3() <= Z7 && Z7 <= end.get_v3())
			{
				parameters.push_back(t7);
			}

			float t8 = (end.get_v2() - origin_y) / dir_y;
			float Z8 = origin_z + dir_z * t8;

			if (t8 > 0 && begin.get_v3() <= Z8 && Z8 <= end.get_v3())
			{
				parameters.push_back(t8);
			}

			float t9 = (begin.get_v3() - origin_z) / dir_z;
			float Y9 = origin_y + dir_y * t9;

			if (t9 > 0 && begin.get_v2() <= Y9 && Y9 <= end.get_v2())
			{
				parameters.push_back(t9);

			}

			float t10 = (end.get_v3() - origin_z) / dir_z;
			float Y10 = origin_y + dir_y * t10;

			if (t10 > 0 && begin.get_v2() <= Y10 && Y10 <= end.get_v2())
			{
				parameters.push_back(t10);

			}

		}
	}

	if (dir_y == 0 && dir_x != 0 && dir_z != 0)
	{
		if (begin.get_v2() <= origin_y && origin_y <= end.get_v2())
		{

			float t11 = (begin.get_v1() - origin_x) / dir_x;
			float Z11 = origin_z + dir_z * t11;

			if (t11 > 0 && begin.get_v3() <= Z11 && Z11 <= end.get_v3())
			{
				parameters.push_back(t11);
			}


			float t12 = (end.get_v1() - origin_x) / dir_x;
			float Z12 = origin_z + dir_z * t12;

			if (t12 > 0 && begin.get_v3() <= Z12 && Z12 <= end.get_v3())
			{
				parameters.push_back(t12);
			}


			float t13 = (begin.get_v3() - origin_z) / dir_z;
			float X13 = origin_x + dir_x * t13;

			if (t13 > 0 && begin.get_v1() <= X13 && X13 <= end.get_v1())
			{
				parameters.push_back(t13);
			}


			float t14 = (end.get_v3() - origin_z) / dir_z;
			float X14 = origin_x + dir_x * t14;

			if (t14 > 0 && begin.get_v1() <= X14 && X14 <= end.get_v1())
			{
				parameters.push_back(t14);
			}


		}


	}

	if (dir_z == 0 && dir_x != 0 && dir_y != 0)
	{

		if (begin.get_v3() <= origin_z && origin_z <= end.get_v3())
		{

			float t15 = (begin.get_v1() - origin_x) / dir_x;
			float Y15 = origin_y + dir_y * t15;

			if (t15 > 0 && begin.get_v2() <= Y15 && Y15 <= end.get_v2())
			{
				parameters.push_back(t15);
			}


			float t16 = (end.get_v1() - origin_x) / dir_x;
			float Y16 = origin_y + dir_y * t16;

			if (t16 > 0 && begin.get_v2() <= Y16 && Y16 <= end.get_v2())
			{
				parameters.push_back(t16);
			}


			float t17 = (begin.get_v2() - origin_y) / dir_y;
			float X17 = origin_x + dir_x * t17;

			if (t17 > 0 && begin.get_v1() <= X17 && X17 <= end.get_v1())
			{
				parameters.push_back(t17);
			}

			float t18 = (end.get_v2() - origin_y) / dir_y;
			float X18 = origin_x + dir_x * t18;

			if (t18 > 0 && begin.get_v1() <= X18 && X18 <= end.get_v1())
			{
				parameters.push_back(t18);
			}

		}

	}

	if (dir_x == 0 && dir_y == 0 && dir_z != 0)
	{
		if (begin.get_v1() <= origin_x && origin_x <= end.get_v1() && begin.get_v2() <= origin_y && origin_y <= end.get_v2())
		{

			float t19 = (begin.get_v3() - origin_z) / dir_z;

			if (t19 > 0)
			{
				parameters.push_back(t19);
			}

			float t20 = (end.get_v3() - origin_z) / dir_z;

			if (t20 > 0)
			{
				parameters.push_back(t20);
			}


		}
	}

	if (dir_y == 0 && dir_z == 0 && dir_x != 0)
	{
		if (begin.get_v2() <= origin_y && origin_y <= end.get_v2() && begin.get_v3() <= origin_z && origin_z <= end.get_v3())
		{
			float t21 = (begin.get_v1() - origin_x) / dir_x;

			if (t21 > 0)
			{
				parameters.push_back(t21);
			}

			float t22 = (end.get_v1() - origin_x) / dir_x;

			if (t22 > 0)
			{
				parameters.push_back(t22);
			}

		}

	}

	if (dir_x == 0 && dir_z == 0 && dir_y != 0)
	{
		if (begin.get_v1() <= origin_x && origin_x <= end.get_v1() && begin.get_v3() <= origin_z && origin_z <= end.get_v3())
		{
			float t23 = (begin.get_v2() - origin_y) / dir_y;

			if (t23 > 0)
			{
				parameters.push_back(t23);
			}

			float t24 = (end.get_v2() - origin_y) / dir_y;

			if (t24 > 0)
			{
				parameters.push_back(t24);
			}
		}
	}


	float min = parameters[0];

	for (auto i : parameters)
	{
		if (abs(i -min)<BIAS)
		{
			min = i;
		}

	}

	Vector_3_float origin(origin_x, origin_y, origin_z);
	Vector_3_float direction(dir_x, dir_y, dir_z);
	Vector_3_float point = origin + direction * min;

	return point;


}

Vector_3_float Box::ret_normal(const float x, const float y, const float z) // Точка Уже на поверхности, но на какой грани7???????
{
	Vector_3_float n1(1.0, 0.0, 0.0);
	Vector_3_float n2(0.0, 1.0, 0.0);
	Vector_3_float n3(0.0, 0.0, 1.0);
	Vector_3_float n4(-1.0, 0.0, 0.0);
	Vector_3_float n5(0.0, -1.0, 0.0);
	Vector_3_float n6(0.0, 0.0, -1.0);

	Vector_3_float diagonal = end - begin;


	if (abs(x - begin.get_v1())<epsi)
	{
		return n4;
	}

	if ((x - end.get_v1())<epsi)
	{
		return n1;
	}


	if (abs(y -begin.get_v2())< epsi)
	{
		return n5;
	}

	if (abs(y - end.get_v2())< epsi)
	{
		return n2;
	}

	if (abs(z - begin.get_v3())< epsi)
	{
		return n6;
	}

	if (abs(z - end.get_v3())< epsi)
	{
		return n3;
	}

}

class Tetrahedron : public Figure
{
private:

	Vector_3_float first_vertex;
	Vector_3_float second_vertex;
	Vector_3_float third_vertex;
	Vector_3_float fourth_vertex;

public:

	Tetrahedron(const float, const float, const float, const float, const float, const float, const float, const float, const float, const float, const float, const float);
	~Tetrahedron();

	bool ray_intersect(const float, const float, const float, const float, const float, const float) override;
	Vector_3_float ret_point(const float, const float, const float, const float, const float, const float) override;
	Vector_3_float ret_normal(const float, const float, const float) override;

	friend bool Triangle_intersection(const float, const float, const float, const float, const float, const float, Vector_3_float, Vector_3_float, Vector_3_float);
	friend float return_param(const float, const float, const float, const float, const float, const float, Vector_3_float, Vector_3_float, Vector_3_float);
	friend bool check_point(const float, const float, const float, Vector_3_float, Vector_3_float, Vector_3_float);


};



bool Triangle_intersection(const float origin_x, const float origin_y, const float origin_z, const float dir_x, const float dir_y, const float dir_z, Vector_3_float v0_, Vector_3_float v1_, Vector_3_float v2_)
{

	Vector_3_float v0 = v0_;
	Vector_3_float v1 = v1_;
	Vector_3_float v2 = v2_;
	Vector_3_float direction(dir_x, dir_y, dir_z);
	Vector_3_float origin(origin_x, origin_y, origin_z);

	Vector_3_float v0v1 = v1 - v0;
	Vector_3_float v0v2 = v2 - v0;
	Vector_3_float pvec = direction.V_product(v0v2);


	float det = v0v1 * pvec;

	if (det < BIAS)
	{
		return false;
	}

	if (std::fabs(det) < BIAS)
	{
		return false;
	}

	float invDet = 1 / det;

	Vector_3_float tvec = origin - v0;
	float u = tvec * pvec * invDet;

	if (u < 0 || u > 1)
	{
		return false;

	}

	Vector_3_float qvec = tvec.V_product(v0v1);

	float v = (direction * qvec) * invDet;

	if (v < 0 || u + v > 1)
	{
		return false;
	}

	float t = (v0v2 * qvec) * invDet;

	if (t > BIAS)
	{
		return true;
	}

	return false;
}


float return_param(const float origin_x, const float origin_y, const float origin_z, const float dir_x, const float dir_y, const float dir_z, Vector_3_float v0_, Vector_3_float v1_, Vector_3_float v2_)
{
	Vector_3_float v0 = v0_;
	Vector_3_float v1 = v1_;
	Vector_3_float v2 = v2_;
	Vector_3_float direction(dir_x, dir_y, dir_z);
	Vector_3_float origin(origin_x, origin_y, origin_z);

	Vector_3_float v0v1 = v1 - v0;
	Vector_3_float v0v2 = v2 - v0;
	Vector_3_float pvec = direction.V_product(v0v2);

	float det = v0v1 * pvec;
	float invDet = 1 / det;

	Vector_3_float tvec = origin - v0;
	Vector_3_float qvec = tvec.V_product(v0v1);

	float t = (v0v2 * qvec) * invDet;

	return t;
}

bool check_point(const float x_, const float y_, const float z_, Vector_3_float A_, Vector_3_float B_, Vector_3_float C_)
{
	Vector_3_float p(x_, y_, z_);
	Vector_3_float a = A_;
	Vector_3_float b = B_;
	Vector_3_float c = C_;

	a = a - p;
	b = b - p;
	c = c - p;

	Vector_3_float u = b.V_product(c);
	Vector_3_float v = c.V_product(a);
	Vector_3_float w = a.V_product(b);

	if (u.magnitude() == 0.0 || v.magnitude() == 0.0 || w.magnitude() == 0.0)
	{
		return true;
	}

	if ((u * v) < 0.0)
	{
		return false;
	}

	if ((u * w) < 0.0)
	{
		return false;
	}

	return true;
}

Tetrahedron::Tetrahedron(const float x1, const float y1, const float z1, const float x2, const float y2, const float z2, const float x3, const float y3, const float z3, const float x4, const float y4, const float z4)
{

	Vector_3_float A(x1, y1, z1);
	Vector_3_float B(x2, y2, z2);
	Vector_3_float C(x3, y3, z3);
	Vector_3_float D(x4, y4, z4);

	first_vertex = A;
	second_vertex = B;
	third_vertex = C;
	fourth_vertex = D;

}


bool Tetrahedron::ray_intersect(const float origin_x, const float origin_y, const float origin_z, const float dir_x, const float dir_y, const float dir_z)
{


	if (Triangle_intersection(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, first_vertex, second_vertex, third_vertex))
	{
		return true;
	}

	if (Triangle_intersection(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, third_vertex, second_vertex, fourth_vertex))
	{
		return true;
	}


	if (Triangle_intersection(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, fourth_vertex, second_vertex, first_vertex))
	{
		return true;
	}

	if (Triangle_intersection(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, first_vertex, third_vertex, fourth_vertex))
	{
		return true;
	}

	return false;


}

Vector_3_float Tetrahedron::ret_point(const float origin_x, const float origin_y, const float origin_z, const float dir_x, const float dir_y, const float dir_z)
{

	Vector_3_float origin(origin_x, origin_y, origin_z);
	Vector_3_float direction(dir_x, dir_y, dir_z);
	Vector_3_float point;

	vector<float> parameters;

	if (Triangle_intersection(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, first_vertex, second_vertex, third_vertex))
	{
		parameters.push_back(return_param(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, first_vertex, second_vertex, third_vertex));
	}

	if (Triangle_intersection(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, third_vertex, second_vertex, fourth_vertex))
	{
		parameters.push_back(return_param(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, third_vertex, second_vertex, fourth_vertex));
	}


	if (Triangle_intersection(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, fourth_vertex, second_vertex, first_vertex))
	{
		parameters.push_back(return_param(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, fourth_vertex, second_vertex, first_vertex));
	}

	if (Triangle_intersection(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, first_vertex, third_vertex, fourth_vertex))
	{
		parameters.push_back(return_param(origin_x, origin_y, origin_z, dir_x, dir_y, dir_z, first_vertex, third_vertex, fourth_vertex));
	}


	if (parameters.capacity() == 0)
	{
		throw std::runtime_error("CERROR\n");
	}

	float min = parameters[0];

	for (auto i : parameters)
	{
		if (i < min)
		{
			min = i;
		}

	}

	point = origin + direction * min;

	return point;


}


Vector_3_float Tetrahedron::ret_normal(const float x, const float y, const float z) // Точка, которая точно лежит на грани тетраэдера, но на какой??
{

	if (check_point(x, y, z, first_vertex, second_vertex, third_vertex))
	{
		Vector_3_float u = second_vertex - first_vertex;
		Vector_3_float v = third_vertex - first_vertex;

		Vector_3_float D = fourth_vertex - first_vertex;



		Vector_3_float Normal = u.V_product(v);

		if ((Normal * D) > 0.0)
		{
			Normal = Normal * -1.0;
			return Normal;
		}

		return Normal;

	}

	if (check_point(x, y, z, third_vertex, second_vertex, fourth_vertex))
	{

		Vector_3_float u = second_vertex - third_vertex;
		Vector_3_float v = fourth_vertex - third_vertex;

		Vector_3_float D = first_vertex - third_vertex;



		Vector_3_float Normal = u.V_product(v);

		if ((Normal * D) > 0.0)
		{
			Normal = Normal * -1.0;
			return Normal;
		}

		return Normal;
	}


	if (check_point(x, y, z, fourth_vertex, second_vertex, first_vertex))
	{
		Vector_3_float u = second_vertex - fourth_vertex;
		Vector_3_float v = first_vertex - fourth_vertex;

		Vector_3_float D = third_vertex - fourth_vertex;



		Vector_3_float Normal = u.V_product(v);

		if ((Normal * D) > 0.0)
		{
			Normal = Normal * -1.0;
			return Normal;
		}

		return Normal;
	}

	if (check_point(x, y, z, first_vertex, third_vertex, fourth_vertex))
	{

		Vector_3_float u = third_vertex - first_vertex;
		Vector_3_float v = fourth_vertex - first_vertex;

		Vector_3_float D = second_vertex - first_vertex;



		Vector_3_float Normal = u.V_product(v);

		if ((Normal * D) > 0.0)
		{
			Normal = Normal * -1.0;
			return Normal;
		}

		return Normal;
	}




}

Tetrahedron::~Tetrahedron()
{
}


class Spectator
{
private:
	float x_;
	float y_;
	float z_;

	float dist_spec_screen_;
	float dist_spec_scene_;
	float angle_of_view_;

public:
	Spectator();

	~Spectator();

	void set_x(const float);
	void set_y(const float);
	void set_z(const float);

	void set_dist_spec_screen(const float);
	void set_dist_spec_scene(const float);
	void set_angle_of_view(const float);

	float get_x(void);
	float get_y(void);
	float get_z(void);


	float get_dist_spec_screen(void);
	float get_dist_spec_scene(void);
	float get_angle_of_view(void);

};

Spectator::Spectator()
{
	x_ = 0.;
	y_ = 0.;
	z_ = 0.;

	dist_spec_screen_ = 5.;
	dist_spec_scene_ = 10.;
	angle_of_view_ = 70.;
}

void Spectator::set_x(const float x)
{
	x_ = x;
}
void Spectator::set_y(const float y)
{
	y_ = y;
}
void Spectator::set_z(const float z)
{
	z_ = z;
}

void Spectator::set_dist_spec_screen(const float v)
{
	dist_spec_screen_ = v;
}
void Spectator::set_dist_spec_scene(const float v)
{
	dist_spec_scene_ = v;
}
void Spectator::set_angle_of_view(const float v)
{
	angle_of_view_ = v;
}

float Spectator::get_x(void)
{
	return x_;
}
float Spectator::get_y(void)
{
	return y_;
}
float Spectator::get_z(void)
{
	return z_;
}

float Spectator::get_dist_spec_screen(void)
{
	return dist_spec_screen_;
}
float Spectator::get_dist_spec_scene(void)
{
	return dist_spec_scene_;
}
float Spectator::get_angle_of_view(void)
{
	return angle_of_view_;
}

Spectator::~Spectator()
{
}

class Screen
{
private:

	float x_;
	float y_;
	float z_;


	float normal[3];
	float up[3];
	float tangent[3];

	int width_;
	int height_;


public:
	Screen();
	~Screen();

	void set_x(const float);
	void set_y(const float);
	void set_z(const float);

	void set_norm_x(const float);
	void set_norm_y(const float);
	void set_norm_z(const float);

	void set_up_x(const float);
	void set_up_y(const float);
	void set_up_z(const float);

	void set_width(const int);
	void set_height(const int);



	float get_x(void);
	float get_y(void);
	float get_z(void);


	float get_norm_x(void);
	float get_norm_y(void);
	float get_norm_z(void);

	float get_up_x(void);
	float get_up_y(void);
	float get_up_z(void);

	int get_width(void);
	int get_height(void);

	void tangent_(void);

	float get_tangent_x(void);
	float get_tangent_y(void);
	float get_tangent_z(void);

};



Screen::Screen()
{

	x_ = 0.;
	y_ = 0.;
	z_ = 0.;

	normal[0] = 1.;
	normal[1] = 0.;
	normal[2] = 0.;

	up[0] = 0.;
	up[1] = 0.;
	up[2] = 1.;

	tangent[0] = 0.;
	tangent[1] = 1.;
	tangent[2] = 0.;

	width_ = 500;
	height_ = 500;
}

int Screen::get_width(void)
{
	return width_;
}
int Screen::get_height(void)
{
	return height_;
}

void Screen::tangent_(void)
{

	tangent[0] = up[1] * normal[2] - up[2] * normal[1];
	tangent[1] = up[2] * normal[0] - up[0] * normal[2];
	tangent[2] = up[0] * normal[1] - up[1] * normal[0];

}

float Screen::get_norm_x(void)
{
	return normal[0];
}
float Screen::get_norm_y(void)
{
	return normal[1];
}
float Screen::get_norm_z(void)
{
	return normal[2];
}

float Screen::get_up_x(void)
{
	return up[0];
}
float Screen::get_up_y(void)
{
	return up[1];
}
float Screen::get_up_z(void)
{
	return up[2];
}

float Screen::get_tangent_x(void)
{
	return tangent[0];
}
float Screen::get_tangent_y(void)
{
	return tangent[1];
}
float Screen::get_tangent_z(void)
{
	return tangent[2];
}

float Screen::get_x(void)
{
	return x_;
}
float Screen::get_y(void)
{
	return y_;
}
float Screen::get_z(void)
{
	return z_;
}

void Screen::set_norm_x(const float x)
{
	normal[0] = x;
}
void Screen::set_norm_y(const float y)
{
	normal[1] = y;

}
void Screen::set_norm_z(const float z)
{
	normal[2] = z;
}

void Screen::set_up_x(const float x)
{
	up[0] = x;

}
void Screen::set_up_y(const float y)
{
	up[1] = y;
}
void Screen::set_up_z(const float z)
{
	up[2] = z;
}

void Screen::set_x(const float x)
{
	x_ = x;
}
void Screen::set_y(const float y)
{
	y_ = y;
}
void Screen::set_z(const float z)
{
	z_ = z;
}

void Screen::set_width(const int v)
{

	width_ = v;
}
void Screen::set_height(const int v)
{
	height_ = v;
}

Screen::~Screen()
{
}

class Light
{
private:

	float x_;
	float y_;
	float z_;

public:
	Light();
	~Light();

	void set_x(const float);
	void set_y(const float);
	void set_z(const float);


	float get_x();
	float get_y();
	float get_z();

};

Light::Light()
{
	x_ = 0.;
	y_ = 0.;
	z_ = 0.;
}

void Light::set_x(const float x)
{
	x_ = x;
}
void Light::set_y(const float y)
{
	y_ = y;
}
void Light::set_z(const float z)
{
	z_ = z;
}

float Light::get_x()
{
	return x_;
}

float Light::get_y()
{
	return y_;
}

float Light::get_z()
{
	return z_;
}


Light::~Light()
{
}



int main()
{

	try
	{
		Spectator camera;
		Screen screen;
		Light lamp;

		vector<Figure*> shapes;

		string file_name;
		string cur_string;

		ifstream file;

		std::cout << "Enter a file name:" << endl;
		cin >> file_name;

		file.open(file_name);

		if (!file.is_open())
		{
			throw std::runtime_error("Failed to open file.\n");
		}


		while (getline(file, cur_string))
		{
			stringstream str_stream;

			str_stream << cur_string;

			string mark;

			str_stream >> mark;

			if (mark == "cam")
			{
				float v1;
				float v2;
				float v3;

				str_stream >> v1;
				str_stream >> v2;
				str_stream >> v3;

				camera.set_x(v1);
				camera.set_y(v2);
				camera.set_z(v3);
			}

			if (mark == "normal")
			{
				float v1;
				float v2;
				float v3;

				str_stream >> v1;
				str_stream >> v2;
				str_stream >> v3;

				screen.set_norm_x(v1);
				screen.set_norm_y(v2);
				screen.set_norm_z(v3);

				Vector_3_float tmp(v1, v2, v3);

				if (tmp.magnitude() == 0)
				{
					throw std::runtime_error("Error.\n");
				}

			}

			if (mark == "up")
			{
				float v1;
				float v2;
				float v3;

				str_stream >> v1;
				str_stream >> v2;
				str_stream >> v3;

				screen.set_up_x(v1);
				screen.set_up_y(v2);
				screen.set_up_z(v3);

			}

			if (mark == "screen")
			{
				float v1;
				str_stream >> v1;
				camera.set_dist_spec_screen(v1);

				if (v1 <= 0.0)
				{
					throw std::runtime_error("Error.\n");
				}
			}

			if (mark == "limit")
			{
				float v1;
				str_stream >> v1;
				camera.set_dist_spec_scene(v1);

			}

			if (mark == "alpha")
			{
				float v1;
				str_stream >> v1;
				camera.set_angle_of_view(v1);

			}

			if (mark == "width")
			{
				int width;
				str_stream >> width;
				screen.set_width(width);
			}

			if (mark == "height")
			{
				int height;
				str_stream >> height;
				screen.set_height(height);
			}

			if (mark == "light")
			{
				float v1;
				float v2;
				float v3;

				str_stream >> v1;
				str_stream >> v2;
				str_stream >> v3;

				lamp.set_x(v1);
				lamp.set_y(v2);
				lamp.set_z(v3);
			}

			if (mark == "sphere")
			{
				float v1;
				float v2;
				float v3;
				float r;

				str_stream >> v1;
				str_stream >> v2;
				str_stream >> v3;
				str_stream >> r;

				shapes.push_back(new Sphere(v1, v2, v3, r));

			}

			if (mark == "box")
			{
				float v1;
				float v2;
				float v3;

				float w1;
				float w2;
				float w3;

				str_stream >> v1;
				str_stream >> v2;
				str_stream >> v3;

				str_stream >> w1;
				str_stream >> w2;
				str_stream >> w3;

				shapes.push_back(new Box(v1, v2, v3, w1, w2, w3));

			}

			if (mark == "tetra")
			{
				float a1;
				float a2;
				float a3;

				float b1;
				float b2;
				float b3;

				float c1;
				float c2;
				float c3;

				float d1;
				float d2;
				float d3;

				str_stream >> a1;
				str_stream >> a2;
				str_stream >> a3;

				str_stream >> b1;
				str_stream >> b2;
				str_stream >> b3;

				str_stream >> c1;
				str_stream >> c2;
				str_stream >> c3;

				str_stream >> d1;
				str_stream >> d2;
				str_stream >> d3;

				shapes.push_back(new Tetrahedron(a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3));
			}

		}

		if (camera.get_dist_spec_scene() <= camera.get_dist_spec_screen())
		{
			throw std::runtime_error("Error1.\n");
		}

		if (camera.get_angle_of_view() <0.0 || camera.get_angle_of_view() >=180.0)
		{
			throw std::runtime_error("Error2.\n");
		}

		if (screen.get_height() <= 0 || screen.get_width() <=0 )
		{
			throw std::runtime_error("Error3.\n");
		}



		screen.tangent_();


		Vector_3_float normal(screen.get_norm_x(), screen.get_norm_y(), screen.get_norm_z());
		Vector_3_float up(screen.get_up_x(), screen.get_up_y(), screen.get_up_z());
		Vector_3_float tangent(screen.get_tangent_x(), screen.get_tangent_y(), screen.get_tangent_z());

		normal.normalize();
		up.normalize();
		tangent.normalize();

		screen.set_x(camera.get_x() + (camera.get_dist_spec_screen() * screen.get_norm_x()) / normal.magnitude()); /// Координаты центра Экрана
		screen.set_y(camera.get_y() + (camera.get_dist_spec_screen() * screen.get_norm_y()) / normal.magnitude());
		screen.set_z(camera.get_z() + (camera.get_dist_spec_screen() * screen.get_norm_z()) / normal.magnitude());

		Vector_3_float center(screen.get_x(), screen.get_y(), screen.get_z());
		Vector_3_float viewer(camera.get_x(), camera.get_y(), camera.get_z());

		if ((up ^ normal) != 90.0)
		{
			throw std::runtime_error("Wrong angle betwen normal and up\n");
		}

		CImg<unsigned char> image(screen.get_width(), screen.get_height(), 1, 3, 0);




		for (auto L : shapes)
		{


			int randomColor[3];

			int g = rand() % 31 + 1;

			randomColor[0] = 36;
			randomColor[1] = 36;
			randomColor[2] = 36;

			//#pragma omp parallel for

			for (int j = 0; j < screen.get_height(); j++)
			{
				for (int i = 0; i < screen.get_width(); i++)
				{

					
					int Cx = 0;
					int Cy = 0;

					if (screen.get_width() % 2 == 0)
					{
						Cx = i - screen.get_width() / 2;
					}

					if (screen.get_width() % 2 != 0)
					{
						Cx = i - (screen.get_width() - 1) / 2;
					}

					if (screen.get_height() % 2 == 0)
					{
						Cy = (screen.get_height() / 2) - j;
					}

					if (screen.get_height() % 2 != 0)
					{
						Cy = ((screen.get_height() - 1) / 2) - j;
					}

					/*	cout << Cx << endl << Cy << endl;*/


						//// Cx, Cy - координаты точки в плоскости экрана, через которую надо пускать луч, которую надо красить

					float Vx = (float)Cx / (float)screen.get_width(); ///// Объемные координаты точки в плоскости, через которую надо пускать луч, но нету Vz??
					float Vy = (float)Cy / (float)screen.get_width();

					Vector_3_float unit_up(up.get_v1(), up.get_v2(), up.get_v3());
					Vector_3_float unit_tangent(tangent.get_v1(), tangent.get_v2(), tangent.get_v3());
					Vector_3_float current;


					unit_up.normalize();
					unit_tangent.normalize();

					unit_tangent = unit_tangent * -1.0; ////////


					current = unit_up * Vy + unit_tangent * Vx;

					Vector_3_float dir = center + current - viewer;

					//cout << "i am here" << endl;


					float origin_x = viewer.get_v1(); /// Точка из которой должен идти луч. dir - Луч
					float origin_y = viewer.get_v2();
					float origin_z = viewer.get_v3();

					float dir1 = dir.get_v1();
					float dir2 = dir.get_v2();
					float dir3 = dir.get_v3();

					Vector_3_float forward = center - viewer;
					Vector_3_float variation = forward + unit_up * Vy;



					if ((forward ^ variation) > (camera.get_angle_of_view() / 2.0))
					{
						continue;
					}

					/*	if (!(L->ray_intersect(origin_x, origin_y, origin_z, dir.get_v1(), dir.get_v2(), dir.get_v3())))
					{
						cout << "what" << endl;
					}*/



					if (L->ray_intersect(origin_x, origin_y, origin_z, dir.get_v1(), dir.get_v2(), dir.get_v3()))
					{



						/*image.draw_point(i, j, randomColor);
						continue;*/

						Vector_3_float point = L->ret_point(origin_x, origin_y, origin_z, dir.get_v1(), dir.get_v2(), dir.get_v3()); // Точка точно на поверхности
						Vector_3_float normal_surface = L->ret_normal(point.get_v1(), point.get_v2(), point.get_v3()); //  Вызываю нормаль с точки, которая уже на поверхности
						Vector_3_float light(lamp.get_x(), lamp.get_y(), lamp.get_z()); // Вектор до света
						Vector_3_float point_to_lamp = light - point; // Вектор от точки на поверхности до лампы
						Vector_3_float reflection;
						Vector_3_float lamp_to_point = point - light;
						Vector_3_float point_to_camera = viewer - point;


						normal_surface.normalize();
						point_to_lamp.normalize();


						unsigned char curcolor[3];
						curcolor[0] = randomColor[0];
						curcolor[1] = randomColor[1];
						curcolor[2] = randomColor[2];


						float light_intense = point_to_lamp * normal_surface;

						if (light_intense <= 0)
						{


							curcolor[0] = (char)0;
							curcolor[1] = (char)0;
							curcolor[2] = (char)0;

							image.draw_point(i, j, curcolor);

							continue;
						}

						reflection = lamp_to_point - normal_surface * 2 * (lamp_to_point * normal_surface);

						reflection.normalize();



						float phong_intense;
						float mirror_coeff = 0.82f;
						float diffuse_coeff = 0.46f;

						float shine = 0.75;

						float R = (float)randomColor[0] * light_intense;
						float G = (float)randomColor[1] * light_intense;
						float B = (float)randomColor[2] * light_intense;

						curcolor[0] = (char)round(R);
						curcolor[1] = (char)round(G);
						curcolor[2] = (char)round(B);

						phong_intense = -diffuse_coeff * (normal_surface * lamp_to_point) + mirror_coeff * pow((normal_surface * reflection), shine);


						curcolor[0] = curcolor[1] = curcolor[2] = curcolor[0] * light_intense + phong_intense * light_intense;




						image.draw_point(i, j, curcolor);

					}
					
						
					


					




				}
			}
		}

		image.display();

		return 0;
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << endl;
		return 1;
	}

}