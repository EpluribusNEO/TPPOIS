#pragma once

#ifndef POINT_H
#define POINT_H
#include <cmath>
#include <iostream>

struct CPoint {
	double x, y, z;
	CPoint(double _x = 0, double _y = 0, double _z = 0):x(_x),y(_y),z(_z) {};
	double distance(const CPoint& P) const;
	double dx(const CPoint& P) const { return fabs(x - P.x); }
	double dy(const CPoint& P) const { return fabs(y - P.y); }
};

CPoint operator + (const CPoint& C1, const CPoint& C2);
CPoint operator -(const CPoint& C1, const CPoint& C2);
CPoint operator / (const CPoint& C1, const double& d);
CPoint operator * (const CPoint& C1, const double& m);

std::ostream& operator << (std::ostream& out, const CPoint& P);



#endif // POINT_H
