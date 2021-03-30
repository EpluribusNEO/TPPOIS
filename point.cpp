#include "point.h"

double CPoint::distance(const CPoint& P) const{
	return (double)sqrt((x - P.x)*(x - P.x) + (y - P.y)*(y - P.y));
}

CPoint operator + (const CPoint& C1, const CPoint& C2){
	CPoint R(C1.x + C2.x, C1.y + C2.y);
	return R;
}

CPoint operator -(const CPoint& C1, const CPoint& C2){
	CPoint R(C1.x - C2.x, C1.y - C2.y);
	return R;
}

CPoint operator * (const CPoint& C1, const double& m){
	CPoint R(C1.x*m, C1.y * m);
	return R;
}

CPoint operator / (const CPoint& C1, const double& d){
	CPoint R(C1.x / d, C1.y / d);
	return R;
}

std::ostream& operator << (std::ostream& out, const CPoint& P){
	out << P.x << " " << P.y << " " << P.z;
	return out;
}

// TODO
