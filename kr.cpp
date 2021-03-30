
#include <iostream>
#include <list>
#include <vector>
#include <set>
#include <fstream>
#include <string>
#include "point.h"
#include "/Lib-test/Matrix.h"
#include "/Lib-test/IsoLines.h"

#define M_PI 3.14159265358979323846

using namespace std;

bool TestTask = 1; // 1 -тестоваязадача, 0 -карта

int N = 100, M = 100; // разбиение исходной области Omega-Квадрат[0, Lx]x[0, Ly]

double x_min(0), x_max(0), y_min(0), y_max(0), z_min(0), z_max(0);
double eps = 0.0001;

//void LoadIsoLineFromFile(char *filename, IsoLines& line);
//void TecPlotV(char *filename, const IsoLines& lines);
void TecPlotU(char *filename, const Matrix& U);
CPoint SearchPointInGrid(constCPoint& P, int& ii, int& jj);
void ApproximateIsoLineInGrid(IsoLines& lines, Matrix& U);
Matrix Init(int n, int m);
void SegmentDecomposition(constCPoint& V1, constCPoint& V2, ListPoints& points);
void Load(string filename, IsoLines& lines);//Загрузка
void Search_min_max(const IsoLines& lines);
void LineSplineX(Matrix& U);
void LineSplineY(Matrix& U);

int main()
{
	setlocale(0, "");
	IsoLines lines;
	if (TestTask) LoadIsoLineFromFile("Example.txt", lines);
	else Load("Borus_full.txt", lines);
	cout << "Изолиниизагружены!\n";

	Search_min_max(lines);
	TecPlotV("Isolines.dat", lines);
	Matrix U = Init(N, M);

	ApproximateIsoLineInGrid(lines, U);
	LineSplineY(U);
	LineSplineX(U);

	x_max = 1.01*(x_max - x_min) * 64000; x_min = -0.01*x_max;
	y_max = 1.01*(y_max - y_min) * 111000; y_min = -0.01*y_max;
	cout << "Загрузка данных для TecPlot...\n";
	TecPlotU("Surface.dat", U);

	system("pause");
	return 0;
}



// Инициализация двумерногомассива


Matrix Init(int n, int m) {
	Matrix M;
	M.resize(n + 1);
	for (inti = 0; i < n + 1; i++) {
		M[i].resize(m + 1);
		for (int j = 0; j < m + 1; j++)
			M[i][j] = z_min;
	}
	return M;
}

void Load(string filename, IsoLines& lines)//Загрузка
{
	fstream Loading(filename.c_str());
	string str, str2;
	char symbol;
	double h, x, y, z;

	while (!Loading.eof())
	{
		do
		{
			getline(Loading, str, '\n');
		} while (str != "[POLYLINE]");
		getline(Loading, str, '\n'); // Type
		getline(Loading, str, '='); // Label
		Loading >> h;
		getline(Loading, str, '='); // Data0
		Loading.get();
		lines.push_back(make_pair(list<CPoint>(), h));

		while (1)
		{
			Loading >> y >> symbol >> x;
			do
			{
				Loading >> symbol;

			} while (!isdigit(Loading.peek()) && symbol != '[');
			lines.back().first.push_back(CPoint(x, y));
			if (symbol == '[' || symbol == 'E') break;


		}
		getline(Loading, str, '\n');
		Loading.get();

	}
}


void LoadIsoLineFromFile(char *filename, IsoLines&lines)
{
	ifstream fin(filename);
	double h, x, y;
	int count;

	while (!fin.eof())
	{
		fin >> h; // высота изолинн
		fin >> count; // количество вершин изолинн
		// вставляем в список новую изолинию

		lines.push_back(make_pair(list<CPoint>(), h));

		for (int i = 0; i < count; i++)
		{
			fin >> x >> y;
			lines.back().first.push_back(CPoint(x, y));
		}
	}
}


void Search_min_max(const IsoLines& lines)
{
	IsoLines::const_iteratoriter_lines = lines.begin();
	list<CPoint>::const_iterator p = iter_lines->first.begin();
	x_min = x_max = p->x;
	y_min = y_max = p->y;
	z_min = z_max = iter_lines->second;
	double h;

	for (iter_lines = lines.begin(); iter_lines != lines.end(); ++iter_lines)
	{
		// Циклпоточкамизолиней
		h = iter_lines->second;
		if (h < z_min) z_min = h;
		if (h > z_max) z_max = h;

		for (p = iter_lines->first.begin(); p != iter_lines->first.end(); ++p)
		{
			if (p->x < x_min) x_min = p->x;
			if (p->x > x_max) x_max = p->x;
			if (p->y < y_min) y_min = p->y;
			if (p->y > y_max) y_max = p->y;
		}
	}
}


void TecPlotV(char *filename, const IsoLines& lines)
{
	ofstream out(filename);
	out << "VARIABLES = \"X\" \"Y\" \"Z\"" << endl;
	// Циклпоизолиниям
	IsoLines::const_iteratoriter_lines;
	int line = 1;
	for (iter_lines = lines.begin(); iter_lines != lines.end(); ++iter_lines)
	{
		out << "ZONE T=\"Polygon_" << line++ << "\"" << endl;
		out << "I=" << iter_lines->first.size() << ", J=" << 1 << ", K=1 ZONETYPE=Ordered" << endl;

		out << "DATAPACKING=POINT" << endl;
		out << "DT=(DOUBLE DOUBLE )" << endl;
		out.setf(ios::scientific, ios::floatfield);

		// Циклпоточкамизолиней
		list<CPoint>::const_iterator p;

		for (p = iter_lines->first.begin(); p != iter_lines->first.end(); ++p)
		{
			out << *p << endl;
		}
		out.setf(ios::fixed, ios::floatfield);
	}
	out.close();
}

void TecPlotU(char *filename, const Matrix& U)
{
	ofstream out(filename);
	// Построениесетки
	double hx = 1.0*(x_max - x_min) / N, hy = 1.0*(y_max - y_min) / M, x, y;

	out << "VARIABLES = \"X\" \"Y\" \"Z\"" << endl;
	out << "ZONE T=\"Surface\"" << endl;
	out << "I=" << M + 1 << ", J=" << N + 1 << ", K=1 ZONETYPE=Ordered" << endl;
	out << "DATAPACKING=POINT" << endl;

	out << "DT=(DOUBLE DOUBLE )" << endl;
	out.setf(ios::scientific, ios::floatfield);
	for (int i = 0; i <= N; i++)
	{
		x = x_min + i * hx;
		for (int j = 0; j <= M; j++)
		{
			y = y_min + j * hy;
			out << x << " " << y << " " << U[i][j] << endl;
		}
	}
	out.close();
}


void ApproximateIsoLineInGrid(IsoLines& lines, Matrix& U)
{
	int i, j;
	// Цикл по изолиниям
	IsoLines::iteratoriter_lines;
	ListPoints points;
	ListPoints::iterator it;

	for (iter_lines = lines.begin(); iter_lines != lines.end(); ++iter_lines)
	{
		// Циклпоточкамизолиней
		list<CPoint>::iterator p = iter_lines->first.begin()++;
		list<CPoint>::iterator q = ++p;
		for (p = iter_lines->first.begin(); p != --iter_lines->first.end())
		{
			points.push_back(*p);
			points.push_back(*q);
			SegmentDecomposition(*p, *q, points);
			// Отобразить множество точек на сетку

			for ( (it = points.begin()); it != points.end(); it++)
			{
				SearchPointInGrid(*it, i, j);
				U[i][j] = iter_lines->second;
			}
			points.clear();
			// Передвижение итераторов
			p = q;
			if (q != iter_lines->first.end()) ++q;

		}
	}
}


CPoint SearchPointInGrid(CPoint& P, int& ii, int& jj)
{
	// i, j -номер ячейки, куда попадает точка P
	double hx = 1.0*(x_max - x_min) / N, hy = 1.0*(y_max - y_min) / M;
	double r1, r2, r3, r4, rr1, rr2;
	int i = (P.x - x_min) / hx;
	int j = (P.y - y_min) / hy;
	if(i >= N) {
		i = N - 1;
	}
	if (i < 0) {
		i = 0;
	}
	if (j >= M) {
		j = M - 1;
	}
	if (j < 0) {
		j = 0;
	}
	r1 = P.distance(CPoint(x_min + i * hx, y_min + j * hy));
	r2 = P.distance(CPoint(x_min + (i + 1)*hx, y_min + j * hy));
	r3 = P.distance(CPoint(x_min + i * hx, y_min + (j + 1)*hy));
	r4 = P.distance(CPoint(x_min + (i + 1)*hx, y_min + (j + 1)*hy));

	if (r1 < r2) {
		rr1 = r1;
		ii = i;
	}else {
		rr1 = r2;
		ii = i + 1;
	}

	if (r3 < r4) {
		rr2 = r3;
		ii = i;
	}else {
		rr2 = r4;
		ii = i + 1;
	}

	if(rr1 < rr2) {
		jj = j;
	}else{
		jj = j + 1;
	}

	return CPoint(ii*hx, jj*hy);
}

void SegmentDecomposition(const CPoint&V1, constCPoint&V2, ListPoints&points)
{
	double hx = 1.0*(x_max - x_min) / N, hy = 1.0*(y_max - y_min) / M;
	if (V1.dx(V2) <= hx && V1.dy(V2) <= hy) return;
	CPoint Center = (V1 + V2) / 2;
	points.push_back(Center);
	SegmentDecomposition(V1, Center, points);
	SegmentDecomposition(V2, Center, points);
}

void LineSplineY(Matrix& U)
{
	double hx = 1.0*(x_max - x_min) / N, hy = 1.0*(y_max - y_min) / M, x, y;
	for (int i = 0; i < N; i++)
	{
		int jp = 0; double zp = U[i][jp];
		for (int j = jp + 1; j < M; j++)
		{
			if (fabs(U[i][j] - z_min) < eps) continue;
			//if (fabs(map[i].z[j]-zp)<eps) continue;
			double hz = (U[i][j] - zp) / (j - jp);
			for (int jj = jp; jj <= j; jj++) {
				U[i][jj] = zp + (jj - jp)*hz;
			}

			zp = U[i][j]; jp = j
		}
	}
}


void LineSplineX(Matrix& U)
{
	for (int j = 0; j < M; j++)
	{
		int ip = 0;
		double zp = U[0][j];
		for (int i = ip + 1; i < N; i++)
		{
			if (fabs(U[i][j] - zp) < eps) {
				continue;
			}
			double hz = (U[i][j] - zp) / (i - ip);
			if (i - ip > 1 && hz > 0) {
				for (int ii = ip; ii <= i; ii++) {
					U[ii][j] = zp + (ii - ip)*hz;
				}
			}
			zp = U[i][j]; ip = i;
		}
		ip = N - 1;
		for (int i = ip; i >= 0; i--)
		{
			if (fabs(U[i][j] - zp) < eps) {
				continue;
			}
			double hz = (U[i][j] - zp) / (ip - i);
			if ((ip - i) > 1 && hz > 0) {
				for (int ii = ip; ii >= i; ii--) {
					U[ii][j] = zp + (ip - ii)*hz;
				}
			}
			zp = U[i][j];
			ip=i;
		}
	}
}
