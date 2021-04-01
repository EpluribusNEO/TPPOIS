
#include <iostream>
#include <list>
#include "point.h"
#include <math.h>
#include <windows.h>
#include <vector>
#include <fstream>
#include <string>
#include <GL/gl.h>
#include <GL/glu.h>
#include <glut.h>
//#include "glaux.h"
#include "Header.h"
using namespace std;
#pragma comment (lib,"opengl32.lib")
#pragma comment (lib,"glu32.lib")
#pragma comment (lib, "glaux.lib")
#pragma comment(lib, "legacy_stdio_definitions.lib")
//Описание типов данных
int N = 150, M = 150; // кол-во разбиений по X & Y
int CountT = 2 * N * M;
double x_min(0), x_max(0), y_min(0), y_max(0), z_min(0),
z_max(0);
double hx, hy;
const int CountColors = 20;
double RGB[CountColors][3] =
{
 {0,0,0},
 {0.15,0.07,0},
 {0.27,0.17,0.06},
 {0.44,0.26,0.17},
 {0.40,0.27,0.19},
 {0.37,0.21,0.08},
 {0.56,0.33,0.17},
 {0.45,0.29,0.17},
 {0.63,0.39,0.2},
 {0.79,0.62,0.25},
 {0.85,0.69,0.3},
 {0.77,0.76,0.38},
 {0.58,0.62,0.27},
 {0.35,0.46,0.15},
 {0.18,0.34,0.04},
 {0.1,0.25,0.01},
 {0.03,0.18,0.04},
 {0.0,0.14,0.06},
 {0,0.13,0.18},
 {0,0.09,0.08}
};
float anglez = 45, anglex = 0, rz = 200;
float ex = 0, ey = 10, ez = 0;
float ax = 0, ay = 0, az = 0;
struct Vershina //випись вершин сетки
{
	float R1, R2, H1, H2, H;
	bool b1, b2;
};
typedef vector<vector<Vershina>> VectorNodes;
typedef vector < vector<bool>> VectorRebro;
typedef list<pair<list<CPoint>, double> > IsoLines;
VectorNodes Nodes;
VectorRebro RebroV, RebroG, RebroN, RebroS;
IsoLines lines;
void Load(string filename, IsoLines& lines);//Ɂɚɝɪɭɡɤɚ
void Search_min_max(const IsoLines& lines);
void statistics(const IsoLines& points);
//стандартная функция для графики (отрисовка поверхностей)
void CALLBACK PressKeyA(void) // шаг по кругу
{
	anglez += 0.1;
}
void CALLBACK PressKeyD(void) // шаг по кругу
{
	anglez -= 0.1;
}
void CALLBACK PressKeyW(void) // ɲшаг вперёд
{
	rz -= 1;
}
void CALLBACK PressKeyS(void)// ɲшаг назад
{
	rz += 1;
}
void CALLBACK PressKeyQ(void) // шаг вверх
{
	ey += 10;
}
void CALLBACK PressKeyE(void) // шаг вниз
{
	ey -= 10;
}
void CALLBACK resize(int width, int height)
{
	glViewport(0, 0, width, height);//область вывода Open-GL граффики, совпадает с размером окна
	glMatrixMode(GL_PROJECTION);//режим работы с матрицей проекции
	glLoadIdentity();//загрузка единичной матрицы --сброс текущей
	gluPerspective(1500, double(width) / height, 1, 1000);//установить тип проекции. Задать перспективную проекцию

	glMatrixMode(GL_MODELVIEW);//включить режим работы с модельно-видовой матрицей

}
void system(int n) //функция - система координат
{
	glLineWidth(2);
	//ось X - красный
	glBegin(GL_LINES);
	glColor3d(1, 0, 0);
	glVertex3d(-n, 0, 0);
	glVertex3d(n, 0, 0);
	glEnd();
	//стрелки на X
	glBegin(GL_LINE_STRIP);
	glVertex3d(n - 2, 2, 0);
	glVertex3d(n, 0, 0);
	glVertex3d(n - 2, -2, 0);
	glEnd();
	//ось Y - синий
	glBegin(GL_LINES);
	glColor3d(0, 0, 1);
	glVertex3d(0, -n, 0);
	glVertex3d(0, n, 0);
	glEnd();
	//стрелки на Y
	glBegin(GL_LINE_STRIP);
	glVertex3d(-2, n - 2, 0);
	glVertex3d(0, n, 0);
	glVertex3d(2, n - 2, 0);
	glEnd();
	//ось я - зелёный
	glBegin(GL_LINES);
	glColor3d(0, 1, 0);
	glVertex3d(0, 0, n);
	glVertex3d(0, 0, -n);
	glEnd();
	//стрелки на Z
	glBegin(GL_LINE_STRIP);
	glVertex3d(0, -2, n - 2);
	glVertex3d(0, 0, n);
	glVertex3d(0, 2, n - 2);
	glEnd();
}
void Load(string filename, IsoLines& lines)//загрузка
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
			if (Loading.eof()) return;
		} while (str != "[POLYLINE]");
		getline(Loading, str, '='); // слдово Type
		getline(Loading, str, '\n'); // значение
		if (!(str == "0x20" || str == "0x21" || str == "0x22"))
		{
			continue;
		}
		getline(Loading, str, '='); // Label
		if (str != "Label") continue;
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
				if (symbol == 'D')
					getline(Loading, str, '=');
			} while (!isdigit(Loading.peek()) && symbol !='[');
			lines.back().first.push_back(CPoint(x * 67000, y * 111000));
			if (symbol == '[' || symbol == 'E') break;
		}
		getline(Loading, str, '\n');
		Loading.get();
	}
}
void Search_min_max(const IsoLines& lines)
{
	IsoLines::const_iterator iter_lines = lines.begin();//IsoLines::const_iterator iter_lines = lines.begin();
	iter_lines++; // <<---------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<------ЭТО-Я-ДОБАВИЛ------------------<<<<<<<<<<<<<<<<<<<<<<<<<
	list<CPoint>::const_iterator p = iter_lines->first.begin(); //first.begin();
	iter_lines = lines.begin(); // <<---------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<------ЭТО-Я-ДОБАВИЛ-----------<<<<<<<<<<<<<<<<<<<<<<<<
	x_min = x_max = p->x;
	y_min = y_max = p->y;
	z_min = z_max = iter_lines->second;
	double h;
	for (iter_lines = lines.begin(); iter_lines != lines.end(); ++iter_lines)
	{
		// ɐɢɤɥ ɩɨ ɬɨɱɤɚɦ ɢɡɨɥɢɧɟɣ
		h = iter_lines->second;
		if (h < z_min)
			z_min = h;
		if (h > z_max)
			z_max = h;
		
		for (p = iter_lines->first.begin(); p != iter_lines->first.end(); ++p)
		{
			if (p->x < x_min) x_min = p->x;
			if (p->x > x_max) x_max = p->x;
			if (p->y < y_min) y_min = p->y;
			if (p->y > y_max) y_max = p->y;
		}
	}
	hx = (x_max - x_min) / N;
	hy = (y_max - y_min) / M;
}
void statistics(const IsoLines& points)
{
	cout << "Ʉɨɥɢɱɟɫɬɜɨ ɢɡɨɥɢɧɢɣ = " << points.size() <<
		endl;
	IsoLines::const_iterator iter_lines;
	long long Count(0);
	for (iter_lines = points.begin(); iter_lines !=
		points.end(); ++iter_lines)
	{
		Count += iter_lines->first.size();
	}
	cout << "Ɉɛɳɟɟ ɤɨɥɢɱɟɫɬɨ ɬɨɱɟɤ = " << Count << endl;
	cout << "Ɇɢɧɢɦɚɥɶɧɚɹ ɜɵɫɨɬɚ = " << z_min << endl;
	cout << "Ɇɚɤɫɢɦɚɥɶɧɚɹ ɜɵɫɨɬɚ = " << z_max << endl;
}
bool Peresechenie(double x1, double y1, double x2, double
	y2, double x3, double y3, double x4, double y4, double &r1,
	double &r2)
{
	double v1, v2, v3, v4;
	v1 = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3);
	v2 = (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3);
	v3 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
	v4 = (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1);
	double Px, Py;
	Px = x1 + (x2 - x1) * fabs(v1) / fabs(v2 - v1);
	Py = y1 + (y2 - y1) * fabs(v1) / fabs(v2 - v1);
	if ((v1*v2 <= 0) && (v3*v4 <= 0))
	{
		r1 = sqrt((Px - x3) * (Px - x3) + (Py - y3) * (Py-y3));
		r2 = sqrt((Px - x4) * (Px - x4) + (Py - y4) * (Py-y4));
		return true;
	}
	else
	{
		return false;
	}
}
bool SetVershina(int i, int j, double distance, double h, VectorRebro& Rebro)
{
	Rebro[i][j] = true;
	if (!Nodes[i][j].b1)// 1 случай
	{
		Nodes[i][j].b1 = true;
		Nodes[i][j].R1 = distance;
		Nodes[i][j].H1 = h;
		return true;
	}
	else if (!Nodes[i][j].b2) // 2 случай
	{
		if (Nodes[i][j].H1 != h)
		{
			Nodes[i][j].b2 = true;
			Nodes[i][j].R2 = distance;
			Nodes[i][j].H2 = h;
			return true;
		}
		else
		{
			Nodes[i][j].b1 = true;
			Nodes[i][j].R1 = min(distance, Nodes[i][j].R1);
			Nodes[i][j].H1 = h;
			return true;
		}
	}
	else
		if (!(distance > Nodes[i][j].R1 && distance >
			Nodes[i][j].R2)) // 3 ɫɥɭɱɚɣ
		{
			if ((Nodes[i][j].H1 != h) && (Nodes[i][j].H2 != h))
				// если все 3 узла разных уровней
			{
				if (Nodes[i][j].R1 < Nodes[i][j].R2) // 
				{
					Nodes[i][j].R2 = distance;
					Nodes[i][j].H2 = h;
					return true;
				}
				else
				{
					Nodes[i][j].R1 = distance;
					Nodes[i][j].H1 = h;
					return true;
				}
			}
			else // если уровень нового узла равен какого-то другого из 2
			{
				if (Nodes[i][j].H1 == h) // если уровень нового узла равен уровню 1 узла
				{
				Nodes[i][j].R1 = min(distance,
			   Nodes[i][j].R1);
				Nodes[i][j].H1 = h;
				return true;
				}
					if (Nodes[i][j].H2 == h)// если уровень нового узла равен уровню 2
					{
					Nodes[i][j].R2 = min(distance,
				   Nodes[i][j].R2);
					Nodes[i][j].H2 = h;
					return true;
					}
			}
		}
	return false;
}
void FirstStage(IsoLines& lines, VectorRebro& Rebro, int
	ykaz)
{
	bool f;
	int i, j, i0, j0, i1, j1;
	double distance1, distance2;
	double h, x_1, x_2, y_1, y_2, x3, x4, y3, y4;
	IsoLines::iterator iter_lines;
	for (iter_lines = lines.begin(); iter_lines != lines.end(); ++iter_lines)
	{
		h = iter_lines->second;
		// цикл по точкам изолиний
		list<CPoint>::iterator p = iter_lines->first.begin()++;
		list<CPoint>::reverse_iterator p_end;
		list<CPoint>::iterator q = ++p;
		for (p = iter_lines->first.begin(); p != --iter_lines->first.end(); )
		{
			x_1 = min(p->x, q->x);
			x_2 = max(p->x, q->x);
			y_1 = min(p->y, q->y);
			y_2 = max(p->y, q->y);
			i0 = (x_1 - x_min) / hx;
			i1 = (x_2 - x_min) / hx;
			j0 = (y_1 - y_min) / hy;
			j1 = (y_2 - y_min) / hy;
			if (ykaz == 1)
			{
				for (i = i0; i <= i1 && i <= N; i++)
					for (j = j0; j <= j1 && j < M; j++)
					{
						x3 = x4 = x_min + i * hx;
						y3 = y_min + j * hy;
						y4 = y_min + (j + 1) * hy;
						f = Peresechenie(p->x, p->y, q->x, q->y,
							x3, y3, x4, y4, distance1, distance2);
						if (f)
						{
							Rebro[i][j] = true;
							SetVershina(i, j, distance1, h, Rebro);
							SetVershina(i, j + 1, distance2, h,
								Rebro);
						}
					}
			}
			else
				if (ykaz == 2)
				{
					for (i = i0; i <= i1 && i < N; i++)
						for (j = j0; j <= j1 && j <= M; j++)
						{
							x3 = x_min + i * hx;
							x4 = x_min + (i + 1) * hx;
							y3 = y4 = y_min + j * hy;
							f = Peresechenie(p->x, p->y, q->x, q->y, x3, y3, x4, y4, distance1, distance2);
							if (f)
							{
								Rebro[i][j] = true;
								SetVershina(i, j, distance1, h,
									Rebro);
								SetVershina(i + 1, j, distance2, h,
									Rebro);
							}
						}
				}
				else
					if (ykaz == 3)
					{
						for (i = i0; i <= i1 && i < N; i++)
							for (j = j0; j <= j1 && j < M; j++)
							{
								x3 = x_min + i * hx;
								x4 = x_min + (i + 1) * hx;
								y3 = y_min + j * hy;
								y4 = y_min + (j + 1) * hy;
								f = Peresechenie(p->x, p->y, q->x, q->y, x3, y3, x4, y4, distance1, distance2);
								if (f)
								{
									Rebro[i][j] = true;
									SetVershina(i, j, distance1, h,
										Rebro);
									SetVershina(i + 1, j + 1,
										distance2, h, Rebro);
								}
							}
					}
					else
						if (ykaz == 4)
						{
							for (i = i0; i <= i1 && i < N; i++)
								for (j = j0; j <= j1 && j < M; j++)
								{
									x3 = x_min + i * hx;
									x4 = x_min + (i + 1) * hx;
									y3 = y_min + (j + 1) * hy;
									y4 = y_min + j * hy;
									f = Peresechenie(p->x, p->y, q->x, q->y, x3, y3, x4, y4, distance1, distance2);
									if (f)
									{
										Rebro[i][j] = true;
											SetVershina(i, j + 1, distance1, h, Rebro);
										SetVershina(i + 1, j,distance2, h, Rebro);
									}
								}
						}
			p = q; // передвижение итератора
			if (q != iter_lines->first.end())
				++q;
		}
		p = iter_lines->first.begin();
		p_end = iter_lines->first.rbegin();
		x_1 = min(p->x, p_end->x);
		x_2 = max(p->x, p_end->x);
		y_1 = min(p->y, p_end->y);
		y_2 = max(p->y, p_end->y);
		i0 = (x_1 - x_min) / hx;
		i1 = (x_2 - x_min) / hx;
		j0 = (y_1 - y_min) / hy;
		j1 = (y_2 - y_min) / hy;
		if (ykaz == 1)
		{
			for (i = i0; i <= i1 && i <= N; i++)
				for (j = j0; j <= j1 && j < M; j++)
				{
					x3 = x4 = x_min + i * hx;
					y3 = y_min + j * hy;
					y4 = y_min + (j + 1) * hy;
					f = Peresechenie(p->x, p->y, p_end->x, p_end->y, x3, y3, x4, y4, distance1, distance2);
					if (f)
					{
						Rebro[i][j] = true;
						SetVershina(i, j, distance1, h, Rebro);
						SetVershina(i, j + 1, distance2, h, Rebro);
					}
				}
		}
		else
			if (ykaz == 2)
			{
				for (i = i0; i <= i1 && i < N; i++)
					for (j = j0; j <= j1 && j <= M; j++)
					{
						x3 = x_min + i * hx;
						x4 = x_min + (i + 1) * hx;
						y3 = y4 = y_min + j * hy;
						f = Peresechenie(p->x, p->y, p_end->x,
							p_end->y, x3, y3, x4, y4, distance1, distance2);
						if (f)
						{
							Rebro[i][j] = true;
							SetVershina(i, j, distance1, h,
								Rebro);
							SetVershina(i + 1, j, distance2, h,
								Rebro);
						}
					}
			}
			else
				if (ykaz == 3)
				{
					for (i = i0; i <= i1 && i < N; i++)
						for (j = j0; j <= j1 && j < M; j++)
						{
							x3 = x_min + i * hx;
							x4 = x_min + (i + 1) * hx;
							y3 = y_min + j * hy;
							y4 = y_min + (j + 1) * hy;
							f = Peresechenie(p->x, p->y, p_end->x, p_end->y, x3, y3, x4, y4, distance1, distance2);
							if (f)
							{
								Rebro[i][j] = true;
								SetVershina(i, j, distance1, h,
									Rebro);
								SetVershina(i + 1, j + 1, distance2,
									h, Rebro);
							}
						}
				}
				else
					if (ykaz == 4)
					{
						for (i = i0; i <= i1 && i < N; i++)
							for (j = j0; j <= j1 && j < M; j++)
							{
								x3 = x_min + i * hx;
								x4 = x_min + (i + 1) * hx;
								y3 = y_min + (j + 1) * hy;
								y4 = y_min + j * hy;
								f = Peresechenie(p->x, p->y,
									p_end->x, p_end->y, x3, y3, x4, y4, distance1, distance2);
								if (f)
								{
									Rebro[i][j] = true;
									SetVershina(i, j + 1, distance1,
										h, Rebro);
									SetVershina(i + 1, j, distance2,
										h, Rebro);
								}
							}
					}
	}
}
void SecondStage(IsoLines& lines, VectorRebro& RebroV,
	VectorRebro& RebroG, VectorRebro& RebroN, VectorRebro&
	RebroS, bool& otmetka)
{
	int i, j, jj, ii;
	double NewDistance0, NewDistance, NewH;
	otmetka = false; // предпологаем, что больше не осталось непомеченных рёбер т.е все имеют репесечения 
		// проход по всем узламɦ
	for (i = 0; i <= N; i++)
		for (j = 0; j <= M; j++)
		{
			if (Nodes[i][j].b1 == true || Nodes[i][j].b2 == true) // если узлу уже приписаны расстояния, хотя бы 1 
			{
				if (Nodes[i][j].b1 && Nodes[i][j].b2)
				{
					if (Nodes[i][j].R1 < Nodes[i][j].R2)
					{
						NewDistance = Nodes[i][j].R1;
						NewH = Nodes[i][j].H1;
					}
					else
					{
						NewDistance = Nodes[i][j].R2;
						NewH = Nodes[i][j].H2;
					}
				}
				else if (Nodes[i][j].b1)
				{
					NewDistance = Nodes[i][j].R1;
					NewH = Nodes[i][j].H1;
				}
				else
				{
					NewDistance = Nodes[i][j].R2;
					NewH = Nodes[i][j].H2;
				}
				//выпускае волну в трёх итерациях
				if (i > 0 && j > 0 && !RebroN[i - 1][j - 1])
					// Ю-С --влево внизу
				{
					ii = i - 1;
					jj = j - 1;
					NewDistance0 = NewDistance;
					while (RebroN[ii][jj] == false && ii > 0 &&	jj > 0)// пока ребро небыло отмечено
					{
						NewDistance0 += sqrt(hy*hy + hx * hx);
						otmetka += SetVershina(ii, jj,
							NewDistance0, NewH, RebroN);
						ii--;
						jj--;
					}
				}
				if (i > 0 && !RebroG[i - 1][j]) //ɜɥɟɜɨ
				{
					ii = i - 1;
					NewDistance0 = NewDistance;
					while (ii > 0 && RebroG[ii][j] == false)
						// пока ребро небыло отмечено
					{
						NewDistance0 += hx;
						otmetka += SetVershina(ii, j,
							NewDistance0, NewH, RebroG);
						ii--;
					}
				}
				if (i > 0 && j < M && !RebroS[i - 1][j + 1])
					// C-Ю --влево вверху
				{
					ii = i - 1;
					jj = j + 1;
					NewDistance0 = NewDistance;
					while (ii > 0 && jj < M && RebroS[ii][jj] == false) //пока ребро небыло отмечено
					{
						NewDistance0 += sqrt(hy * hy + hx * hx);
						otmetka += SetVershina(ii, jj,
							NewDistance0, NewH, RebroS);
						ii--;
						jj++;
					}
				}
				if (j < M && !RebroV[i][j + 1])// вверх
				{
					jj = j + 1;
					NewDistance0 = NewDistance;
					while (jj < M && RebroV[i][jj] == false)
						// пока ребро небыло отмечено
					{
						NewDistance0 += hy;
						otmetka += SetVershina(i, jj,
							NewDistance0, NewH, RebroV);
						jj++;
					}
				}
				if (i < N && j < M && !RebroN[i + 1][j + 1])// вверх
				{
					ii = i + 1;
					jj = j + 1;
					NewDistance0 = NewDistance;
					while (ii < N && jj < M && RebroN[ii][jj]
						== false) // пока ребро небыло отмечено
					{
						NewDistance0 += sqrt(hy * hy + hx * hx);
						otmetka += SetVershina(ii, jj,
							NewDistance0, NewH, RebroN);
						ii++;
						jj++;
					}
				}
				if (i < N && !RebroG[i + 1][j])// вправо
				{
					ii = i + 1;
					NewDistance0 = NewDistance;
					while (ii < N && RebroG[ii][j] == false) //пока ребро небыло отмечено
					{
					NewDistance0 += hx;
					otmetka += SetVershina(ii, j,
				   NewDistance0, NewH, RebroG);
					ii++;
					}
				}
				if (i < N && j>0 && !RebroS[i + 1][j - 1])//
				{
					ii = i + 1;
					jj = j - 1;
					NewDistance0 = NewDistance;
					while (ii < N && jj > 0 && RebroS[ii][jj] == false) // пока ребро небыло отмечено
					{
						NewDistance0 += sqrt(hy * hy + hx * hx);
						otmetka += SetVershina(ii, jj,
							NewDistance0, NewH, RebroS);
						ii++;
						jj--;
					}
				}
				if (j > 0 && !RebroV[i][j - 1])// вниз
				{
					jj = j - 1;
					NewDistance0 = NewDistance;
					while (jj > 0 && RebroV[i][jj] == false)
						// пока ребро небыло отмечено
					{
						NewDistance0 += hy;
						otmetka = SetVershina(i, jj,
							NewDistance0, NewH, RebroV);
						jj--;
					}
				}
			}
		}
}
void DrawTriangle(GLdouble x1, GLdouble y1, GLdouble z1,
	GLdouble x2, GLdouble y2, GLdouble z2, GLdouble x3,
	GLdouble y3, GLdouble z3, GLdouble r1, GLdouble g1,
	GLdouble b1, GLdouble r2, GLdouble g2, GLdouble b2,
	GLdouble r3, GLdouble g3, GLdouble b3)
{ // передаются координаты и цвета трёх точнк треугольника
	glBegin(GL_TRIANGLES);
	glColor3d(r1, g1, b1); glVertex3d(x1, y1, z1);
	glColor3d(r2, g2, b2); glVertex3d(x2, y2, z2);
	glColor3d(r3, g3, b3); glVertex3d(x3, y3, z3);
	glEnd();
}
void TRIANGLES_Kontur(float x1, float y1, float z1, float
	x2, float y2, float z2, float x3, float y3, float z3, float
	r, float b, float g)
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//режим
	glBegin(GL_TRIANGLES);
	glColor3d(r, g, b);
	glVertex3f(x1, y1, z1);
	glVertex3f(x2, y2, z2);
	glVertex3f(x3, y3, z3);
	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//восстановление режима поумолчанию - заливка
}
void drawtriangle()
{
	//цикл по треугольникам
	float H;
	int k1, k2, k3;
	float r1, r2, r3, g1, g2, g3, b1, b2, b3;
	H = (z_max - z_min) / CountColors;
	//масштабирование
	float DX = (x_max + x_min) / 2;
	float DY = (y_max + y_min) / 2;
	float DZ = (z_max + z_min) / 2;
	float MAXXX = max((z_max - z_min) / 2, max((x_max - x_min) / 2, (y_max - y_min) / 2));
	float K = hy / hx;
	for (int i = 1; i < N - 1; i++)
		for (int j = 1; j < M - 1; j++)
		{
			// нижний треугольник
			k1 = (z_max - Nodes[i][j].H) / H;
			if (k1 == CountColors) k1 = CountColors - 1;
			r1 = RGB[k1][0]; g1 = RGB[k1][1]; b1 = RGB[k1][2];
			k2 = (z_max - Nodes[i + 1][j].H) / H;
			if (k2 == CountColors) k2 = CountColors - 1;
			r2 = RGB[k2][0]; g2 = RGB[k2][1]; b2 = RGB[k2][2];
			k3 = (z_max - Nodes[i + 1][j + 1].H) / H;
			if (k3 == CountColors) k3 = CountColors - 1;
			r3 = RGB[k3][0]; g3 = RGB[k3][1]; b3 = RGB[k3][2];
			//нарисовать треугольник
			DrawTriangle(
				(x_min + i * hx - DX) / MAXXX * 100 * K, (y_min + j * hy-DY) / MAXXX * 100, (Nodes[i][j].H - DZ) / MAXXX * 100 * 2,
				(x_min + (i + 1) * hx - DX) / MAXXX * 100 * K,
				(y_min + j * hy - DY) / MAXXX * 100, (Nodes[i + 1][j].H - DZ)
				/ MAXXX * 100 * 2,
				(x_min + (i + 1) * hx - DX) / MAXXX * 100 * K,
				(y_min + (j + 1) * hy - DY) / MAXXX * 100, (Nodes[i + 1][j + 1].H
					- DZ) / MAXXX * 100 * 2,
				r1, g1, b1,
				r2, g2, b2,
				r3, g3, b3);
			// верхний треугольник
			k1 = (z_max - Nodes[i][j].H) / H;
			if (k1 == CountColors) k1 = CountColors - 1;
			r1 = RGB[k1][0]; g1 = RGB[k1][1]; b1 = RGB[k1][2];
			k2 = (z_max - Nodes[i][j + 1].H) / H;
			if (k2 == CountColors) k2 = CountColors - 1;
			r2 = RGB[k2][0]; g2 = RGB[k2][1]; b2 = RGB[k2][2];
			k3 = (z_max - Nodes[i + 1][j + 1].H) / H;
			if (k3 == CountColors) k3 = CountColors - 1;
			r3 = RGB[k3][0]; g3 = RGB[k3][1]; b3 = RGB[k3][2];
			//нарисовать треугольник
			DrawTriangle(
				(x_min + i * hx - DX) / MAXXX * 100 * K, (y_min
					+ j * hy - DY) / MAXXX * 100, (Nodes[i][j].H - DZ) / MAXXX
				* 100 * 2,
				(x_min + i * hx - DX) / MAXXX * 100 * K, (y_min
					+ (j + 1) * hy - DY) / MAXXX * 100, (Nodes[i][j + 1].H - DZ) /
				MAXXX * 100 * 2,
				(x_min + (i + 1) * hx - DX) / MAXXX * 100 * K,
				(y_min + (j + 1) * hy - DY) / MAXXX * 100, (Nodes[i + 1][j
					+ 1].H - DZ) / MAXXX * 100 * 2,
				r1, g1, b1,
				r2, g2, b2,
				r3, g3, b3);
		}
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
void scene()
{
	glRotated(-90, 1, 0, 0);
	glRotatef(0.45f, 0.0f, 0.0f, 1.0f); //поворачивает по шкале z на 45
	glPushMatrix(); //сохранть текущую систему соординат
	drawtriangle();
	glPopMatrix();
}


void CALLBACK display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// Очистить буфер-вета и буфер-глубины
	glLoadIdentity();
	ex = rz * sin(anglez);
	ez = rz * cos(anglez);
	gluLookAt(ex, ey, ez, ax, ay, az, 0, 1, 0);
	scene();
	auxSwapBuffers();//вывести срдержимое буфера на экран
}


void RunOpenGL()
{
	auxKeyFunc(AUX_w, PressKeyW);
	auxKeyFunc(AUX_s, PressKeyS);
	auxKeyFunc(AUX_a, PressKeyA);
	auxKeyFunc(AUX_d, PressKeyD);
	auxKeyFunc(AUX_q, PressKeyQ);
	auxKeyFunc(AUX_e, PressKeyE);
	auxInitPosition(0, 0, 1000, 1000);//--------------------------------
	auxInitDisplayMode(AUX_RGB | AUX_DEPTH | AUX_DOUBLE);//
	auxInitWindow(L"Glaux Template");// ---------------------------------
	// ----------------------------------------------------
		//-----------------------------------------
		auxIdleFunc(display);
	//------------------------------------------
	auxReshapeFunc(resize);
	glClearColor(1, 1, 1, 0); //-----------------------------
	glEnable(GL_DEPTH_TEST);//--------------------------
}


//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
void main()
{
	setlocale(LC_ALL, "RUS");
	int ykaz; Search_min_max(lines); statistics(lines); 
	int i, j;
	Load("borus_full_2.txt", lines); //--------------------------------

	cout << "Изолинии загружены!\n";
	cout << "___________________ЭТАП №1___________________" << endl;
	Nodes.resize(N + 1);
	RebroV.resize(N + 1);
	RebroG.resize(N + 1);
	RebroN.resize(N + 1);
	RebroS.resize(N + 1);
	for (i = 0; i <= N; i++)
	{
		Nodes[i].resize(M + 1);
		RebroV[i].resize(M + 1);
		RebroG[i].resize(M + 1);
		RebroN[i].resize(M + 1);
		RebroS[i].resize(M + 1);
		for (j = 0; j <= M; j++)
		{
			Nodes[i][j].b1 = Nodes[i][j].b2 = false;
			Nodes[i][j].R1 = Nodes[i][j].R2 = 0;
			RebroV[i][j] = false;
			RebroG[i][j] = false;
			RebroN[i][j] = false;
			RebroS[i][j] = false;
		}
	}
	cout << "_________Пересечения с вертикальными рёбрами сетки_____________ " << endl;
	FirstStage(lines, RebroV, 1);
	
	cout << "________Горизонтальными с вертикальными рёбрами сетки______________ " << endl;
	FirstStage(lines, RebroG, 2);
	
	cout << "______Пересечения с диагональными Ю-С рёбрами сетки________________ " << endl;
	FirstStage(lines, RebroN, 3);
	
	cout << "_______Пересечения с диагональными С-Ю рёбрами сетки_______________ " << endl;
	FirstStage(lines, RebroS, 4);
	
	cout << "________________Этап_#2___________________" << endl;
	
	bool otmetka = true;
	while (otmetka)
	{
		SecondStage(lines, RebroV, RebroG, RebroN,RebroS, otmetka);
	}
	//вычисляем ɜвысоту в узлах с двумя значениями
	for (i = 0; i <= N; i++)
		for (j = 0; j <= M; j++)
		{
			if (Nodes[i][j].b1 && Nodes[i][j].b2)
			{
				Nodes[i][j].H = (Nodes[i][j].H1*
					Nodes[i][j].R2 + Nodes[i][j].H2 *
					Nodes[i][j].R1) / (Nodes[i][j].R1 + Nodes[i][j].R2);
			}
			else if (Nodes[i][j].b1)
			{
				Nodes[i][j].H = Nodes[i][j].H1;
			}
			else
			{
				Nodes[i][j].H = Nodes[i][j].H2;
			}
		}
	RunOpenGL();
	auxMainLoop(display); //главцый цикл обработки событий
	system("pause");
}

/* end of line */
