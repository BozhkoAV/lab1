#include <stdio.h>
#include "lagrange.h"
#include "quanc8.h"
#include "splines.h"
#include <math.h>
#include <iostream>

double x = 0;
double t;

// Подынтергальная функция
double function(double t)
{
    return 1 + pow(2 * t * x, 2);
}


// пространство имен std - включаючает в себя большинство стандартных функций
using namespace std;
void main( void )
{
    double a = 0;             // нижний предел интегрирования
    double b = 1;             // верхний предел интегрирования
    double abserr = 1.0e-14;  // абсолютная погрешность
    double relerr = 0;  // относительная погрешность
    
    double errest;            // оценка погрешности
    int nofun;                // количество вычислений подынтегральной функции
    double flag;              // индикатор надежности

    // считаем в точках 0.1 + 0.2 * k при помощи QUANC8
    x = 0.1;
    double results1[11];
    for (int i = 0; i <= 10; i++)
    {
        quanc8(function, a, b, abserr, relerr, &results1[i], &errest, &nofun, &flag);
        x += 0.2;
    }
    x = 0;

    // вычислим X[i], Y[i] для полинома Лагранжа
    double lagrangeX[11];
    double lagrangeY[11];

    // вычислим X[i], Y[i] для сплайна
    double splineX[12];
    double splineY[12];

    double j = 0;

    for (int i = 0; i <= 10; i++)
    {
        lagrangeX[i] = j;
        splineX[i + 1] = j;
        j += 0.2;
    }

    for (int i = 0; i <= 10; i++)
    {
        quanc8(function, a, b, abserr, relerr, &lagrangeY[i], &errest, &nofun, &flag);
        quanc8(function, a, b, abserr, relerr, &splineY[i + 1], &errest, &nofun, &flag);
        x += 0.2;
    }

    // точка, в которой необходимо вычислить значения полином Лагранжа и сплайна, берём между узлами
    double tr = 0.1;
    double results2[11];

    // считаем в точках 0.1 + 0.2 * k при помощи полинома Лагранжа
    for (int i = 0; i <= 10; i++)
    {
        results2[i] = lagrange(11, lagrangeX, lagrangeY, tr);
        tr += 0.2;
    }

    tr = 0.1;
    double results3[11];

    // вычислим коэффициенты для сплайна
    double splineB[12], splineC[12], splineD[12];
    spline(11, splineX, splineY, splineB, splineC, splineD);

    // считаем в точках 0.1 + 0.2 * k при помощи полинома Лагранжа
    for (int i = 0; i <= 10; i++)
    {
        results3[i] = seval(11, &tr, splineX, splineY, splineB, splineC, splineD);
        tr += 0.2;
    }

    // вывод
    cout << "   x           QUANC8         Lagrange        Spline" << endl;
    for (int i = 0; i <= 10; i++) {
        cout.width(5);
        cout << 0.1 + 0.2 * i << " ";
        cout.width(15);
        cout << results1[i] << " ";
        cout.width(15);
        cout << results2[i] << " ";
        cout.width(15);
        cout << results3[i] << endl;
    }
}