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
    return sin(t) / (0.1 + x * x);
}


// пространство имен std - включаючает в себя большинство стандартных функций
using namespace std;
void main( void )
{
    double a = 0;             // нижний предел интегрирования
    double b = 1;             // верхний предел интегрирования
    double abserr = 1.0e-14;  // абсолютная погрешность
    double relerr = 0;        // относительная погрешность
    
    double errest;            // оценка погрешности
    int nofun;                // количество вычислений подынтегральной функции
    double flag;              // индикатор надежности

    // считаем в точках 0.1 + 0.2 * k при помощи QUANC8
    x = 0.1;
    double results1[10];
    for (int i = 0; i <= 9; i++)
    {
        quanc8(function, a, b, abserr, relerr, &results1[i], &errest, &nofun, &flag);
        x += 0.2;
    }
    x = 0;

    // вычислим X[i], Y[i] для полинома Лагранжа
    double lagrangeX[10];
    double lagrangeY[10];

    // вычислим X[i], Y[i] для сплайна
    double splineX[11];
    double splineY[11];

    double j = 0;

    for (int i = 0; i <= 9; i++)
    {
        lagrangeX[i] = j;
        splineX[i + 1] = j;
        j += 0.2;
    }

    for (int i = 0; i <= 9; i++)
    {
        quanc8(function, a, b, abserr, relerr, &lagrangeY[i], &errest, &nofun, &flag);
        quanc8(function, a, b, abserr, relerr, &splineY[i + 1], &errest, &nofun, &flag);
        x += 0.2;
    }

    // точка, в которой необходимо вычислить значения полином Лагранжа и сплайна, берём между узлами
    double tr = 0.1;
    double results2[10];

    // считаем в точках 0.1 + 0.2 * k при помощи полинома Лагранжа
    for (int i = 0; i <= 9; i++)
    {
        results2[i] = lagrange(10, lagrangeX, lagrangeY, tr);
        tr += 0.2;
    }

    tr = 0.1;
    double results3[10];

    // вычислим коэффициенты для сплайна
    double splineB[11], splineC[11], splineD[11];
    spline(10, splineX, splineY, splineB, splineC, splineD);

    // считаем в точках 0.1 + 0.2 * k при помощи сплайна
    for (int i = 0; i <= 9; i++)
    {
        results3[i] = seval(10, &tr, splineX, splineY, splineB, splineC, splineD);
        tr += 0.2;
    }

    // вывод
    printf("   x       ");
    printf("QUANC8           ");
    printf("Lagrange           ");
    printf("Spline        ");
    printf("Lagrange Error     ");
    printf("Spline Error \n");
    for (int i = 0; i <= 9; i++) {
        printf("%5.1f ", (0.1 + 0.2 * i));
        printf("%14.10f ", results1[i]);
        printf("%17.10f ", results2[i]);
        printf("%17.10f ", results3[i]);
        printf("%17.10f ", abs(results2[i] - results1[i]));
        printf("%17.10f \n", abs(results3[i] - results1[i]));
    }
}