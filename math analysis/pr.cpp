#include <iostream>
#include <string>
#include <cmath>
#include <sstream>

// Lagrange polynom
using namespace std;
const int sym = 3; // количество знаков после запятой

template <typename Type> // преобразование в строку
string to_str(const Type &t)
{
    ostringstream os;
    os << t;
    return os.str();
}

void after_point(double &num)
{
    int ar = pow(10, sym);
    int modulo = (int)(num * ar) % ar;    // остаток от деления - знаки после запятой
    num = (int)num + (double)modulo / ar; // число с заданной точностью
}

// вычисление полинома
double *calc_func(double *lag, double *save_lag, double *clone_lag, double coeff_sum_el, double *x, double *y, int size, int lag_pos)
{
    for (int q = 0; q < size; q++)
        lag[q] = 0;

    for (int i = 0; i < size; i++)
    {
        for (int q = 0; q < size; q++)
            save_lag[q] = lag[q]; // копируем

        for (int q = 0; q < size - 1; q++)
            lag[q] = 0;    // обнуляем всё, кроме крайнего элемента
        lag[size - 1] = 1; // это изначальный многочлен f(x)=1

        coeff_sum_el = 1;
        lag_pos = size - 1;
        for (int j = 0; j < size; j++)
        {
            // на каждой итераци умножается каждая новая компонента
            // сначала двигается полином, который был получен в результате умножения
            // на предыдущей итерации
            // итерации: 1 двигаем * (x-3)
            //           (x-3) двигаем * (x-4)
            //           (x^2-7x+12) двигаем * (x-5)
            if (i != j)
            {
                lag_pos--;
                coeff_sum_el *= (x[i] - x[j]);
                for (int q = lag_pos; q < size - 1; q++)
                    lag[q] = lag[q + 1]; // двигаем

                // обнуляем крайний после сдвига
                lag[size - 1] = 0;

                for (int q = 0; q < size; q++)
                    clone_lag[q] = lag[q];
                // получаем новый полином
                for (int q = 1; q < size; q++)
                    lag[q] += -x[j] * clone_lag[q - 1];
            }
        }
        coeff_sum_el = y[i] / coeff_sum_el;
        // полный новый полином
        for (int q = 0; q < size; q++)
        {
            lag[q] *= coeff_sum_el;
            lag[q] += save_lag[q];
        }
        // сохраняем для дальнейших действий
        for (int q = 0; q < size; q++)
            save_lag[q] = lag[q];
    }
    return lag;
}

// проверка полинома на нулевой
bool check_func(double *lag, int size)
{
    for (int i = 0; i < size; i++)
        if (lag[i] != 0)
            return false;
    return true;
}

// вычисление значения полинома
double value_func(double *lag, double arg, int size)
{
    double value = 0;
    for (int i = 0; i < size; i++)
        value += lag[i] * pow(arg, size - i - 1);
    return value;
}

// придаём полиному красивый вид
string view_func(double *lag, int size)
{
    string polynom;
    for (int i = 0; i < size; i++)
        // если коэффиент ненулевой
        if (lag[i] != 0)
        {
            // степень полинома выше 1
            if (i < size - 2)
            {
                if (lag[i] == 1)
                    polynom += "+x^" + to_str(size - i - 1);
                else if (lag[i] == -1)
                    polynom += "x^" + to_str(size - i - 1);
                else if (lag[i] > 0)
                    polynom += "+" + to_str(lag[i]) + "*x^" + to_str(size - i - 1);
                else if (lag[i] < 0)
                    polynom += to_str(lag[i]) + "*x^" + to_str(size - i - 1);
            }
            // равна 1
            else if (i == size - 2)
            {
                if (lag[i] == 1)
                    polynom += "+x";
                else if (lag[i] == -1)
                    polynom += "x";
                else if (lag[i] > 0)
                    polynom += "+" + to_str(lag[i]) + "*x";
                else if (lag[i] < 0)
                    polynom += to_str(lag[i]) + "*x";
            }
            // равна 0
            else if (i == size - 1)
            {
                if (lag[i] > 0)
                    polynom += "+" + to_str(lag[i]);
                else if (lag[i] < 0)
                    polynom += to_str(lag[i]);
            }
        }
    // удаляем лишний плюс, если он есть
    if (polynom[0] == '+')
        polynom.erase(0, 1);
    return polynom;
}

int main()
{
    system("cls");
    string polynom;   // для построения полинома вида a0*x^n+a1*x^n-1+...+an*x^0
    double value = 0; // значение полинома

    int size; // степень многочлена
    cout << "Polinomial degree: ";
    cin >> size;
    size++;

    double arg; // аргумент полинома
    cout << "Polinomial argument: ";
    cin >> arg;

    double *x = new double[size]; // аргументы функций
    bool flag = true;

    while (flag)
    {
        flag = false;
        cout << "\nVector x: ";
        for (int i = 0; i < size; i++)
            cin >> x[i];

        // проверка на одинаковые элементы x
        for (int i = 0; i < size - 1; i++)
            for (int j = i + 1; j < size; j++)
                if (x[i] == x[j])
                    flag = true;

        if (flag)
            cout << "\n||| Try again |||\n";
    }

    double *y = new double[size]; // значения функций
    cout << "\nVector y: ";
    for (int i = 0; i < size; i++)
        cin >> y[i];

    double coeff_sum_el;                  // коэффициент каждого компонента суммы
    int lag_pos;                          // начальная позиция многочлена до сдвига
    double *lag = new double[size];       // многочлен лагранжа
    double *clone_lag = new double[size]; // клон для вычисления произведения многочленов
    double *save_lag = new double[size];  // клон для суммирования компонентов

    calc_func(lag, save_lag, clone_lag, coeff_sum_el, x, y, size, lag_pos);

    if (check_func(lag, size))
        polynom = "0";
    else
    {
        value = value_func(lag, arg, size);
        for (int i = 0; i < size; i++)
            after_point(lag[i]);

        polynom = view_func(lag, size);
    }

    cout << '\n';
    cout << "Lagrange polynom: ";
    cout << polynom << '\n';

    after_point(value);
    cout << "Polinomial value: ";
    cout << value << '\n';

    // очищаем динамическую память
    delete x;
    delete y;
    delete lag;
    delete clone_lag;
    delete save_lag;
}