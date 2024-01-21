#include <iomanip>
#include "MethodMult_DualLipsh_p.h" // include MethodMult_p

bool Checked_method(double val_meth, double val_true, double eps)
{
    if (fabs(val_true - val_meth) <= eps)
        return true;
    else
        return false;
}

bool Checked_method_mult(int n, std::vector<double> val_meth, std::vector<double> val_true, double eps)
{
    if (n == 2) {
        if (fabs(val_true[0] - val_meth[0]) <= eps && fabs(val_true[1] - val_meth[1]) <= eps)
            return true;
        else
            return false;
    }
    else {
        if (fabs(val_true[0] - val_meth[0]) <= eps)
            return true;
        else
            return false;
    }
}


int main(int argc, char** argv)
{

    double E = 0.01, r = 6.0;
    int n = 2, m = 10;
    int count_true = 0;
    double average = 0;
    int task = 1; // 0 - Grishagin, 1 - GKLS
    double r_loc = 1.5;
    double* y = new double[n];

    int p;
    int size = 180;
    std::cout << "Enter p\n";
    std::cin >> p;
    std::vector<std::vector<int>> matrix1(size, std::vector<int>(size));
    std::vector<std::vector<int>> matrix2(size, std::vector<int>(size));
    std::vector<std::vector<std::vector<int>>> matrix_res(p,
        std::vector<std::vector<int>>(size, std::vector<int>(size)));

    std::mt19937 gen(time(NULL));
    std::uniform_int_distribution<> distrib(-100, 100);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix1[i][j] = distrib(gen);
            matrix2[i][j] = distrib(gen);
        }
    }

    double average_time = 0;
    MethodMult_DualLipsh_p met(p);

    std::ofstream out;
    out.open("Result.txt", std::ofstream::ios_base::app);
    out << "Task " << task << "\n";
    out << "Mult_DualLipsh_p " << p << " E = " << E << " r = " << r << " r_loc = " << r_loc << std::endl;


    for (int index = 0; index < 100; index++)
    {
        met.Init_Dual_p(task, index, y, -1, 1, E, r_loc, r, n, m, p);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        met.SolveMult_DualLipsh_p(y, matrix1, matrix2, matrix_res, size);
        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        auto sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        average_time += sec.count();

        if (Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), E))
        {
            average += met.GetBestIndex();
            count_true++;
            out << met.GetBestIndex() << std::endl;
        }
        else
            out << "wrong" << std::endl;
        met.ClearMethod_p();
    }
    out << "Average Time = " << average_time / 100 << "\n";
    out << "count_true " << count_true << " average " << average / count_true << std::endl;
    std::cout << "count_true " << count_true << " average " << average / count_true << std::endl;

    out << "Mult_DualLipsh_p_finish" << std::endl;
    out.close();
    delete[]y;
}
