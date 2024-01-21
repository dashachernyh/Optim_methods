#include "MethodMult_DualLipsh_p.h"

double MethodMult_DualLipsh_p::Funk_multMat(double* y, std::vector<std::vector<int>>& matrix1,
    std::vector<std::vector<int>>& matrix2,
    std::vector<std::vector<std::vector<int>>>& matrix_res,
    int size, int index)
{
    //std::this_thread::sleep_for(std::chrono::milliseconds(100));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix_res[index][i][j] = 0;
            for (int k = 0; k < size; ++k) {
                matrix_res[index][i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    vector<double> val(2);
    val[0] = y[0];
    val[1] = y[1];
    if (task == 0)
        return grish_fam[index_problem]->ComputeFunction(val);
    else
        return gkls_fam[index_problem]->ComputeFunction(val);
}

void MethodMult_DualLipsh_p::Init_Dual_p(int _task, int _index_problem, double* y, double _a, double _b, double _e,
    double _r_l, double _r_g, int _n, int _m, const int _p) {
    r_loc = _r_l;
    r_glob = _r_g;
    a = _a;
    b = _b;
    p = _p;
    eps = _e;
    n = _n;
    m = _m;
    task = _task;
    index_problem = _index_problem;

    Trial_thread current, first, second;

    best_i = 0;
    out_optimal = { 0, 0, 0 };

    first.x = 0;
    first.thread = 0;
    mapd(0, m, y, n, 1);
    ScaleFunc(y[0]);
    InsertScale(y);
    first.z = Funk_mult(y);
    trials_thread.push_back(first);
    out_optimal[0] = y[0];
    out_optimal[1] = y[1];
    out_optimal[2] = first.z;



    /*std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    double x_next;
    std::vector<double> next_point(p - size_t(1));
    for (int j = 0; j < p - 1; j++) {
        x_next = gen() % 100 / (100 * 1.0);
        while (x_next == 1 || x_next == 0)
        {
            x_next = gen() % 100 / (100 * 1.0);
        }
        next_point[j] = x_next;
    }
    sort(next_point.begin(), next_point.end());
*/

    for (int j = 0; j < p - 1; j++) {
        current.x = (j + 1) / double(p); // next_point[j];
        current.thread = 0;
        mapd(current.x, m, y, n, 1);
        InsertScale(y);
        current.z = Funk_mult(y);
        trials_thread.push_back(current);
        if (out_optimal[2] > current.z) {
            out_optimal[0] = y[0];
            out_optimal[1] = y[1];
            out_optimal[2] = current.z;
        }
    }

    second.x = 1;
    second.thread = 0;
    mapd(1, m, y, n, 1);
    InsertScale(y);
    second.z = Funk_mult(y);
    trials_thread.push_back(second);
    if (out_optimal[2] > second.z) {
        out_optimal[0] = y[0];
        out_optimal[1] = y[1];
        out_optimal[2] = second.z;
    }
}

void MethodMult_DualLipsh_p::SolveMult_DualLipsh_p(double* y, std::vector<std::vector<int>>& matrix1,
    std::vector<std::vector<int>>& matrix2,
    std::vector<std::vector<std::vector<int>>>& matrix_res,
    int size)
{
    double M, R_loc, R_glob;
    int itr = 0;
    std::vector<Characteristic_dual> charact(p);

    double z_min = out_optimal[2];

    double power = 1 / double(n);
    double curr_eps = pow(trials_thread[1].x - trials_thread[0].x, power);
    /*out_optimal[2] = z_min; //trials_thread[0].z; //z_min;*/

    double ro = (1 - 1 / r_glob) / (1 - 1 / r_loc) * (1 - 1 / r_glob) / (1 - 1 / r_loc);
    std::vector<double> true_opt = GetTrueOpt();

    /*std::ofstream out1;
    out1.open("Grishagin_p.txt", std::ofstream::ios_base::app);*/

    //while (curr_eps > eps)
    while (itr <= 2000 && (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps)) {

        double d_z = fabs(trials_thread[1].z - trials_thread[0].z);
        double d_x = fabs(trials_thread[1].x - trials_thread[0].x);
        d_x = pow(d_x, power);
        M = d_z / d_x;  //M = 50;

        for (size_t i = 2; i < trials_thread.size(); ++i)
        {
            double max = fabs((trials_thread[i].z - trials_thread[i - 1].z))
                / pow(trials_thread[i].x - trials_thread[i - 1].x, power);
            if (max > M)
                M = max;
        }

        if (M == 0)
            M = 1;

        double R, k;
        for (size_t i = 1; i < trials_thread.size(); ++i)
        {
            k = pow(trials_thread[i].x - trials_thread[i - 1].x, power);

            R_loc = k + (pow((trials_thread[i].z - trials_thread[i - 1].z) / (M * r_loc), 2) / k) -
                2 * (trials_thread[i].z + trials_thread[i - 1].z - 2 * z_min) / (r_loc * M);
            R_glob = k + (pow((trials_thread[i].z - trials_thread[i - 1].z) / (M * r_glob), 2) / k) -
                2 * (trials_thread[i].z + trials_thread[i - 1].z - 2 * z_min) / (r_glob * M);

            if (ro * R_loc > R_glob)
            {
                R = ro * R_loc;
                r = r_loc;
            }
            else
            {
                R = R_glob;
                r = r_glob;
            }

            if (i <= p) {
                charact[i - 1].R = R;
                charact[i - 1].r = r;
                charact[i - 1].pos = i;
                charact[i - 1].right = trials_thread[i];

                if (i == p) {
                    sort(charact.begin(), charact.end());
                }
            }
            else if (R > charact[0].R) {
                charact[0].R = R;
                charact[0].r = r;
                charact[0].pos = i;
                charact[0].right = trials_thread[i];
                sort(charact.begin(), charact.end());
            }
        }

        double delta_z_pos, sgn;
        std::vector<std::pair<Trial_thread, double[2]>> vect_current(p);
        int iter = 0;

#pragma omp parallel num_threads(p) shared(vect_current) private(delta_z_pos, sgn, iter) // proc_bind(spread) 
        {
#pragma omp for schedule(static)
            for (iter = 0; iter < p; ++iter) {
                // ������ ��������� �� ��������, ��� ����� ���������� ��������� �� ����� �����
                delta_z_pos = trials_thread[charact[iter].pos].z - trials_thread[charact[iter].pos - 1].z;

                sgn = 0;
                if (delta_z_pos < 0)
                    sgn = -1;
                if (delta_z_pos > 0)
                    sgn = 1;

                // ������ ����� ��������� ���������� ����� � ����� ������
                vect_current[iter].first.x = (trials_thread[charact[iter].pos].x +
                    trials_thread[charact[iter].pos - 1].x) / 2
                    - sgn / (2 * charact[iter].r) * pow(delta_z_pos / M, n);
                vect_current[iter].first.thread = omp_get_thread_num();
            }
        }

        // �������� ���������� �������� ������ 0
        for (int j = 0; j < p; ++j) {
            mapd(vect_current[j].first.x, m, vect_current[j].second, n, 1);
            InsertScale(vect_current[j].second);
        }
        iter = 0;
#pragma omp parallel num_threads(p) shared(vect_current) private(iter) // proc_bind(spread)
        {
#pragma omp for schedule(static)
            for (iter = 0; iter < p; ++iter) {
                vect_current[iter].first.z = Funk_multMat(vect_current[iter].second,
                    matrix1, matrix2, matrix_res, size, iter);
            }
        }

        for (int j = 0; j < p; ++j) {
            auto it = find(trials_thread.begin(), trials_thread.end(), charact[j].right);
            trials_thread.insert(it, vect_current[j].first);

            if (z_min > vect_current[j].first.z)
            {
                z_min = vect_current[j].first.z;
                best_i = itr;
                out_optimal[0] = vect_current[j].second[0];
                out_optimal[1] = vect_current[j].second[1];
                out_optimal[2] = z_min;
            }
        }
        ++itr;
    }

    //std::cout << "itr = " << itr << std::endl;
    //out1 << trials_thread.size() << std::endl;
    //out1.close();
}