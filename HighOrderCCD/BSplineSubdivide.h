#ifndef BSPLINESUBDIVIDE_H
#define BSPLINESUBDIVIDE_H

#include "HighOrderCCD/Config/Config.h"


PRJ_BEGIN

using namespace std;
class BSplineSubdivide
{
public:
	static int u_k_position(double u, vector<double> u_knot)
    {
        for (int i = 0; i < u_knot.size(); i++)
        {
            if (u == u_knot[0])
            {
                return 0;
            }
            else if (u < u_knot[i])
            {
                return i - 1;
            }
            else if (u == u_knot.back())
            {
                return u_knot.size() - 1;
            }
        }
    }

    static void divide(int poly_degree, vector<vector<Eigen::VectorXd>>& divide_ctrpoint, vector<vector<vector<double>>>& divide_u_knot, std::vector<std::vector<Eigen::MatrixXd>>& divide_basis)
    {
        cout << "µÚ" << subdivision_number + 1 << "´ÎÏ¸·Ö£º" << endl;
        for (int l = 0; l < pow(2, subdivision_number); l++)
        {
            double u = (divide_u_knot[subdivision_number][l][0] + divide_u_knot[subdivision_number][l][divide_u_knot[subdivision_number][l].size()-1]) / 2;
            //u = (divide_u_knot[subdivision_number][l][(int)divide_u_knot[subdivision_number][l].size()/2] + divide_u_knot[subdivision_number][l][(int)(divide_u_knot[subdivision_number][l].size()) / 2-1]) / 2;
            auto BSpline = divide_ctrpoint[subdivision_number][l];
            auto u_knot_temp1 = divide_u_knot[subdivision_number][l];
            auto temp_ctrpoint = BSpline;
            int k = u_k_position(u, u_knot_temp1);
            //cout << "u=" << u << ";k=" << k << endl;
            vector<vector<double>> a;
            a.resize(poly_degree);
            for (int i = 0; i < poly_degree; i++)
            {
                a[i].resize(poly_degree - i);
                for (int j = 0; j < a[i].size(); j++)
                {
                    a[i][j] = (u - u_knot_temp1[j + k - poly_degree + i + 1]) / (u_knot_temp1[j + k + 1] - u_knot_temp1[j + k - poly_degree + i + 1]);
                    //cout << a[i][j] << " ";
                }
                //cout << endl;
            }

            vector<vector<vector<double>>> ctrpoint_new;
            ctrpoint_new.resize(poly_degree + 1);
            for (int i = 0; i < poly_degree + 1; i++)
            {
                ctrpoint_new[i].resize(poly_degree - i + 1);
                if (i == 0)
                {
                    for (int j = 0; j < ctrpoint_new[i].size(); j++)
                    {
                        ctrpoint_new[i][j].resize(3);
                        ctrpoint_new[i][j][0] = BSpline[3 * (k - poly_degree + j) + 0];
                        ctrpoint_new[i][j][1] = BSpline[3 * (k - poly_degree + j) + 1];
                        ctrpoint_new[i][j][2] = BSpline[3 * (k - poly_degree + j) + 2];
                        //cout << "P_" << j + k - poly_degree << "," << i << "=" << endl;
                    }
                }
                else
                {
                    for (int j = 0; j < ctrpoint_new[i].size(); j++)
                    {
                        ctrpoint_new[i][j].resize(3);
                        ctrpoint_new[i][j][0] = (1 - a[i - 1][j]) * ctrpoint_new[i - 1][j][0] + a[i - 1][j] * ctrpoint_new[i - 1][j + 1][0];
                        ctrpoint_new[i][j][1] = (1 - a[i - 1][j]) * ctrpoint_new[i - 1][j][1] + a[i - 1][j] * ctrpoint_new[i - 1][j + 1][1];
                        ctrpoint_new[i][j][2] = (1 - a[i - 1][j]) * ctrpoint_new[i - 1][j][2] + a[i - 1][j] * ctrpoint_new[i - 1][j + 1][2];
                        //cout << "P_" << j + k - poly_degree << "," << i << "=" << (1 - a[i - 1][j])<<"*"<<"P_"<< j + k - poly_degree << "," << i-1 <<"+" << (a[i - 1][j]) << "*" << "P_" << j+1 + k - poly_degree << "," << i - 1 << endl;
                    }
                }
            }
            vector<vector<vector<double>>> ctrpoint_new2;
            ctrpoint_new2.resize(poly_degree + 1);
            for (int i = 0; i < poly_degree + 1; i++)
            {
                ctrpoint_new2[i].resize(poly_degree - i + 1);
                if (i == 0)
                {
                    for (int j = 0; j < ctrpoint_new2[i].size(); j++)
                    {
                        ctrpoint_new2[i][j].resize(3);
                        ctrpoint_new2[i][j][0] = BSpline[3 * (k - poly_degree + j) + BSpline.size() / 2 + 0];
                        ctrpoint_new2[i][j][1] = BSpline[3 * (k - poly_degree + j) + BSpline.size() / 2 + 1];
                        ctrpoint_new2[i][j][2] = BSpline[3 * (k - poly_degree + j) + BSpline.size() / 2 + 2];
                        //cout << "P_" << j + k - poly_degree << "," << i << "=" << endl;
                    }
                }
                else
                {
                    for (int j = 0; j < ctrpoint_new2[i].size(); j++)
                    {
                        ctrpoint_new2[i][j].resize(3);
                        ctrpoint_new2[i][j][0] = (1 - a[i - 1][j]) * ctrpoint_new2[i - 1][j][0] + a[i - 1][j] * ctrpoint_new2[i - 1][j + 1][0];
                        ctrpoint_new2[i][j][1] = (1 - a[i - 1][j]) * ctrpoint_new2[i - 1][j][1] + a[i - 1][j] * ctrpoint_new2[i - 1][j + 1][1];
                        ctrpoint_new2[i][j][2] = (1 - a[i - 1][j]) * ctrpoint_new2[i - 1][j][2] + a[i - 1][j] * ctrpoint_new2[i - 1][j + 1][2];
                        //cout << "P_" << j + k - poly_degree << "," << i << "=" << (1 - a[i - 1][j])<<"*"<<"P_"<< j + k - poly_degree << "," << i-1 <<"+" << (a[i - 1][j]) << "*" << "P_" << j+1 + k - poly_degree << "," << i - 1 << endl;
                    }
                }
            }

            temp_ctrpoint.resize(BSpline.size() + 3 * 2 * (poly_degree + 1));

            divide_ctrpoint[subdivision_number + 1][2 * l].resize(3 * 2 * (k + 1));
            divide_ctrpoint[subdivision_number + 1][2 * l + 1].resize(3 * 2 * (BSpline.size() / 6 - k + poly_degree));
            Eigen::MatrixXd A1(k + 1, BSpline.size() / 6);
            Eigen::MatrixXd A2(BSpline.size() / 6 - k + poly_degree, BSpline.size() / 6);
            A1.setZero();
            A2.setZero();

            for (int i = 0; i < k - poly_degree + 1; i++)
            {
                A1(i, i) = 1;
            }

            A1(k - poly_degree + 1, k - poly_degree) = 1 - a[0][0];
            A1(k - poly_degree + 1, k - poly_degree + 1) = a[0][0];
            A1(k - poly_degree + 2, k - poly_degree) = (1 - a[0][0]) * (1 - a[1][0]);
            A1(k - poly_degree + 2, k - poly_degree + 1) = (a[0][0]) * (1 - a[1][0]) + (1 - a[0][1]) * a[1][0];
            A1(k - poly_degree + 2, k - poly_degree + 2) = (a[0][1]) * a[1][0];
            A1(k - poly_degree + 3, k - poly_degree) = (1 - a[0][0]) * (1 - a[1][0]) * (1 - a[2][0]);
            A1(k - poly_degree + 3, k - poly_degree + 1) = (a[0][0]) * (1 - a[1][0]) * (1 - a[2][0]) + (1 - a[0][1]) * a[1][0] * (1 - a[2][0]) + (1 - a[0][1]) * (1 - a[1][1]) * a[2][0];
            A1(k - poly_degree + 3, k - poly_degree + 2) = (a[0][1]) * a[1][0] * (1 - a[2][0]) + (a[0][1]) * (1 - a[1][1]) * (a[2][0]) + (1 - a[0][2]) * (a[1][1]) * a[2][0];
            A1(k - poly_degree + 3, k - poly_degree + 3) = (a[0][2]) * (a[1][1]) * (a[2][0]);



            divide_basis[subdivision_number + 1][2 * l] = A1 * divide_basis[subdivision_number][l];

            for (int i = BSpline.size() / 6 - k + poly_degree - 1; i >= poly_degree; i--)
            {
                A2(i, i + k - poly_degree) = 1;
            }
            A2(2, k) = a[0][2];
            A2(2, k - 1) = 1 - a[0][2];
            A2(1, k) = a[0][2] * a[1][1];
            A2(1, k - 1) = (1 - a[0][2]) * a[1][1] + (a[0][1]) * (1 - a[1][1]);
            A2(1, k - 2) = (1 - a[0][1]) * (1 - a[1][1]);
            A2(0, k) = a[0][2] * a[1][1] * a[2][0];
            A2(0, k - 1) = A2(1, k - 1) * a[2][0] + a[0][1] * a[1][0] * (1 - a[2][0]);
            A2(0, k - 2) = (1 - a[0][1]) * (1 - a[1][1]) * a[2][0] + (1 - a[0][1]) * a[1][0] * (1 - a[2][0]) + a[0][0] * (1 - a[1][0]) * (1 - a[2][0]);
            A2(0, k - 3) = (1 - a[0][0]) * (1 - a[1][0]) * (1 - a[2][0]);

            divide_basis[subdivision_number + 1][2 * l + 1] = A2 * divide_basis[subdivision_number][l];
            for (int i = 0; i < (k - poly_degree); i++)
            {
                temp_ctrpoint[3 * i] = BSpline[3 * i];
                temp_ctrpoint[3 * i + 1] = BSpline[3 * i + 1];
                temp_ctrpoint[3 * i + 2] = BSpline[3 * i + 2];

                divide_ctrpoint[subdivision_number + 1][2 * l][3 * i] = BSpline[3 * i];
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * i + 1] = BSpline[3 * i + 1];
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * i + 2] = BSpline[3 * i + 2];
                //cout << i << endl;
            }
            //cout << endl;
            for (int i = 0; i < (poly_degree + 1); i++)
            {
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * (k - poly_degree + i)] = ctrpoint_new[i][0][0];
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * (k - poly_degree + i) + 1] = ctrpoint_new[i][0][1];
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * (k - poly_degree + i) + 2] = ctrpoint_new[i][0][2];
                temp_ctrpoint[3 * (k - poly_degree + i)] = ctrpoint_new[i][0][0];
                temp_ctrpoint[3 * (k - poly_degree + i) + 1] = ctrpoint_new[i][0][1];
                temp_ctrpoint[3 * (k - poly_degree + i) + 2] = ctrpoint_new[i][0][2];
                //cout << k - poly_degree + i << endl;
            }
            //cout << endl;
            for (int i = poly_degree; i >= 0; i--)
            {
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (poly_degree - i)] = ctrpoint_new[i][ctrpoint_new[i].size() - 1][0];
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (poly_degree - i) + 1] = ctrpoint_new[i][ctrpoint_new[i].size() - 1][1];
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (poly_degree - i) + 2] = ctrpoint_new[i][ctrpoint_new[i].size() - 1][2];
                //cout << ( poly_degree - i) << endl;
            }
            //cout << endl;
            for (int i = k + poly_degree + 2; i < temp_ctrpoint.size() / 6; i++)
            {
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (i - k - 1)] = BSpline[3 * (i - (poly_degree + 1))];
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (i - k - 1) + 1] = BSpline[3 * (i - (poly_degree + 1)) + 1];
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (i - k - 1) + 2] = BSpline[3 * (i - (poly_degree + 1)) + 2];
                //cout << (i - k - 1) << endl;
            }

            //cout << endl;
            for (int i = 0; i < (k - poly_degree); i++)
            {
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * i + divide_ctrpoint[subdivision_number + 1][2 * l].size() / 2] = BSpline[3 * i + BSpline.size() / 2];
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * i + 1 + divide_ctrpoint[subdivision_number + 1][2 * l].size() / 2] = BSpline[3 * i + 1 + BSpline.size() / 2];
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * i + 2 + divide_ctrpoint[subdivision_number + 1][2 * l].size() / 2] = BSpline[3 * i + 2 + BSpline.size() / 2];
                //cout << i + divide_ctrpoint[subdivision_number + 1][2 * l].size() / 6 << endl;
            }
            //cout << endl;
            for (int i = 0; i < (poly_degree + 1); i++)
            {
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * (k - poly_degree + i) + divide_ctrpoint[subdivision_number + 1][2 * l].size() / 2] = ctrpoint_new2[i][0][0];
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * (k - poly_degree + i) + 1 + divide_ctrpoint[subdivision_number + 1][2 * l].size() / 2] = ctrpoint_new2[i][0][1];
                divide_ctrpoint[subdivision_number + 1][2 * l][3 * (k - poly_degree + i) + 2 + divide_ctrpoint[subdivision_number + 1][2 * l].size() / 2] = ctrpoint_new2[i][0][2];
                //cout << k - poly_degree + i + divide_ctrpoint[subdivision_number + 1][2 * l].size() / 6 << endl;
            }
            //cout << endl;
            for (int i = poly_degree; i >= 0; i--)
            {
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (poly_degree - i) + divide_ctrpoint[subdivision_number + 1][2 * l + 1].size() / 2] = ctrpoint_new2[i][ctrpoint_new2[i].size() - 1][0];
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (poly_degree - i) + 1 + divide_ctrpoint[subdivision_number + 1][2 * l + 1].size() / 2] = ctrpoint_new2[i][ctrpoint_new2[i].size() - 1][1];
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (poly_degree - i) + 2 + divide_ctrpoint[subdivision_number + 1][2 * l + 1].size() / 2] = ctrpoint_new2[i][ctrpoint_new2[i].size() - 1][2];
                //cout << ( poly_degree - i) + divide_ctrpoint[subdivision_number + 1][2 * l + 1].size() / 6 << endl;
            }
            //cout << endl;
            for (int i = k + poly_degree + 2; i < temp_ctrpoint.size() / 6; i++)
            {
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (i - k - 1) + divide_ctrpoint[subdivision_number + 1][2 * l + 1].size() / 2] = BSpline[3 * (i - (poly_degree + 1)) + BSpline.size() / 2];
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (i - k - 1) + 1 + divide_ctrpoint[subdivision_number + 1][2 * l + 1].size() / 2] = BSpline[3 * (i - (poly_degree + 1)) + 1 + BSpline.size() / 2];
                divide_ctrpoint[subdivision_number + 1][2 * l + 1][3 * (i - k - 1) + 2 + divide_ctrpoint[subdivision_number + 1][2 * l + 1].size() / 2] = BSpline[3 * (i - (poly_degree + 1)) + 2 + BSpline.size() / 2];
                //cout << (i-k-1) + divide_ctrpoint[subdivision_number + 1][2 * l + 1].size() / 6 << endl;
            }

            //BSpline = temp_ctrpoint;
            auto knot_temp = u_knot_temp1;
            divide_u_knot[subdivision_number + 1][2 * l].resize(k + poly_degree + 2);
            divide_u_knot[subdivision_number + 1][2 * l + 1].resize(BSpline.size() / 6 - k + 2 * (poly_degree + 1) - 1);
            knot_temp.resize(u_knot_temp1.size() + 2 * (poly_degree + 1));
            for (int i = 0; i <= k; i++)
            {
                knot_temp[i] = u_knot_temp1[i];
                //cout << knot_temp[i] << " ";
            }
            for (int i = 0; i < 2 * (poly_degree + 1); i++)
            {
                knot_temp[i + k + 1] = u;
                //cout << knot_temp[i + k + 1] << " ";
            }
            for (int i = k + 2 * (poly_degree + 1) + 1; i < knot_temp.size(); i++)
            {
                knot_temp[i] = u_knot_temp1[i - (2 * poly_degree + 2)];
                //cout << knot_temp[i] << " ";
            }
            for (int i = 0; i < divide_u_knot[subdivision_number + 1][2 * l].size(); i++)
            {
                divide_u_knot[subdivision_number + 1][2 * l][i] = knot_temp[i];
                //cout << divide_u_knot[subdivision_number + 1][2 * l][i] << "  ";
            }
            //cout << endl;
            for (int i = 0; i < divide_u_knot[subdivision_number + 1][2 * l + 1].size(); i++)
            {
                divide_u_knot[subdivision_number + 1][2 * l + 1][i] = knot_temp[i + divide_u_knot[subdivision_number + 1][2 * l].size()];
                //cout << divide_u_knot[subdivision_number + 1][2 * l+1][i] << "  ";
            }
            //cout << endl;
        }
        subdivision_number++;
    }

};

PRJ_END

#endif