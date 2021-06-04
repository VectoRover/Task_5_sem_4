#include "read.h"
#include <chrono>
#include <omp.h>


vector<figure*> CRead::ReadFigure(ifstream& File) {
    vector<figure*> figures;
    string str;
    int n = 0;
    double tetr[4][3];
    while (!File.eof())
    {
        getline(File, str);
        n++;
    }
    File.seekg(0, ios_base::beg);
    File.close();
    File.open("datadat.txt");
    figures.resize(n);
    double M = 0;
    double m = 1000000;
    int iii = 0;
    while (getline(File, str))
    {
        vector<string> vec;
        size_t t = 0;
        string res = "";
        while (t < str.size())
        {
            if (str[t] != ' ') {
                res += str[t];
            }
            else {
                vec.push_back(res);
                res = "";
            }
            t++;
        }
        vec.push_back(res);
        if (vec[0] == "sphere") {
            figures[iii] = new sphere(stod(vec[1]), stod(vec[2]), stod(vec[3]), stod(vec[4]));
        }
        int ii = 1;
        if (vec[0] == "tetra") {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 3; ++k) {
                    tetr[j][k] = stod(vec[ii]);
                    ii++;
                }
            }
            figures[iii] = new tetra(tetr);
        }
        iii++;
    }
    for (size_t i = 0; i < figures.size(); i++) {
        if ((*figures[i]).distance(cam) > M) {
            M = (*figures[i]).distance(cam);
        }
    }
    for (size_t i = 0; i < figures.size(); i++) {
        if ((*figures[i]).distance(cam) < m) {
            m = (*figures[i]).distance(cam);
        }
    }
    double clr[3];
    for (size_t i = 0; i < figures.size(); ++i) {
        double dist = (*figures[i]).distance(cam);
        clr[0] = (M - dist) / (M - m) * 255;
        clr[2] = (dist - m) / (M - m) * 255;
        figures[i]->SetColour(clr);
    }
    return figures;
}

CRead::CRead()
{
    SetZero();
}

CRead::CRead(const string data)
{
    ifstream file(data);
    if (!file.is_open()) {
        cout << "Error! Cannot open\n";
        throw 1;
    }
    double tmp = 0;
    string str;
    vector<string> vec;
    while (getline(file, str))
    {

        size_t t = 0;
        string res = "";
        while (t < str.size())
        {
            if (str[t] != ' ') {
                res += str[t];
            }
            else {
                vec.push_back(res);
                res = "";
            }
            t++;
        }
        vec.push_back(res);
    }
    for (int i = 0; i < 3; ++i) { cam[i] = stod(vec[1 + i]); }
    for (int i = 0; i < 3; ++i) { normal[i] = stod(vec[5 + i]); }
    tmp = sqrt(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2));
    for (int i = 0; i < 3; ++i) { normal[i] = normal[i] / tmp; }
    for (int i = 0; i < 3; ++i) { up[i] = stod(vec[9 + i]); }
    tmp = sqrt(pow(up[0], 2) + pow(up[1], 2) + pow(up[2], 2));
    for (int i = 0; i < 3; ++i) { up[i] = up[i] / tmp; }
    screen = stod(vec[13]);
    depth = stod(vec[15]);
    alpha = stod(vec[17]);
    wp = stoi(vec[19]);
    hp = stoi(vec[21]);
    for (int i = 0; i < 3; ++i) { light[i] = stod(vec[23 + i]); }
    height = 2 * screen * tan((alpha / 2) * M_PI / 180);
    pix = height / hp;
    width = pix * wp;
    down[0] = up[1] * normal[2] - up[2] * normal[1];
    down[1] = up[2] * normal[0] - up[0] * normal[2];
    down[2] = up[0] * normal[1] - up[1] * normal[0];
    for (int i = 0; i < 3; ++i) {
        topleft[i] = cam[i] + screen * normal[i] + height / 2 * up[i] + width / 2 * down[i];
    }
    PIX = new double** [hp];
    for (int i = 0; i < hp; i++) {
        PIX[i] = new double* [wp];
        for (int j = 0; j < wp; j++) {
            PIX[i][j] = new double[3];
        }
    }
    for (int i = 0; i < 3; ++i) {
        PIX[0][0][i] = topleft[i] - down[i] * pix - up[i] * pix;
    }
    for (int m = 1; m < hp; ++m) {
        for (int i = 0; i < 3; ++i) {
            PIX[m][0][i] = PIX[m - 1][0][i] - up[i] * pix;
        }
    }
    for (int m = 0; m < hp; ++m) {
        for (int n = 1; n < wp; ++n) {
            for (int i = 0; i < 3; ++i) {
                PIX[m][n][i] = PIX[m][n - 1][i] - down[i] * pix;
            }
        }
    }
    file.close();
}

CRead::~CRead()
{
    for (int i = 0; i < hp; i++) {
        for (int j = 0; j < wp; j++) {
            delete[] PIX[i][j];
        }
        delete[] PIX[i];
    }
    delete[] PIX;

}

void CRead::SetZero()
{
    alpha = 0;
    screen = 0;
    wp = 0;
    hp = 0;
    height = 0;
    width = 0;
    depth = 0;
    pix = 0;
    for (int i = 0; i < 3; ++i) {
        cam[i] = 0;
        normal[i] = 0;
        up[i] = 0;
        light[i] = 0;
        down[i] = 0;
        topleft[i] = 0;
    }
    PIX = nullptr;
}

bool In(double A[3], double B[3], double C[3], double D[3]);

CImg<unsigned char> CRead::Image(vector<figure*> figures)
{
    CImg<unsigned char> img(wp, hp, 1, 3);
    img.fill(0);
    double a[3];
    double c[3];
    double n[3];
    double nl[3];
    double D = 0;
    double t = 0;
    double f[3];
    double tmpclr[3];
    double tmp = 0;
    double* clr;
    double A, B, C;
    double PQRS[4][3];
    double x0 = 0;
    double y0 = 0;
    double z0 = 0;

    auto start = std::chrono::system_clock::now();
#pragma omp parallel for private(x0, y0, z0, clr, tmpclr, a, c, n, nl, D, f, tmp, A, B, C, PQRS, t)
    for (int i = 0; i < hp; ++i) {
        for (int j = 0; j < wp; ++j) {
            t = 0;
            for (size_t m = 0; m < figures.size(); ++m) {
                x0 = (*(figures[m])).XYZ[0];
                y0 = (*(figures[m])).XYZ[1];
                z0 = (*(figures[m])).XYZ[2];
                for (int k = 0; k < 3; ++k) {
                    a[k] = PIX[i][j][k];
                    c[k] = PIX[i][j][k] - cam[k];
                }
                tmp = sqrt(pow(c[0], 2) + pow(c[1], 2) + pow(c[2], 2));
                for (int k = 0; k < 3; ++k) {
                    c[k] = c[k] / tmp;
                }
                sphere* sph = dynamic_cast<sphere*>(figures[m]);
                if (sphere* sph = dynamic_cast<sphere*>(figures[m])) {
                    clr = sph->RGB;
                    D = pow((c[0] * cam[0] + c[1] * cam[1] + c[2] * cam[2] - c[0] * x0 - c[1] * y0 - c[2] * z0), 2) - (pow(c[0], 2) + pow(c[1], 2) + pow(c[2], 2)) * (pow(cam[0], 2) + pow(cam[1], 2) + pow(cam[2], 2) + pow(x0, 2) + pow(y0, 2) + pow(z0, 2) - 2 * (cam[0] * x0 + cam[1] * y0 + cam[2] * z0) - pow(sph->Rad, 2));
                    if (D < 0) {
                    }
                    else {
                        tmp = (-(c[0] * cam[0] + c[1] * cam[1] + c[2] * cam[2] - c[0] * x0 - c[1] * y0 - c[2] * z0) - sqrt(D)) / (pow(c[0], 2) + pow(c[1], 2) + pow(c[2], 2));
                        if ((tmp < t || t == 0) && tmp <= depth) {
                            t = tmp;
                            for (int k = 0; k < 3; ++k) {
                                f[k] = c[k] * t + cam[k];
                            }
                            tmp = sqrt(pow(f[0] - x0, 2) + pow(f[1] - y0, 2) + pow(f[2] - z0, 2));
                            n[0] = (f[0] - x0) / tmp;
                            n[1] = (f[1] - y0) / tmp;
                            n[2] = (f[2] - z0) / tmp;
                            for (int k = 0; k < 3; ++k) {
                                nl[k] = light[k] - f[k];
                            }
                            tmp = sqrt(pow(nl[0], 2) + pow(nl[1], 2) + pow(nl[2], 2));
                            for (int k = 0; k < 3; ++k) {
                                nl[k] = nl[k] / tmp;
                            }
                            tmp = n[0] * nl[0] + n[1] * nl[1] + n[2] * nl[2];
                            if (tmp < 0) { tmp = 0; }
                            for (int k = 0; k < 3; ++k) {
                                tmpclr[k] = clr[k] * tmp;
                            }
                            img.draw_point(j, i, tmpclr);
                        }
                    }
                }
                if (tetra* tet = dynamic_cast<tetra*>(figures[m])) {
                    clr = tet->RGB;
                    int i1 = 1;
                    int i2 = 2;
                    int i3 = 3;
                    for (int i0 = 0; i0 < 4; i0++) {
                        for (int p = 0; p < 3; ++p) { PQRS[0][p] = tet->a[i0][p]; }
                        for (int p = 0; p < 3; ++p) { PQRS[1][p] = tet->a[i1][p]; }
                        for (int p = 0; p < 3; ++p) { PQRS[2][p] = tet->a[i2][p]; }
                        A = (PQRS[1][1] - PQRS[0][1]) * (PQRS[2][2] - PQRS[0][2]) - (PQRS[2][1] - PQRS[0][1]) * (PQRS[1][2] - PQRS[0][2]);
                        B = (PQRS[2][0] - PQRS[0][0]) * (PQRS[1][2] - PQRS[0][2]) - (PQRS[1][0] - PQRS[0][0]) * (PQRS[2][2] - PQRS[0][2]);
                        C = (PQRS[1][0] - PQRS[0][0]) * (PQRS[2][1] - PQRS[0][1]) - (PQRS[2][0] - PQRS[0][0]) * (PQRS[1][1] - PQRS[0][1]);
                        D = -PQRS[0][0] * A - PQRS[0][1] * B - PQRS[0][2] * C;
                        if (A * c[0] + B * c[1] + C * c[2] != 0) {
                            tmp = (-(A * cam[0] + B * cam[1] + C * cam[2] + D)) / (A * c[0] + B * c[1] + C * c[2]);
                            if ((tmp < t || t == 0) && tmp <= depth) {
                                for (int k = 0; k < 3; ++k) { f[k] = c[k] * tmp + cam[k]; }
                                if (In(PQRS[0], PQRS[1], PQRS[2], f)) {
                                    t = tmp;
                                    for (int i = 0; i < 3; ++i) {
                                        n[0] = A;
                                        n[1] = B;
                                        n[2] = C;
                                    }
                                    if ((A * (tet->a[i3][0] - f[0]) + B * (tet->a[i3][1] - f[1]) + C * (tet->a[i3][2] - f[2])) < 0) {
                                    }
                                    else {
                                        for (int i = 0; i < 3; ++i) {
                                            n[i] = -n[i];
                                        }
                                    }
                                    tmp = sqrt(pow(n[0], 2) + pow(n[1], 2) + pow(n[2], 2));
                                    for (int i = 0; i < 3; ++i) { n[i] = n[i] / tmp; }
                                    for (int k = 0; k < 3; ++k) { nl[k] = light[k] - f[k]; }
                                    tmp = sqrt(pow(nl[0], 2) + pow(nl[1], 2) + pow(nl[2], 2));
                                    for (int k = 0; k < 3; ++k) { nl[k] = nl[k] / tmp; }
                                    tmp = n[0] * nl[0] + n[1] * nl[1] + n[2] * nl[2];
                                    if (tmp < 0) { tmp = 0; }
                                    for (int k = 0; k < 3; ++k) { tmpclr[k] = clr[k] * tmp; }
                                    img.draw_point(j, i, tmpclr);

                                }
                            }
                        }
                        i1++;
                        i2++;
                        i3++;
                        if (i1 == 4) { i1 = 0; }
                        if (i2 == 4) { i2 = 0; }
                        if (i3 == 4) { i3 = 0; }
                    }
                }
            }
            t = 0;
        }
    }
    auto end = std::chrono::system_clock::now();
    int elapsed_ms = static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "Addition operator runtime is " << elapsed_ms << " ms\n";
    for (size_t i = 0; i < figures.size(); ++i) { delete figures[i]; }
    return img;
}





bool In(double a1[3], double a2[3], double a3[3], double d[3]) {
    double a[3];
    double b[3];
    for (int i = 0; i < 3; ++i) { a[i] = a2[i] - a1[i]; }
    for (int i = 0; i < 3; ++i) { b[i] = d[i] - a1[i]; }
    double t[3];
    double c[3];
    t[0] = a[1] * b[2] - a[2] * b[1];
    t[1] = a[2] * b[0] - a[0] * b[2];
    t[2] = a[0] * b[1] - a[1] * b[0];
    if (t[0] != 0 || t[1] != 0 || t[2] != 0) {
        double tmp;
        tmp = sqrt(pow(t[0], 2) + pow(t[1], 2) + pow(t[2], 2));
        for (int i = 0; i < 3; ++i) { t[i] = t[i] / tmp; }
        for (int i = 0; i < 3; ++i) { a[i] = a3[i] - a2[i]; }
        for (int i = 0; i < 3; ++i) { b[i] = d[i] - a2[i]; }
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
        if (c[0] != 0 || c[1] != 0 || c[2] != 0) {
            tmp = sqrt(pow(c[0], 2) + pow(c[1], 2) + pow(c[2], 2));
            for (int i = 0; i < 3; ++i) { c[i] = c[i] / tmp; }
            for (int i = 0; i < 3; ++i) { if (fabs(t[i] - c[i]) >= 0.000001) { return false; } }
            for (int i = 0; i < 3; ++i) { a[i] = a1[i] - a3[i]; }
            for (int i = 0; i < 3; ++i) { b[i] = d[i] - a3[i]; }
            c[0] = a[1] * b[2] - a[2] * b[1];
            c[1] = a[2] * b[0] - a[0] * b[2];
            c[2] = a[0] * b[1] - a[1] * b[0];
            if (c[0] != 0 || c[1] != 0 || c[2] != 0) {
                tmp = sqrt(pow(c[0], 2) + pow(c[1], 2) + pow(c[2], 2));
                for (int i = 0; i < 3; ++i) { c[i] = c[i] / tmp; }
                for (int i = 0; i < 3; ++i) { if (fabs(t[i] - c[i]) >= 0.000001) { return false; } }
                return true;
            }
            else return false;
        }
        else return false;
    }
    else return false;
}
