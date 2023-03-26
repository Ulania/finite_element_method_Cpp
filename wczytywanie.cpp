#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

#include "wczytywanie.h"

using namespace std;

void Grid::Read(GlobalData* global_data, Grid* grid, Schemat_Calkowania schemat_calkowania) {
    schemat = schemat_calkowania;

    ifstream data_file; //ZMIANA PLIKU
    data_file.open("Test2_4_4_MixGrid.txt");
    // Test1_4_4.txt      Test2_4_4_MixGrid.txt     Test3_31_31_kwadrat.txt

    string bin;

    if (!data_file) {
        cerr << "File can't be opened! " << endl;
        system("PAUSE");
        exit(1);
    }

    data_file >> bin >> global_data->t >> bin >> global_data->dt >> bin >> global_data->k_t >> bin >>
        global_data->alfa >> bin >> global_data->t8 >> bin >> global_data->t0 >> bin >> global_data->ro >> bin >>
        global_data->cp >> bin >> bin >> nN >> bin >> bin >> nE >> bin;

    for (int i = 0; i < nN; i++) {
        ND.push_back(Node());
        data_file >> bin >> ND[i].x >> bin >> ND[i].y;
    }

    data_file >> bin >> bin;

    for (int i = 0; i < nE; i++) {
        EL.push_back(Element());
        data_file >> bin;
        data_file >> EL[i].ID[0] >> bin >> EL[i].ID[1] >> bin >> EL[i].ID[2] >> bin >> EL[i].ID[3];
    }

    data_file >> bin;


    int counter = 0;
    int check = 0;

    for (int i = 0; i < nN; i++) {
        data_file >> check >> bin;
        ND[check - 1].BC = 1;
    }

}

void Grid::Write(GlobalData* global_data, Grid* grid) {
    cout << "-----------------  MES -----------------" << endl;
    cout << "  Data:" << endl;

    cout << "    SimulationTime  " << global_data->t << endl << "    SimulationStepTime  " << global_data->dt << endl << "    Conductivity  " << global_data->k_t << endl << "    Alfa  " <<
        global_data->alfa << endl << "    Tot  " << global_data->t8 << endl << "    InitialTemp  " << global_data->t0 << endl << "    Density  " << global_data->ro << endl << "    SpecificHeat  " << global_data->cp << endl
        << "    nodeNumber  " << nN << endl << "    elementNumber  " << nE
        << endl << "    *Node";

    cout.precision(9);

    for (int i = 0; i < nN; i++) {

        //cout << setprecision(10) << fixed;
        cout << endl << "       ID: " << i + 1 << "  x:" << ND[i].x;
        //cout << setprecision(10) << fixed;
        cout << "   y:" << ND[i].y << "  BC: ";
        cout << ND[i].BC << "  ";
        //cout << setprecision(2) << fixed;
    }

    cout << endl << "     *Element, type=DC2D4" << endl;

    for (int i = 0; i < nE; i++) {
        cout << "         " << i + 1 << ".  ";
        for (int j = 0; j < 4; j++) {
            cout << EL[i].ID[j] << "  ";
        }
        cout << endl;
    }

    cout << endl;
}

void Element4::policz_funkcje_ksztaltu_dNdksi_dNdeta(Schemat_Calkowania schemat) // obliczanie dNdeta dNdksi
{
    this->schemat = schemat;

    // uklad lokalny - obliczanie dN/deta, dN/eksi
    int sch_pnk = schemat.sch_pnk;
    int pnk_calk = sch_pnk * sch_pnk;
    this->dNdeta = new double* [pnk_calk];
    this->dNdksi = new double* [pnk_calk];
    this->funkcje_ksztaltu = new double* [pnk_calk];

    for (int i = 0; i < pnk_calk; i++) {
        dNdeta[i] = new double[4];
        dNdksi[i] = new double[4];
        funkcje_ksztaltu[i] = new double[4];
    }

    for (int i = 0; i < sch_pnk; i++) {
        for (int j = 0; j < sch_pnk; j++) {
            int index = (i * sch_pnk) + j;

            double p_j = schemat.pc[j];
            double p_i = schemat.pc[i];

            //funkcje ksztaltu
            funkcje_ksztaltu[index][0] = 0.25 * (1 - p_j) * (1 - p_i);
            funkcje_ksztaltu[index][1] = 0.25 * (1 + p_j) * (1 - p_i);
            funkcje_ksztaltu[index][2] = 0.25 * (1 + p_j) * (1 + p_i);
            funkcje_ksztaltu[index][3] = 0.25 * (1 - p_j) * (1 + p_i);

            //pochodne po eta -> 1 pc, 2 pc, ...
            dNdeta[index][0] = -0.25 * (1 - p_j);
            dNdeta[index][1] = -0.25 * (1 + p_j);
            dNdeta[index][2] = 0.25 * (1 + p_j);
            dNdeta[index][3] = 0.25 * (1 - p_j);

            //pochodne po ksi
            dNdksi[index][0] = -0.25 * (1 - p_i);
            dNdksi[index][1] = 0.25 * (1 - p_i);
            dNdksi[index][2] = 0.25 * (1 + p_i);
            dNdksi[index][3] = -0.25 * (1 + p_i);
        }
    }

    //tworzenie sciany do HBC
    for (int i = 0; i < 4; i++) {
        sciany[i].eta = new double[sch_pnk]; 
        sciany[i].ksi = new double[sch_pnk]; 
        sciany[i].wagi = new double[sch_pnk]; 
        sciany[i].funkcje_ksztaltu = new double* [sch_pnk]; 

        for (int j = 0; j < sch_pnk; j++) { 
            sciany[i].funkcje_ksztaltu[j] = new double[4];
        }
    }

    //1 sciana punkty calkowania
    for (int i = 0; i < sch_pnk; i++) { // i tu
        //1 sciana punkty calkowani
        sciany[0].ksi[i] = schemat.pc[i];
        sciany[0].eta[i] = -1;
        sciany[0].wagi[i] = schemat.pc[i];
        //2 sciana punkty calkowania
        sciany[1].ksi[i] = 1;
        sciany[1].eta[i] = schemat.pc[i];
        sciany[1].wagi[i] = schemat.pc[i];
        //3sciana punkty calkwoania
        sciany[2].ksi[i] = schemat.pc[i];
        sciany[2].eta[i] = 1;
        sciany[2].wagi[i] = schemat.pc[i];
        //4 sciana punkty calkwoania
        sciany[3].ksi[i] = -1;
        sciany[3].eta[i] = schemat.pc[i];
        sciany[3].wagi[i] = schemat.pc[i];
    }

    //obliczanie na scianie 
    for (int f = 0; f < 4; f++) {
        for (int m = 0; m < sch_pnk; m++) { // i tu
            sciany[f].funkcje_ksztaltu[m][0] = 0.25 * (1.0 - sciany[f].ksi[m]) * (1.0 - sciany[f].eta[m]);
            sciany[f].funkcje_ksztaltu[m][1] = 0.25 * (1.0 + sciany[f].ksi[m]) * (1.0 - sciany[f].eta[m]);
            sciany[f].funkcje_ksztaltu[m][2] = 0.25 * (1.0 + sciany[f].ksi[m]) * (1.0 + sciany[f].eta[m]);
            sciany[f].funkcje_ksztaltu[m][3] = 0.25 * (1.0 - sciany[f].ksi[m]) * (1.0 + sciany[f].eta[m]);
        }
    }


}

void Element4::wypisz_funkcje_ksztaltu_dNdksi_dNdeta() {
    int lp = schemat.sch_pnk * schemat.sch_pnk;

    cout << endl << "tabelka dla funkcji ksztaltu:" << endl;
    for (int j = 0; j < lp; j++) {
        cout << funkcje_ksztaltu[j][0] << "  " << funkcje_ksztaltu[j][1] << "  " << funkcje_ksztaltu[j][2] << "  " << funkcje_ksztaltu[j][3] << endl;
    }

    cout << endl << "tabelka dla dN/dksi:" << endl;
    for (int j = 0; j < lp; j++) {
        cout << dNdksi[j][0] << "  " << dNdksi[j][1] << "  " << dNdksi[j][2] << "  " << dNdksi[j][3] << endl;
    }

    cout << endl << "tabelka dla dN/deta:" << endl;
    for (int j = 0; j < lp; j++) {
        cout << dNdeta[j][0] << "  " << dNdeta[j][1] << "  " << dNdeta[j][2] << "  " << dNdeta[j][3] << endl;
    }
    cout << endl;
}


void Jakobian::obliczanie_jakobianu(Element4 elem4, vector <Node> wezly, vector <Element> elementy, int i, Jakobian* j, int element_index) {
    double x[4];
    double y[4];
    double* pointKsi = elem4.dNdksi[i];
    double* pointEta = elem4.dNdeta[i];

    for (int k = 0; k < 4; k++) {
        x[k] = wezly[elementy[element_index].ID[k] - 1].x;
        y[k] = wezly[elementy[element_index].ID[k] - 1].y;
    }

    //interpolacja
    jakobian[0][0] = pointKsi[0] * x[0] + pointKsi[1] * x[1] + pointKsi[2] * x[2] + pointKsi[3] * x[3];
    jakobian[1][1] = pointEta[0] * y[0] + pointEta[1] * y[1] + pointEta[2] * y[2] + pointEta[3] * y[3];
    jakobian[1][0] = pointKsi[0] * y[0] + pointKsi[1] * y[1] + pointKsi[2] * y[2] + pointKsi[3] * y[3];
    jakobian[0][1] = pointEta[0] * x[0] + pointEta[1] * x[1] + pointEta[2] * x[2] + pointEta[3] * x[3];
    //wyznacznik jakobianu
    wyznacznik = jakobian[0][0] * jakobian[1][1] - jakobian[0][1] * jakobian[1][0];
    //jakobian odwrotny
    jakobian_odwrotny[0][0] = jakobian[1][1] / wyznacznik;
    jakobian_odwrotny[0][1] = -jakobian[1][0] / wyznacznik;
    jakobian_odwrotny[1][0] = -jakobian[0][1] / wyznacznik;
    jakobian_odwrotny[1][1] = jakobian[0][0] / wyznacznik;

    /*// Wypisywanie liczenia jakobianu
    cout << endl << "-----------   " << i << " punkt calkowania   -----------" << endl;
    cout << "jakobian" << endl;
    cout << jakobian[0][0] << "  ";
    cout << jakobian[0][1] << endl;
    cout << jakobian[1][0] << "  ";
    cout << jakobian[1][1] << endl << endl;
    //wyznacznik jakobianu
    cout << "wyznacznik jakobianu:" << endl;
    cout << wyznacznik << endl << endl;
    //jakobian odwrotny
    cout << "jakobian odwrotny" << endl;
    cout << jakobian_odwrotny[0][0] << "  ";
    cout << jakobian_odwrotny[0][1] << endl;
    cout << jakobian_odwrotny[1][0] << "   ";
    cout << jakobian_odwrotny[1][1] << endl << endl;*/
}

void Grid::stworz_H() {
    for (int i = 0; i < nE; i++) {
        EL[i].H = new double* [4];

        for (int j = 0; j < 4; j++) {
            EL[i].H[j] = new double[4];
            for (int k = 0; k < 4; k++) {
                EL[i].H[j][k] = 0;
            }
        }
    }
};

void Grid::oblicz_dN_dy(double* dN_dy, Jakobian jakobian, Element4 element, int j) {
    for (int i = 0; i < 4; i++) {
        dN_dy[i] = jakobian.jakobian_odwrotny[1][1] * element.dNdeta[j][i] + jakobian.jakobian_odwrotny[1][0] * element.dNdksi[j][i];
        //cout << dN_dy[i] << "  ";
    }
    //cout << endl;
}

void Grid::oblicz_dN_dx(double* dN_dx, Jakobian jakobian, Element4 element, int j) {
    for (int i = 0; i < 4; i++) {
        dN_dx[i] = jakobian.jakobian_odwrotny[0][0] * element.dNdksi[j][i] + jakobian.jakobian_odwrotny[0][1] * element.dNdeta[j][i];
        //cout << dN_dx[i] << "  ";
    }
    //cout << endl;
}

void Grid::macierz_H(Element4 element, GlobalData global_data) {
    int pc = schemat.sch_pnk * schemat.sch_pnk;

    double dN_dx[4];
    double dN_dy[4];

    Jakobian jakobian;
    //liczenie jakobianu
    for (int i = 0; i < nE; i++) {
        for (int j = 0; j < pc; j++) {
            jakobian.obliczanie_jakobianu(element, ND, EL, j, &jakobian, i);
            //cout << endl << "dn_dx";
            oblicz_dN_dx(dN_dx, jakobian, element, j);
            //cout << endl << "dn_dy";
            oblicz_dN_dy(dN_dy, jakobian, element, j);

            for (int c = 0; c < 4; c++) {
                for (int h = 0; h < 4; h++) {
                    EL[i].H[c][h] += (dN_dx[c] * dN_dx[h] + dN_dy[c] * dN_dy[h]) * (schemat.wagi[j / schemat.sch_pnk] * schemat.wagi[j % schemat.sch_pnk]) * global_data.k_t * jakobian.wyznacznik; // CALKOWANIE
                    //cout << "EL["<<i<<"].H["<<c<<"]["<<h<<"]" << endl;// << " += " << "(" << dN_dx[c] << " * " << dN_dx[h] << " + " << dN_dy[c] << " * " << dN_dy[h] << ")* (" << schemat.wagi[j / schemat.sch_pnk] << " * " << schemat.wagi[j % schemat.sch_pnk] << " * " << global_data.k_t << " * " << jakobian.wyznacznik << endl;
                }
            }

        }
    }
}

void Grid::wypisz_H()
{
    cout << endl;
    for (int n = 0; n < nE; ++n) {
        cout << "MACIERZ H  " << endl;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                cout << EL[n].H[i][j] << "  ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

void Grid::stworz_H_globalne() {
    H = new double* [nN];
    for (int i = 0; i < nN; i++) {
        H[i] = new double[nN];
    }
    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            H[i][j] = 0.0;
        }
    }
}

void Grid::oblicz_H_globalne() {
    for (int i = 0; i < nE; i++) {
        for (int x = 0; x < 4; x++) {
            for (int y = 0; y < 4; y++) {
                H[EL[i].ID[x] - 1][EL[i].ID[y] - 1] += EL[i].H[x][y] + EL[i].HBC[x][y];
            }
        }
    }
}

void Grid::wypisz_H_globalne() {
    cout << endl << "MACIERZ H GLOBALNA " << endl << endl;
    cout << setprecision(2) << fixed;
    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            cout << H[i][j] << "   ";
        }
        cout << endl;
    }
}

void Grid::stworz_Hbc() {
    for (int i = 0; i < nE; i++) {
        EL[i].HBC = new double* [4];
        for (int j = 0; j < 4; j++) {
            EL[i].HBC[j] = new double[4];
            for (int k = 0; k < 4; k++) {
                EL[i].HBC[j][k] = 0;
            }
        }
    }
}

double Grid::dlugosc_L(Node node1, Node node2) { // twierdzenie pitagorasa
    return sqrt((node2.x - node1.x) * (node2.x - node1.x) + (node2.y - node1.y) * (node2.y - node1.y));
}

void Grid::oblicz_Hbc(Element4 element, GlobalData global_data) {
    double wyznacznik_jakobian;
    for (int e = 0; e < nE; e++) {
        // 1 sciana
        if (ND[EL[e].ID[0] - 1].BC == 1 && ND[EL[e].ID[1] - 1].BC == 1) {
            wyznacznik_jakobian = dlugosc_L(ND[EL[e].ID[0] - 1], ND[EL[e].ID[1] - 1]) / 2.0;
            for (int i = 0; i < element.schemat.sch_pnk; i++) {
                for (int x = 0; x < 4; x++) {
                    for (int y = 0; y < 4; y++) {
                        EL[e].HBC[x][y] += element.schemat.wagi[i] * (element.sciany[0].funkcje_ksztaltu[i][x]) * (element.sciany[0].funkcje_ksztaltu[i][y]) * wyznacznik_jakobian * global_data.alfa;
                    }

                }
            }
        }

        // 2 sciana
        if (ND[EL[e].ID[1] - 1].BC == 1 && ND[EL[e].ID[2] - 1].BC == 1) {
            wyznacznik_jakobian = dlugosc_L(ND[EL[e].ID[1] - 1], ND[EL[e].ID[2] - 1]) / 2.0;
            for (int i = 0; i < element.schemat.sch_pnk; i++) {
                for (int x = 0; x < 4; x++) {
                    for (int y = 0; y < 4; y++) {
                        EL[e].HBC[x][y] += element.schemat.wagi[i] * (element.sciany[1].funkcje_ksztaltu[i][x]) * (element.sciany[1].funkcje_ksztaltu[i][y]) * wyznacznik_jakobian * global_data.alfa;
                    }
                }
            }
        }
        // 3 sciana
        if (ND[EL[e].ID[2] - 1].BC == 1 && ND[EL[e].ID[3] - 1].BC == 1) {
            wyznacznik_jakobian = dlugosc_L(ND[EL[e].ID[2] - 1], ND[EL[e].ID[3] - 1]) / 2.0;
            for (int i = 0; i < element.schemat.sch_pnk; i++) {
                for (int x = 0; x < 4; x++) {
                    for (int y = 0; y < 4; y++) {
                        EL[e].HBC[x][y] += element.schemat.wagi[i] * (element.sciany[2].funkcje_ksztaltu[i][x]) * (element.sciany[2].funkcje_ksztaltu[i][y]) * wyznacznik_jakobian * global_data.alfa;
                    }
                }
            }
        }
        // 4 sciana
        if (ND[EL[e].ID[3] - 1].BC == 1 && ND[EL[e].ID[0] - 1].BC == 1) {
            wyznacznik_jakobian = dlugosc_L(ND[EL[e].ID[0] - 1], ND[EL[e].ID[3] - 1]) / 2.0;
            for (int i = 0; i < element.schemat.sch_pnk; i++) {
                for (int x = 0; x < 4; x++) {
                    for (int y = 0; y < 4; y++) {
                        EL[e].HBC[x][y] += element.schemat.wagi[i] * (element.sciany[3].funkcje_ksztaltu[i][x]) * (element.sciany[3].funkcje_ksztaltu[i][y]) * wyznacznik_jakobian * global_data.alfa;
                    }

                }
            }
        }


    }
}

void Grid::wypisz_Hbc() {
    for (int index = 0; index < nE; ++index) {
        cout << endl << " MACIERZ HBC " << endl;
        for (int x = 0; x < 4; x++) {
            for (int y = 0; y < 4; y++) {
                cout << EL[index].HBC[x][y] << "  ";
            }
            cout << endl;
        }
    }
}

void Grid::stworz_P() {
    for (int i = 0; i < nE; i++) {
        EL[i].P = new double[4];
        for (int j = 0; j < 4; j++) {
            EL[i].P[j] = 0.0;
        }
    }
}

void Grid::oblicz_P(Element4 element, GlobalData global_data) {
    double jakobian_wyznacznik;
    for (int p = 0; p < nE; p++) {
        // 1 sciana
        if (ND[EL[p].ID[0] - 1].BC == 1 && ND[EL[p].ID[1] - 1].BC == 1) {
            jakobian_wyznacznik = dlugosc_L(ND[EL[p].ID[0] - 1], ND[EL[p].ID[1] - 1]) / 2.0;
            for (int i = 0; i < element.schemat.sch_pnk; i++) {
                for (int x = 0; x < 4; x++) {
                    EL[p].P[x] += (jakobian_wyznacznik * element.schemat.wagi[i] * (element.sciany[0].funkcje_ksztaltu[i][x]) * global_data.alfa * global_data.t8);
                }
            }
        }
        // 2 sciana
        if (ND[EL[p].ID[1] - 1].BC == 1 && ND[EL[p].ID[2] - 1].BC == 1) {
            jakobian_wyznacznik = dlugosc_L(ND[EL[p].ID[1] - 1], ND[EL[p].ID[2] - 1]) / 2.0;
            for (int i = 0; i < element.schemat.sch_pnk; i++) {
                for (int x = 0; x < 4; x++) {
                    EL[p].P[x] += (jakobian_wyznacznik * element.schemat.wagi[i] * (element.sciany[1].funkcje_ksztaltu[i][x]) * global_data.alfa * global_data.t8);
                }
            }
        }
        // 3 sciana
        if (ND[EL[p].ID[2] - 1].BC == 1 && ND[EL[p].ID[3] - 1].BC == 1) {
            jakobian_wyznacznik = dlugosc_L(ND[EL[p].ID[2] - 1], ND[EL[p].ID[3] - 1]) / 2.0;
            for (int i = 0; i < element.schemat.sch_pnk; i++) {
                for (int x = 0; x < 4; x++) {
                    EL[p].P[x] += (jakobian_wyznacznik * element.schemat.wagi[i] * (element.sciany[2].funkcje_ksztaltu[i][x]) * global_data.alfa * global_data.t8);
                }
            }
        }
        // 4 sciana
        if (ND[EL[p].ID[3] - 1].BC == 1 && ND[EL[p].ID[0] - 1].BC == 1) {
            jakobian_wyznacznik = dlugosc_L(ND[EL[p].ID[0] - 1], ND[EL[p].ID[3] - 1]) / 2.0;
            for (int i = 0; i < element.schemat.sch_pnk; i++) {
                for (int x = 0; x < 4; x++) {
                    EL[p].P[x] += (jakobian_wyznacznik * element.schemat.wagi[i] * (element.sciany[3].funkcje_ksztaltu[i][x]) * global_data.alfa * global_data.t8);
                }
            }
        }
    }
}

void Grid::wypisz_P() {
    cout << endl << endl;
    for (int n = 0; n < nE; ++n) {
        cout << "WEKTOR P " << endl;
        for (int i = 0; i < 4; i++) {
            cout << EL[n].P[i] << "   ";
        }
        cout << endl << endl;
    }
    cout << endl;
}

void Grid::stworz_C() {
    int pc = schemat.sch_pnk * schemat.sch_pnk;
    for (int i = 0; i < nE; i++) {
        EL[i].C = new double* [pc];
        for (int j = 0; j < pc; j++) {
            EL[i].C[j] = new double[pc];
            for (int k = 0; k < pc; k++) {
                EL[i].C[j][k] = 0;
            }
        }
    }
}

void Grid::oblicz_C(Element4 element, GlobalData global_data) {
    int pc = schemat.sch_pnk * schemat.sch_pnk;

    for (int k = 0; k < nE; k++) {
        Jakobian jakobian;

        for (int l = 0; l < pc; l++) {
            jakobian.obliczanie_jakobianu(element, ND, EL, l, &jakobian, k);

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    EL[k].C[i][j] += jakobian.wyznacznik * (schemat.wagi[l / schemat.sch_pnk] * schemat.wagi[l % schemat.sch_pnk]) * element.funkcje_ksztaltu[l][i] * element.funkcje_ksztaltu[l][j] * global_data.cp * global_data.ro;
                    //cout<< "jakobian.wyznacznik * (schemat.wagi["<<l  <<" / "<< schemat.sch_pnk<<"] * schemat.wagi[" << l<<" % "<< schemat.sch_pnk << "])* element.funkcje_ksztaltu[" << l << "][" << i << "] * element.funkcje_ksztaltu[" << l << "][" << j << "] * " << global_data.cp << " * " << global_data.ro << endl;
                }
            }

        }
    }
}

void Grid::wypisz_C() {
    for (int n = 0; n < nE; ++n) {
        cout << endl << " MACIERZ C " << endl;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                cout << EL[n].C[i][j] << "   ";
            }
            cout << endl;
        }
    }
}

void Grid::stworz_P_globalne() {
    P = new double[nN];
    for (int i = 0; i < nN; i++) {
        P[i] = 0;
    }
}

void Grid::oblicz_P_globalne() {
    for (int i = 0; i < nE; i++) {
        for (int x = 0; x < 4; x++) {
            P[EL[i].ID[x] - 1] += EL[i].P[x];
        }
    }
}

void Grid::wypisz_P_globalne() {
    cout << endl << "WEKTOR P GLOBALNY" << endl;
    for (int i = 0; i < nN; i++) {
        cout << P[i] << endl;
    }
}

void Grid::stworz_C_globalne() {
    C = new double* [nN];
    CdT = new double* [nN];
    for (int i = 0; i < nN; i++) {
        C[i] = new double[nN];
        CdT[i] = new double[nN];
    }
    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            C[i][j] = 0.0;
            CdT[i][j] = 0.0;
        }
    }
}

void Grid::oblicz_C_globalne() {
    for (int i = 0; i < nE; i++) {
        for (int x = 0; x < 4; x++) {
            for (int y = 0; y < 4; y++) {
                cout << setprecision(2) << fixed;
                C[EL[i].ID[x] - 1][EL[i].ID[y] - 1] += EL[i].C[x][y];
            }
        }
    }
}

void Grid::wypisz_C_globalne() {
    cout << endl << "MACIERZ C GLOBALNA " << endl << endl;
    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            cout << C[i][j] << "  ";
        }
        cout << endl;
    }
}


void Grid::oblicz_kopie(double** H, double* P) {
    double zmienna;
    for (int i = 0; i < (nN - 1); i++) {
        for (int k = (i + 1); k < nN; k++) {
            zmienna = H[k][i] / H[i][i];
            for (int j = 0; j < nN; j++) {
                H[k][j] -= H[i][j] * zmienna;
            }
            P[k] -= zmienna * P[i];
        }
    }
}

void Grid::Gauss_rozwiazanie(double* P, double* z) {
    double** kopiaH = new double* [nN];
    double* kopiaP = new double[nN];
    for (int i = 0; i < nN; i++) {
        kopiaH[i] = new double[nN];
        kopiaP[i] = P[i];
        for (int j = 0; j < nN; j++) {
            kopiaH[i][j] = H[i][j];
        }
    }

    oblicz_kopie(kopiaH, kopiaP);

    for (int i = (nN - 1); i >= 0; i--) { // UKLAD ROWNAN
        z[i] = kopiaP[i];
        for (int j = i; j < (nN - 1); j++) {
            z[i] -= kopiaH[i][j + 1] * z[j + 1];
        }
        z[i] /= kopiaH[i][i]; // temperatury na wszytkich wezlach -> wektor
    }
}

void Grid::oblicz_temperature(double* temperatury, double* zmienna, GlobalData global_data) {

    double minimalna;
    double maksymalna;
    cout << endl << endl;
    int j = 0;
    do {
        j++;

        for (int i = 0; i < nN; i++) {
            zmienna[i] = 0;
            for (int j = 0; j < nN; j++) {
                zmienna[i] += temperatury[j] * CdT[i][j];
            }
        }

        for (int i = 0; i < nN; i++) {

            zmienna[i] += P[i]; // otrzymujemy drugi czlon rownania C/delta_tal * t0 + wektor P -> wektor
        }

        Gauss_rozwiazanie(zmienna, temperatury);


        minimalna = temperatury[0];
        maksymalna = temperatury[0];
        for (int i = 1; i < nN; i++) {
            if (minimalna > temperatury[i]) {
                minimalna = temperatury[i];
            }
            if (maksymalna < temperatury[i]) {
                maksymalna = temperatury[i];
            }
        }
        cout << "Czas [s]: " << (j * global_data.dt) << "\t" << "\tMinimalna temperatura\t " << minimalna << "\t\t" << "\t Max temperatura\t " << maksymalna << endl;
    } while (j < (global_data.t / global_data.dt));
}

void Grid::oblicz_CdT(GlobalData global_data) { // dostajemy H * C/delta tal -> pierwszy czlon rownania przez ktory przemnazamy nieznana temp. -> macierz
    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            CdT[i][j] = C[i][j] / global_data.dt;
        }
    }


    for (int i = 0; i < nN; i++) {
        for (int j = 0; j < nN; j++) {
            this->H[i][j] += CdT[i][j];
        }
    }
}

