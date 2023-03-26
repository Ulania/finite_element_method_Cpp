#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

#include "wczytywanie.h"

using namespace std;

int main()
{
    //KWADRATURY GAUSSA - calkowanie numeryczne

    double* s2p = new double[2];
    double* s3p = new double[3];
    double* s4p = new double[4];
    double* s2wagi = new double[2];
    double* s3wagi = new double[3];
    double* s4wagi = new double[4];

    //punkty calkowania - wartosc punktu x na osi X
    s2p[0] = -1.0 / sqrt(3.0);
    s2p[1] = 1.0 / sqrt(3.0);
    s3p[0] = -sqrt(3.0 / 5.0);
    s3p[1] = 0.0;
    s3p[2] = sqrt(3.0 / 5.0);
    s4p[0] = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)); // -0.861136;
    s4p[1] = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)); // -0.339981;
    s4p[2] = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));  // 0.339981;
    s4p[3] = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));  // 0.861136;

    //wagi - dlugosc od 0 do konca funkcji
    s2wagi[0] = 1.0;
    s2wagi[1] = 1.0;
    s3wagi[0] = 5.0 / 9.0;
    s3wagi[1] = 8.0 / 9.0;
    s3wagi[2] = 5.0 / 9.0;
    s4wagi[0] = (18.0 + sqrt(30)) / 36; // 0.347855;
    s4wagi[1] = (18.0 - sqrt(30)) / 36; // 0.652145;
    s4wagi[2] = (18.0 - sqrt(30)) / 36; // 0.652145;
    s4wagi[3] = (18.0 + sqrt(30)) / 36; // 0.347855;

    //schematy calkowania dla 2, 3, 4 punktow
    Schemat_Calkowania sch_2p = { 2, s2p, s2wagi };
    Schemat_Calkowania sch_3p = { 3, s3p, s3wagi };
    Schemat_Calkowania sch_4p = { 4, s4p, s4wagi };

    Schemat_Calkowania wybrany_schemat = sch_2p; //WYBIERZ SCHEMAT

    Element4 elem4;
    elem4.policz_funkcje_ksztaltu_dNdksi_dNdeta(wybrany_schemat);
    elem4.wypisz_funkcje_ksztaltu_dNdksi_dNdeta();

    GlobalData* global_data = new GlobalData;
    Grid* grid = new Grid;    

    grid->Read(global_data, grid, wybrany_schemat);
    grid->Write(global_data, grid);

    grid->stworz_H();
    grid->macierz_H(elem4, *global_data);
    grid->wypisz_H();

    grid->stworz_Hbc();
    grid->oblicz_Hbc(elem4, *global_data);
    grid->wypisz_Hbc();

    grid->stworz_P();
    grid->oblicz_P(elem4, *global_data);
    grid->wypisz_P();

    grid->stworz_C();
    grid->oblicz_C(elem4, *global_data);
    grid->wypisz_C();

    grid->stworz_H_globalne();
    grid->oblicz_H_globalne(); // AGREGACJA H i HBC
    grid->wypisz_H_globalne();

    grid->stworz_P_globalne();
    grid->oblicz_P_globalne();
    grid->wypisz_P_globalne();

    grid->stworz_C_globalne();
    grid->oblicz_C_globalne();
    grid->wypisz_C_globalne();
        

    double* temperatury = new double[grid->nN];
    double* zmienna = new double[grid->nN]; // Cdt+P
    for (int i = 0; i < grid->nN; i++) {
        temperatury[i] = global_data->t0;
        zmienna[i] = 0.0;
    }
    grid->oblicz_CdT(*global_data); // dostajemy H * C/delta tal -> pierwszy czlon rownania przez ktory przemnazamy nieznana temp. -> macierz
    grid->oblicz_temperature(temperatury, zmienna, *global_data); 

}