#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>


using namespace std;

/////////////////   STRUKTURY   /////////////////

struct GlobalData {
    int t;    // czas symulacji
    int dt;   // krok czasowy
    int k_t;  // wspolczynnik przewodnosci cieplnej
    int alfa; // do warunku brzegowego konwekcji
    int t8;   // temp. otoczenia
    int t0;   // temp. poczatkowa
    int ro;   // gestosc
    int cp;   // cieplo wlasciwe
};

struct Schemat_Calkowania {
    int sch_pnk;
    double* pc;
    double* wagi;
};

struct Sciana {
    double* eta;
    double* ksi;
    double* wagi;
    double** funkcje_ksztaltu;
};

struct Node {
    double x = 0;
    double y = 0;
    double t = 0;
    int BC = 0;
};

struct Element {
    int ID[4]; // id wezlow
    double** H;
    double** HBC;
    double* P;
    double** C;
};


struct Element4 {
    Schemat_Calkowania schemat;
    Sciana* sciany = new Sciana[4];
    double** dNdeta;
    double** dNdksi;
    double** funkcje_ksztaltu;

    void policz_funkcje_ksztaltu_dNdksi_dNdeta(Schemat_Calkowania schemat);
    void wypisz_funkcje_ksztaltu_dNdksi_dNdeta();
};

struct Jakobian {
    double jakobian[2][2];
    double jakobian_odwrotny[2][2];
    double wyznacznik;
     
    void obliczanie_jakobianu(Element4 element, vector <Node> wezly, vector <Element> elementy, int i, Jakobian* j, int element_index);
};

struct Grid {
    Schemat_Calkowania schemat;
    int nN = 0; //liczba wezlow
    int nE = 0; // liczba elementow
    std::vector<Node> ND;
    std::vector<Element> EL;
    double** H;
    double* P;
    double** C;
    double** CdT;

    void stworz_H();
    void macierz_H(Element4 element, GlobalData global_data);
    void oblicz_dN_dx(double* dN_dx, Jakobian jakobian, Element4 element, int j);
    void oblicz_dN_dy(double* dN_dy, Jakobian jakobian, Element4 element, int j);
    void wypisz_H();

    void stworz_Hbc();
    double dlugosc_L(Node node1, Node node2);
    void oblicz_Hbc(Element4 elem4, GlobalData global_data);
    void wypisz_Hbc();

    void stworz_P();
    void oblicz_P(Element4 element, GlobalData global_data);
    void wypisz_P();

    void stworz_C();
    void oblicz_C(Element4 element, GlobalData global_data);
    void wypisz_C();

    void stworz_H_globalne();
    void oblicz_H_globalne();
    void wypisz_H_globalne();

    void stworz_P_globalne();
    void oblicz_P_globalne();
    void wypisz_P_globalne();

    void stworz_C_globalne();
    void oblicz_C_globalne();
    void wypisz_C_globalne();


    void oblicz_kopie(double** H, double* P);
    void Gauss_rozwiazanie(double* P, double* z);
    void oblicz_temperature(double* temperatury, double* zmienna, GlobalData global_data);
    void oblicz_CdT(GlobalData global_data);

    void Read(GlobalData* global_data, Grid* grid, Schemat_Calkowania schemat_calkowania);
    void Write(GlobalData* global_data, Grid* grid);
};
