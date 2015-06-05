#include "newGA.hh"

#include <algorithm>
#include <iostream>
#include <set>

using namespace std;

using skeleton newGA;

Problem pbm;

int desde;

// aux -------------------------------------------------

int modulo_positivo(int a, int b) {
	int res = a % b;
	if (res < 0)
		res += b;
	return res;
}

bool mas_cercano_origen(int i, int j) {
	return pbm.duraciones_desde_origen(i) < pbm.duraciones_desde_origen(j);
}

bool mas_cercano(int i, int j) {
	return pbm.duraciones(desde, i) < pbm.duraciones(desde, j);
}

bool mas_cercano_destino(int i, int j) {
	return pbm.duraciones_hasta_destino(i) < pbm.duraciones_hasta_destino(j);
}

int main (int argc, char** argv) {
	ifstream f2(argv[1]);
	if (!f2)
		show_message(12);

	f2 >> pbm;

	Solution sol(pbm);

	for (int contenedor=0; contenedor<pbm.cant_contenedores(); contenedor++) {
		int dia_que_falla;

		while ((dia_que_falla = sol.cumple_demanda(contenedor)) > -1) {
			int primer_dia = dia_que_falla - (int)floor(pbm.tiempo_en_llenarse(contenedor));

			assert(primer_dia < dia_que_falla);

			int dia_menos_lleno = modulo_positivo(primer_dia, pbm.dimension());
			float volumen_dia_menos_lleno = sol.volumen(dia_menos_lleno);

			for (int i=primer_dia+1; i<dia_que_falla; i++) {
				int dia = modulo_positivo(i, pbm.dimension());
				float volumen = sol.volumen(dia);
				if (volumen < volumen_dia_menos_lleno) {
					dia_menos_lleno = dia;
					volumen_dia_menos_lleno = volumen;
				}
			}

			sol.var(dia_menos_lleno).push_back(contenedor);
			sol.set_se_pasa_por(dia_menos_lleno, contenedor, true);

			assert(sol.cumple_capacidad(dia_menos_lleno));
		}
	}

	for (int dia=0; dia<pbm.dimension(); dia++) {
		if (sol.var(dia).size() > 0) {
			set<int> contenedores_disponibles(sol.var(dia).begin(), sol.var(dia).end());

			sol.var(dia).clear();

			int mejor = *min_element(contenedores_disponibles.begin(), contenedores_disponibles.end(), mas_cercano_origen);

			sol.var(dia).push_back(mejor);
			contenedores_disponibles.erase(mejor);

			while (contenedores_disponibles.size() > 1) {
				desde = mejor;
				mejor = *min_element(contenedores_disponibles.begin(), contenedores_disponibles.end(), mas_cercano);

				sol.var(dia).push_back(mejor);
				contenedores_disponibles.erase(mejor);				
			}

			if (contenedores_disponibles.size() > 0) {
				mejor = *min_element(contenedores_disponibles.begin(), contenedores_disponibles.end(), mas_cercano_destino);
				sol.var(dia).push_back(mejor);
			}
		}
	}

	cout << "Solucion: " << sol << endl;
	cout << "Fitness: " << sol.fitness() << endl;
}
