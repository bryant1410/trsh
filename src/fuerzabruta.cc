#include "newGA.hh"

#include <iostream>
#include <deque>

using skeleton newGA;

deque<int> contenedores_totales(const Problem &pbm) {
	deque<int> res;

	for (int contenedor=0; contenedor<pbm.cant_contenedores(); contenedor++)
		res.push_back(contenedor);

	return res;
}

/*void optimo(Solution sol, int dia, deque<int> contenedoresDisponiblesParaElDia, Solution &solOptima) {
	if (sol.fitness() < solOptima.fitness()) {
		if (dia == sol.pbm().dimension()) {
			if (sol.cumple_demanda())
				solOptima = sol;
		} else if (contenedoresDisponiblesParaElDia.size() == 0) {
			optimo(sol, dia+1, contenedores(sol.pbm()), solOptima);
		} else {
			int contenedor = contenedoresDisponiblesParaElDia.front();
			contenedoresDisponiblesParaElDia.pop_front();

			// Pruebo no ponerlo:
			optimo(sol, dia, contenedoresDisponiblesParaElDia, solOptima);

			// Pruebo ponerlo
			sol.var(dia).push_back(contenedor);
			sol.set_se_pasa_por(dia, contenedor, true);
			optimo(sol, dia, contenedoresDisponiblesParaElDia, solOptima);
		}
	}
}*/

void optimo_tsp(Solution sol, deque<int> contenedoresDisponibles, Solution &solOptima) {
	if (sol.fitness() < solOptima.fitness()) {
		if (contenedoresDisponibles.size() == 0)
			solOptima = sol;
		else {
			int contenedor = contenedoresDisponibles.front();
			contenedoresDisponibles.pop_front();

			// Pruebo no ponerlo:
			optimo_tsp(sol, contenedoresDisponibles, solOptima);

			// Pruebo ponerlo:
			sol.var(0).push_back(contenedor);
			sol.set_se_pasa_por(0, contenedor, true);
			optimo_tsp(sol, contenedoresDisponibles, solOptima);
		}
	}
}

void probar_contenedores(Solution **solAux, deque<int> contenedoresDisponibles, deque<int> contenedores, Problem &pbm) {
	if (contenedoresDisponibles.size() == 0) {
		int i = 0;

		for (deque<int>::iterator it=contenedores.begin(); it!=contenedores.end(); it++)
			i += 1 << *it;

		Solution sol(pbm);

		solAux[i] = new Solution(pbm);

		optimo_tsp(sol, contenedores, *solAux[i]);
	} else {
		int contenedor = contenedoresDisponibles.front();
		contenedoresDisponibles.pop_front();

		// Pruebo no ponerlo:
		probar_contenedores(solAux, contenedoresDisponibles, contenedores, pbm);

		// Pruebo ponerlo:
		contenedores.push_back(contenedor);
		probar_contenedores(solAux, contenedoresDisponibles, contenedores, pbm);
	}
}

int main (int argc, char** argv) {
	ifstream f2(argv[1]);
	if (!f2)
		show_message(12);

	Problem pbm;
	f2 >> pbm;

	Solution sol(pbm);
	Solution solOptima(pbm);

	deque<int> contenedores;

	Solution **solAux = new Solution*[1 << 5];

	// Resuelvo TSP primero.
	probar_contenedores(solAux, contenedores_totales(pbm), contenedores, pbm);

	for (int i=(1 << 5)-1; i>=0; i++)
		cout << *solAux[i] << endl;

	// Inizializo el óptimo con la peor solución.
	/*for (int dia=0; dia<pbm.dimension(); dia++)
		for (int contenedor=0; contenedor<pbm.cant_contenedores(); contenedor++) {
			solOptima.var(dia).push_back(contenedor);
			solOptima.set_se_pasa_por(dia, contenedor, true);
		}

	optimo(sol, 0, contenedores(pbm), solOptima);

	cout << "Solucion: " << sol << endl;
	cout << "Fitness: " << sol.fitness() << endl;*/
}
