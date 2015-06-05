#ifndef INC_REQ_newGA
#define INC_REQ_newGA
#include "newGA.hh"
#include <math.h>

#include <algorithm>
#include <assert.h>
#include <cstring>
#include <deque>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

skeleton newGA
{

	// aux -------------------------------------------------

	int modulo_positivo(int a, int b) {
		int res = a % b;
		if (res < 0)
			res += b;
		return res;
	}

	// Problem ---------------------------------------------------------------

	Problem::Problem ():_dimension(7)
	{}

	Problem::Problem(const Problem& pbm)
	{
		*this=pbm;
	}

	ostream& operator<< (ostream& os, const Problem& pbm)
	{
		os << endl << endl << "Number of containers: " << pbm._cant_contenedores << endl;
		os << "Time to handle a container: " << pbm._tiempo_recoger << endl;

		return os;
	}

	istream& operator>> (istream& is, Problem& pbm)
	{
		string coordenadas_origen, coordenadas_destino;
		is >> coordenadas_origen >> coordenadas_destino;

		is >> pbm._cant_contenedores;

		pbm._tiempo_en_llenarse = Rarray<float>(pbm._cant_contenedores);

		pbm._duraciones_desde_origen = Rarray<int>(pbm._cant_contenedores);
		pbm._duraciones_hasta_destino = Rarray<int>(pbm._cant_contenedores);

		pbm._duraciones = Rarray<Rarray<int> >(pbm._cant_contenedores);

		pbm._capacidad_contenedor = Rarray<int>(pbm._cant_contenedores);

		for (int i=0; i<pbm._cant_contenedores; i++) {
			pbm._duraciones[i] = Rarray<int>(pbm._cant_contenedores);

			string coordenadas;
			is >> coordenadas;
		}

		for (int i=0; i<pbm._cant_contenedores; i++) {
			is >> pbm._tiempo_en_llenarse[i];
			pbm._tiempo_en_llenarse[i] = min(13.0f, pbm._tiempo_en_llenarse[i]);
		}

		for (int i=0; i<pbm._cant_contenedores; i++)
			is >> pbm._duraciones_desde_origen[i];

		for (int i=0; i<pbm._cant_contenedores; i++)
			is >> pbm._duraciones_hasta_destino[i];

		for (int i=0; i<pbm._cant_contenedores; i++)
			for (int j=0; j<pbm._cant_contenedores; j++)
				if (i != j) {
					int _i, _j;
					is >> _i >> _j >> pbm._duraciones[i][j];
				}

		is >> pbm._tiempo_recoger;

		for (int i=0; i<pbm._cant_contenedores; i++)
			is >> pbm._capacidad_contenedor[i];

		is >> pbm._capacidad_camion;

		return is;
	}

 	Problem& Problem::operator= (const Problem &pbm)
	{
		_duraciones_desde_origen = pbm._duraciones_desde_origen;
		_duraciones_hasta_destino = pbm._duraciones_hasta_destino;
		_dimension = pbm._dimension;
		_cant_contenedores = pbm._cant_contenedores;
		_demanda = pbm._demanda;
		_duraciones = pbm._duraciones;
		_tiempo_recoger = pbm._tiempo_recoger;
		_tiempo_en_llenarse = pbm._tiempo_en_llenarse;
		_capacidad_contenedor = pbm._capacidad_contenedor;
		_capacidad_camion = pbm._capacidad_camion;
		return *this;
	}

	bool Problem::operator== (const Problem& pbm) const
	{
		if (_dimension!=pbm.dimension()) return false; // TODO
		return true;
	}

	bool Problem::operator!= (const Problem& pbm) const
	{
		return !(*this == pbm);
	}

	Direction Problem::direction() const
	{
		return minimize;
	}

	int Problem::dimension() const
	{
		return _dimension;
	}

	int Problem::cant_contenedores() const
	{
		return _cant_contenedores;
	}

	Rarray<float>& Problem::demanda()
	{
		return _demanda;
	}

	int Problem::duraciones(int i, int j) const
	{
		return _duraciones[i][j];
	}

	int Problem::duraciones_desde_origen(int i) const
	{
		return _duraciones_desde_origen[i];
	}

	int Problem::duraciones_hasta_destino(int i) const
	{
		return _duraciones_hasta_destino[i];
	}

	int Problem::tiempo_recoger() const
	{
		return _tiempo_recoger;
	}

	float Problem::tiempo_en_llenarse(int i) const
	{
		return _tiempo_en_llenarse[i];
	}

	int Problem::capacidad_contenedor(int i) const
	{
		return _capacidad_contenedor[i];
	}

	int Problem::capacidad_camion() const
	{
		return _capacidad_camion;
	}

	Problem::~Problem()
	{
	}

	// Solution --------------------------------------------------------------

	Solution::Solution (const Problem& pbm):_pbm(pbm),_var(pbm.dimension()),_se_pasa_por(pbm.dimension())
	{
		for (int dia=0; dia<pbm.dimension(); dia++) {
			_se_pasa_por[dia] = Rarray<bool>(pbm.cant_contenedores());
			
			for (int contenedor=0; contenedor<pbm.cant_contenedores(); contenedor++)
				_se_pasa_por[dia][contenedor] = false;
		}

	}

	const Problem& Solution::pbm() const
	{
		return _pbm;
	}

	Solution::Solution(const Solution& sol):_pbm(sol.pbm()),_var(sol.pbm().dimension()),_se_pasa_por(sol.pbm().dimension())
	{
		*this=sol;
	}

	istream& operator>> (istream& is, Solution& sol)
	{
		for (int dia=0; dia<sol._pbm.dimension(); dia++) {
			int contenedor;
			is >> contenedor;
			while (contenedor != -1) {
				sol._var[dia].push_back(contenedor);
				sol._se_pasa_por[dia][contenedor] = true;
				is >> contenedor;
			}
		}

		return is;
	}

	void imprimir_dia(ostream& os, const Solution& sol, int dia) {
		switch(dia) {
			case 0:
				os << "L\t";
				break;
			case 1:
				os << "M\t";
				break;
			case 2:
				os << "X\t";
				break;
			case 3:
				os << "J\t";
				break;
			case 4:
				os << "V\t";
				break;
			case 5:
				os << "S\t";
				break;
			case 6:
				os << "D\t";
				break;
		}
		if (sol._var[dia].size() > 0) {
			deque<int>::iterator it;
			for (it = sol._var[dia].begin(); it != sol._var[dia].end()-1; it++)
				os << *it << '|';
			os << *it;
		}
	}

	ostream& operator<< (ostream& os, const Solution& sol)
	{
		os << endl << '\t';
		int dia;
		for (dia=0; dia<sol._pbm.dimension()-1; dia++) {
			imprimir_dia(os, sol, dia);
			os << endl << '\t';
		}
		imprimir_dia(os, sol, dia);
		os << endl;

		return os;
	}

	NetStream& operator << (NetStream& ns, const Solution& sol)
	{
		// Escribir en un archivo distinto para cada ns.pid
		//cout << endl << "COMIENZO SALIDA" << endl;
		//cout << sol;
		for (int dia=0; dia<sol._pbm.dimension(); dia++) {

			ns << sol._var[dia].size();
			//int i = 10;

			for (deque<int>::iterator it = sol._var[dia].begin(); it != sol._var[dia].end(); it++) {
				ns << *it;
				//i--;
				//cout << *it << ' ';
			}

			/*while(i>0) {
				ns << 9;
				i--;
			}*/
			//ns << -1;
			//cout << -1 << ' ';
		}

		//cout << sol;

		//cout << endl << "FIN SALIDA" << endl;

		return ns;
	}

	NetStream& operator >> (NetStream& ns, Solution& sol)
	{
		//cout << endl << "COMIENZO ENTRADA" << endl;
		for (int dia=0; dia<sol._pbm.dimension(); dia++) {

			/*int i = 10;

			int contenedor;
			ns >> contenedor;
			i--;
			//cout << contenedor << ' ';
			while (contenedor != 9) {
				sol._var[dia].push_back(contenedor);
				assert(!sol._se_pasa_por[dia][contenedor]);
				sol._se_pasa_por[dia][contenedor] = true;
				ns >> contenedor;
				i--;
				//cout << contenedor << ' ';
			}

			while (i>0) {
				ns >> contenedor;
				i--;
			}*/

			int size;
			ns >> size;

			for (int i=0; i<size; i++) {
				int contenedor;
				ns >> contenedor;

				sol._var[dia].push_back(contenedor);
				assert(!sol._se_pasa_por[dia][contenedor]);
				sol._se_pasa_por[dia][contenedor] = true;
			}

			/*int contenedor;
			ns >> contenedor;
			//cout << contenedor << ' ';
			while (contenedor != -1) {
				sol._var[dia].push_back(contenedor);
				assert(!sol._se_pasa_por[dia][contenedor]);
				sol._se_pasa_por[dia][contenedor] = true;
				ns >> contenedor;
				//cout << contenedor << ' ';
			}*/
		}

		//cout << sol;

		//cout << endl << "FIN ENTRADA" << endl;

		return ns;
	}

 	Solution& Solution::operator= (const Solution &sol)
	{
		_var = sol._var;
		_se_pasa_por = sol._se_pasa_por;
		return *this;
	}

	bool Solution::operator== (const Solution& sol) const
	{
		if (sol.pbm() != _pbm) return false;
		for(int i = 0; i < _var.size(); i++) {
			deque<int>::iterator it1 = _var[i].begin();
			deque<int>::iterator it2 = sol._var[i].begin();

			while (it1 != _var[i].end() && it2 != sol._var[i].end()) {
				if (*it1 != *it2) return false;
			}

			if (it1 != _var[i].end() || it2 != sol._var[i].end()) return false;
		}
		return true;
	}

	bool Solution::operator!= (const Solution& sol) const
	{
		return !(*this == sol);
	}

	int Solution::cumple_demanda(int contenedor) {
		int dias_sin_pasar = 0;

		for (int i=0; i<2; i++) // Hago 2 vueltas para mirarlo circularmente.
			for (int dia=0; dia<_pbm.dimension(); dia++) {
				dias_sin_pasar++;
				if (dias_sin_pasar > _pbm.tiempo_en_llenarse(contenedor))
					return dia*(i+1);
				if (_se_pasa_por[dia][contenedor])
					dias_sin_pasar = 0;
			}

		return -1;
	}

	bool Solution::cumple_demanda() {
		for (int contenedor=0; contenedor<_pbm.cant_contenedores(); contenedor++) {
			if (cumple_demanda(contenedor) > -1)
				return false;
		}
		return true;
	}

	int Solution::dias_desde_ultima_pasada(int contenedor, int dia) {
		int dia_anterior;

		do {
			dia_anterior = modulo_positivo(dia_anterior-1, _pbm.dimension());
		} while (!_se_pasa_por[dia_anterior][contenedor] && dia_anterior != dia);

		if (dia_anterior == dia)
			return _pbm.dimension();
		else
			return modulo_positivo(dia - dia_anterior, _pbm.dimension());
	}

	float Solution::llenado(int contenedor, int dia) {
		return min(1.0f, dias_desde_ultima_pasada(contenedor, dia)/_pbm.tiempo_en_llenarse(contenedor));
	}

	float Solution::volumen(int dia) {
		float volumen = 0.0f;

		for (deque<int>::iterator it=_var[dia].begin(); it!=_var[dia].end(); it++)
			volumen += _pbm.capacidad_contenedor(*it)*llenado(*it, dia);

		return volumen;
	}

	bool Solution::cumple_capacidad(int dia) {
		return volumen(dia) <= _pbm.capacidad_camion();
	}

	bool Solution::cumple_capacidad() {
		for (int dia=0; dia<_pbm.dimension(); dia++) {
			if (!cumple_capacidad(dia))
				return false;
		}
		return true;
	}

	void Solution::initialize() {
		deque<int> dias;
		for (int dia=0; dia<_pbm.dimension(); dia++)
			dias.push_back(dia);

		bool pude_crear = false;
		int intento = 0;

		while (!pude_crear && intento<_pbm.cant_contenedores()*100) {
			bool todo_bien = true;

			for (int contenedor=0; contenedor<_pbm.cant_contenedores(); contenedor++) {
				deque<int> shuffle = dias;
				random_shuffle(shuffle.begin(), shuffle.end());

				while (cumple_demanda(contenedor) > -1) {
					if (shuffle.size() == 0) {
						todo_bien = false;

						// Dejo la solución vacía.
						for (int dia=0; dia<_pbm.dimension(); dia++) {
							_var[dia].clear();

							for (int contenedor=0; contenedor<_pbm.cant_contenedores(); contenedor++)
								_se_pasa_por[dia][contenedor] = false;
						}

						intento++;
						break;
					}

					int dia = shuffle.front();
					shuffle.pop_front();
					_se_pasa_por[dia][contenedor] = true;
					_var[dia].push_back(contenedor);

					if (!cumple_capacidad(dia)) {
						_se_pasa_por[dia][contenedor] = false;
						_var[dia].pop_back();
					}
				}

				if (!todo_bien)
					break;
			}

			if (!todo_bien)
				continue;

			for (int dia=0; dia<_pbm.dimension(); dia++)
				random_shuffle(_var[dia].begin(), _var[dia].end());

			pude_crear = true;
		}

		//cout << intento << "------" << *this << endl;

		assert(pude_crear);
	}

	double Solution::fitness () {
		double fitness = 0.0;

		for (int dia=0;dia<_pbm.dimension();dia++) {
			if (_var[dia].size() > 0) {

				fitness += _pbm.duraciones_desde_origen(_var[dia].front());


				for (deque<int>::iterator it=_var[dia].begin(); it!=_var[dia].end(); it++) {
					deque<int>::iterator sig = it+1;
					if (sig != _var[dia].end()) {
						fitness += _pbm.duraciones(*it, *sig) + _pbm.tiempo_recoger();
					}
				}

				fitness += _pbm.tiempo_recoger() + _pbm.duraciones_hasta_destino(_var[dia].back());
			}
		}

		return fitness;
	}

	char *Solution::to_String() const
	{
		const char separador = '|';

		stringstream ss;

		for (int dia=0; dia<_pbm.dimension(); dia++) {
			ss << _var[dia].size();
			ss << separador;

			for (int j=0; j<_var[dia].size(); j++) {
				ss << _var[dia][j];
				ss << separador;
			}
		}

		string str;

		ss >> str;

		char* res = new char[str.length()+1];

		strcpy (res, str.c_str());

		return res;
	}

	void Solution::to_Solution(char *_string_)
	{
		const char separador = '|';

		string str(_string_);

		stringstream ss;

		ss << str;

		char _separador;

		int dia = 0;

		while (dia<_pbm.dimension()) {
			int cant;
			ss >> cant >> _separador;

			assert(_separador == separador);

			for (int i=0; i<cant; i++) {
				int contenedor;
				ss >> contenedor >> _separador;

				assert(_separador == separador);

				_var[dia].push_back(contenedor);
				_se_pasa_por[dia][contenedor] = true;
			}
			
			dia++;
		}
	}

	unsigned int Solution::size() const
	{
		int res = 0;

		for (int dia=0; dia<_pbm.dimension(); dia++) {
			res += 4 + _var[dia].size()*4; // 4 = 3 + 1. 3 por la máxima cantidad de dígitos de los contenedores, 1 por el pipe ('|').
		}

		return res;
	}


	deque<int>& Solution::var(const int index)
	{
		return _var[index];
	}


	Rarray<deque<int> >& Solution::array_var()
	{
		return _var;
	}

	bool Solution::se_pasa_por(int dia, int contenedor) const
	{
		return _se_pasa_por[dia][contenedor];
	}

	Rarray<bool>& Solution::se_pasa_por(const int dia)
	{
		return _se_pasa_por[dia];
	}

	void Solution::set_se_pasa_por(int dia, int contenedor, bool se_pasa)
	{
		_se_pasa_por[dia][contenedor] = se_pasa;
	}

	void Solution::poner_contenedores_para_cumplir_demanda() {
		for (int contenedor=0; contenedor<_pbm.cant_contenedores(); contenedor++) {
			int dia_que_falla;
			while ((dia_que_falla = cumple_demanda(contenedor)) > -1) {
				int primer_dia = dia_que_falla - (int)floor(_pbm.tiempo_en_llenarse(contenedor)); // primer dia del intervalo que falla.

				assert(primer_dia < dia_que_falla);

				deque<int> dias_candidatos;
				for (int i=primer_dia; i!=dia_que_falla; i++)
					dias_candidatos.push_front(modulo_positivo(i, _pbm.dimension()));

				assert(dias_candidatos.size() > 0);

				random_shuffle(dias_candidatos.begin(), dias_candidatos.end());

				int dia_elegido = dias_candidatos.front();

				int lugar_elegido = rand_int(0, _var[dia_elegido].size());

				_var[dia_elegido].insert(_var[dia_elegido].begin()+lugar_elegido, contenedor);
				_se_pasa_por[dia_elegido][contenedor] = true;
			}
		}
	}

	Solution::~Solution()
	{}

	// UserStatistics -------------------------------------------------------

	UserStatistics::UserStatistics ()
	{}

	ostream& operator<< (ostream& os, const UserStatistics& userstat)
	{
		os << "\n---------------------------------------------------------------" << endl;
		os << "                   STATISTICS OF TRIALS                   	 " << endl;
		os << "------------------------------------------------------------------" << endl;

		for (int i=0;i< userstat.result_trials.size();i++)
		{
			os << endl
			   << userstat.result_trials[i].trial
			   << "\t" << userstat.result_trials[i].best_cost_trial
			   << "\t\t" << userstat.result_trials[i].worst_cost_trial
			   << "\t\t\t" << userstat.result_trials[i].nb_evaluation_best_found_trial
			   << "\t\t\t" << userstat.result_trials[i].nb_iteration_best_found_trial
			   << "\t\t\t" << userstat.result_trials[i].time_best_found_trial
			   << "\t\t" << userstat.result_trials[i].time_spent_trial;
		}
		os << endl << "------------------------------------------------------------------" << endl;
		return os;
	}

	UserStatistics& UserStatistics::operator= (const UserStatistics& userstats)
	{
		result_trials=userstats.result_trials;
		return (*this);
	}

	void UserStatistics::update(const Solver& solver)
	{
		if( (solver.pid()!=0) || (solver.end_trial()!=true)
		  || ((solver.current_iteration()!=solver.setup().nb_evolution_steps())
		       && !terminateQ(solver.pbm(),solver,solver.setup())))
			return;

		struct user_stat *new_stat;

		if ((new_stat=(struct user_stat *)malloc(sizeof(struct user_stat)))==NULL)
			show_message(7);
		new_stat->trial         		 		 = solver.current_trial();
		new_stat->nb_evaluation_best_found_trial = solver.evaluations_best_found_in_trial();
		new_stat->nb_iteration_best_found_trial  = solver.iteration_best_found_in_trial();
		new_stat->worst_cost_trial     		 	 = solver.worst_cost_trial();
		new_stat->best_cost_trial     		 	 = solver.best_cost_trial();
		new_stat->time_best_found_trial		 	 = solver.time_best_found_trial();
		new_stat->time_spent_trial 		 		 = solver.time_spent_trial();

		result_trials.append(*new_stat);
	}

	void UserStatistics::clear()
	{
		result_trials.remove();
	}

	UserStatistics::~UserStatistics()
	{
		result_trials.remove();
	}

// Intra_operator  --------------------------------------------------------------

	Intra_Operator::Intra_Operator(const unsigned int _number_op):_number_operator(_number_op),probability(NULL)
	{}

	unsigned int Intra_Operator::number_operator() const
	{
		return _number_operator;
	}

	Intra_Operator *Intra_Operator::create(const unsigned int _number_op)
	{
		switch (_number_op)
		{
			case 0: return new Crossover;break;
			case 1: return new Mutation();break;
		}
	}

	ostream& operator<< (ostream& os, const Intra_Operator& intra)
	{
		switch (intra.number_operator())
		{
			case 0: os << (Crossover&)intra;break;
			case 1: os << (Mutation&)intra;break;
		}
		return os;
	}

	Intra_Operator::~Intra_Operator()
	{}

//  Crossover:Intra_operator -------------------------------------------------------------

	Crossover::Crossover():Intra_Operator(0)
	{
		probability = new float[1];
	}

	void borrar_elemento_de_deque(deque<int> &dq, int elem) {
		for (deque<int>::iterator it=dq.begin(); it!=dq.end(); it++)
			if (*it == elem) {
				dq.erase(it);
				return;
			}
		assert(false);
	}

	// Intercambiamos el día del corte hasta el corte dentro del día.
	void cruzar_en_dia_de_corte(Solution &padre1, Solution &sol2, int corte_dia_semana, int corte_en_dia_1, int corte_en_dia_2) {
		for (int i=0; i<corte_en_dia_2; i++) {
			sol2.set_se_pasa_por(corte_dia_semana, sol2.var(corte_dia_semana).front(), false);
			sol2.var(corte_dia_semana).pop_front();
		}

		deque<int>::iterator it1 = padre1.var(corte_dia_semana).begin();
		advance(it1, corte_en_dia_1);
		
		for (; it1 != padre1.var(corte_dia_semana).begin(); it1--) {
			if (sol2.se_pasa_por(corte_dia_semana, *it1)) { // Si ya estaba el contenedor que voy a meter.
				if (rand_int(0, 1) == 1) { // Si no entro en el if, elijo no insertar y que quede el que ya estaba.
					borrar_elemento_de_deque(sol2.var(corte_dia_semana), *it1);
					sol2.var(corte_dia_semana).push_front(*it1);
				}
			} else
				sol2.var(corte_dia_semana).push_front(*it1);
			sol2.set_se_pasa_por(corte_dia_semana, *it1, true);
		}

		sol2.poner_contenedores_para_cumplir_demanda();
	}

	void hacer_cumplir_capacidad(Solution &sol, Solution &padre) {
		for (int dia=0; dia<sol.pbm().dimension(); dia++) {
			bool cumple_capacidad;

			if (!(cumple_capacidad = sol.cumple_capacidad(dia))) {

				deque<int> contenedores_del_dia = sol.var(dia);

				while (contenedores_del_dia.size() > 0) {
					Solution copia = sol;

					deque<int> contenedores_disponibles = contenedores_del_dia;
					random_shuffle(contenedores_disponibles.begin(), contenedores_disponibles.end());

					while (!(cumple_capacidad = copia.cumple_capacidad(dia)) && contenedores_disponibles.size() > 0) {
						int contenedor = contenedores_disponibles.front();
						contenedores_disponibles.pop_front();

						borrar_elemento_de_deque(copia.var(dia), contenedor);
						copia.set_se_pasa_por(dia, contenedor, false);

						if (copia.cumple_demanda(contenedor) > -1) {
							borrar_elemento_de_deque(contenedores_del_dia, contenedor);
							break;
						}
					}

					if (cumple_capacidad) {
						sol = copia;
						break;
					} else if (contenedores_disponibles.size() == 0)
						break;
				}
			}

			if (!cumple_capacidad) {
				sol = padre;
				break;
			}
		}
	}

	void Crossover::cross(Solution& sol1, Solution& sol2) const // dadas dos soluciones de la poblacion, las cruza
	{
		Solution padre1 = sol1;
		Solution padre2 = sol2;

		int corte_dia_semana = rand_int(0, sol1.pbm().dimension()-1);

		int corte_en_dia_1, corte_en_dia_2;

		if (sol1.var(corte_dia_semana).size() == 0) // TODO: revisar
			corte_en_dia_1 = 0;
		else
			corte_en_dia_1 = rand_int(0, sol1.var(corte_dia_semana).size()-1);

		if (sol2.var(corte_dia_semana).size() == 0)
			corte_en_dia_2 = 0;
		else
			corte_en_dia_2 = rand_int(0, sol2.var(corte_dia_semana).size()-1);

		// Intercambiamos los días anteriores al corte

		for (int dia=0; dia<corte_dia_semana; dia++) {
			sol2.var(dia) = sol1.var(dia);
			sol2.se_pasa_por(dia) = sol1.se_pasa_por(dia);
			sol1.var(dia) = padre2.var(dia);
			sol1.se_pasa_por(dia) = padre2.se_pasa_por(dia);
		}

		cruzar_en_dia_de_corte(sol1, sol2, corte_dia_semana, corte_en_dia_1, corte_en_dia_2);
		cruzar_en_dia_de_corte(padre2, sol1, corte_dia_semana, corte_en_dia_2, corte_en_dia_1);

		hacer_cumplir_capacidad(sol1, padre1);
		hacer_cumplir_capacidad(sol2, padre2);
	}

	void Crossover::execute(Rarray<Solution*>& sols) const
	{
		for (int i=0;i+1<sols.size();i=i+2)
		 	if (rand01()<=probability[0]) cross(*sols[i],*sols[i+1]);
	}

	ostream& operator<< (ostream& os, const Crossover&  cross)
	{
		 os << "Crossover." << " Probability: "
                    << cross.probability[0]
		    << endl;
		 return os;
	}

	void Crossover::RefreshState(const StateCenter& _sc) const
	{
		_sc.set_contents_state_variable("_crossover_probability",(char *)probability,1,sizeof(float));
	}

	void Crossover::UpdateFromState(const StateCenter& _sc)
	{
		 unsigned long nbytes,length;
		 _sc.get_contents_state_variable("_crossover_probability",(char *)probability,nbytes,length);
	}

	void Crossover::setup(char line[MAX_BUFFER])
	{
		int op;
		sscanf(line," %d %f ",&op,&probability[0]);
		assert(probability[0]>=0);
	}

	Crossover::~Crossover()
	{
		delete [] probability;
	}

	//  Mutation: Sub_operator -------------------------------------------------------------

	Mutation::Mutation():Intra_Operator(1)
	{
		probability = new float[2];
	}

	void Mutation::mutate(Solution& sol) const
	{
		int tiro_moneda = rand_int(0, 1);

		if (tiro_moneda == 0) {
			for (int dia=0; dia<sol.pbm().dimension(); dia++) {
				for (deque<int>::iterator it=sol.var(dia).begin(); it!=sol.var(dia).end(); ) {
					if (rand01()<=probability[1]) {
						tiro_moneda = rand_int(0, 2);
						switch(tiro_moneda) {
							case 0: { // Estas llaves son porque se declara una variable dentro del case, y se verían desde afuera sino pongo las llaves.
								int contenedor_al_azar = rand_int(0, sol.pbm().cant_contenedores()-1);
								if (!sol.se_pasa_por(dia, contenedor_al_azar)) {
									sol.set_se_pasa_por(dia, contenedor_al_azar, true);
									it = sol.var(dia).insert(it, contenedor_al_azar);

									if (!sol.cumple_capacidad(dia)) {
										sol.set_se_pasa_por(dia, contenedor_al_azar, false);
										it = sol.var(dia).erase(it);
									}
								} else
									it++;
								break;
							} case 1: {
								sol.set_se_pasa_por(dia, *it, false);
								if (sol.cumple_demanda(*it) == -1)
									it = sol.var(dia).erase(it);
								else {
									sol.set_se_pasa_por(dia, *it, true);
									it++;
								}
								break;
							}
							default: {
								int otro_dia = rand_int(0, sol.pbm().dimension()-1);
								if (sol.var(otro_dia).size() > 0) {
									int otra_posicion = rand_int(0, sol.var(otro_dia).size()-1);
									deque<int>::iterator it2 = sol.var(otro_dia).begin();
									advance(it2, otra_posicion);

									if (!sol.se_pasa_por(dia, *it2) && !sol.se_pasa_por(otro_dia, *it)) {
										sol.set_se_pasa_por(dia, *it, false);
										sol.set_se_pasa_por(dia, *it2, true);
										sol.set_se_pasa_por(otro_dia, *it2, false);
										sol.set_se_pasa_por(otro_dia, *it, true);
										swap(*it, *it2);
										if (sol.cumple_demanda(*it) > -1 || sol.cumple_demanda(*it2) > -1
												|| !sol.cumple_capacidad(dia) || !sol.cumple_capacidad(otro_dia)) {
											swap(*it, *it2);
											sol.set_se_pasa_por(dia, *it2, false);
											sol.set_se_pasa_por(dia, *it, true);
											sol.set_se_pasa_por(otro_dia, *it, false);
											sol.set_se_pasa_por(otro_dia, *it2, true);
										}
									}
								}
								it++;
							}
						}
					} else {
						it++;
					}
				}
			}
		} else {
			for (int dia=0; dia<sol.pbm().dimension(); dia++) {
				if (rand01()<=probability[1]*sol.var(dia).size()) { // Lo hago con probabilidad tanto como la probabilidad de mutación por el largo de la lista (cantidad de posiciones).
					int a, b;
					a = rand_int(0, sol.var(dia).size()-1);
					b = rand_int(a, sol.var(dia).size()-1);

					deque<int>::iterator ita = sol.var(dia).begin();
					advance(ita, a);

					deque<int>::iterator itb = sol.var(dia).begin();
					advance(itb, b);

					while (a <= b) {
						swap(*ita, *itb); // Como se intercambian dentro del mismo día, no cambia la variable _se_pasa_por, ni se mira la restricción demanda ni capacidad.

						if (rand01()<=probability[1]) {
							int i = rand_int(0, sol.var(dia).size()-1);

							deque<int>::iterator iti = sol.var(dia).begin();
							advance(iti, i);

							swap(*ita, *iti);
						}

						a++;
						b--;
						ita++;
						itb--;
					}
				}
			}
		}
	}

	void Mutation::execute(Rarray<Solution*>& sols) const
	{
		for (int i=0;i<sols.size();i++)
			if(rand01() <= probability[0])	mutate(*sols[i]);
	}

	ostream& operator<< (ostream& os, const Mutation&  mutation)
	{
		os << "Mutation." << " Probability: " << mutation.probability[0]
		   << " Probability1: " << mutation.probability[1]
		   << endl;
		return os;
	}

	void Mutation::setup(char line[MAX_BUFFER])
	{
		int op;
		sscanf(line," %d %f %f ",&op,&probability[0],&probability[1]);
		assert(probability[0]>=0);
		assert(probability[1]>=0);
	}

	void Mutation::RefreshState(const StateCenter& _sc) const
	{
		_sc.set_contents_state_variable("_mutation_probability",(char *)probability,2,sizeof(probability));
	}

	void Mutation::UpdateFromState(const StateCenter& _sc)
	{
		unsigned long nbytes,length;
		_sc.get_contents_state_variable("_mutation_probability",(char *)probability,nbytes,length);
	}

	Mutation::~Mutation()
	{
		delete [] probability;
	}

// StopCondition_1 -------------------------------------------------------------------------------------

	StopCondition_1::StopCondition_1():StopCondition()
	{}

	bool StopCondition_1::EvaluateCondition(const Problem& pbm,const Solver& solver,const SetUpParams& setup)
	{
		//return ((int)solver.best_cost_trial() == pbm.dimension());
		return false;//setup.nb_evolution_steps() == 10000;
	}

	StopCondition_1::~StopCondition_1()
	{}

	//------------------------------------------------------------------------
	// Specific methods ------------------------------------------------------
	//------------------------------------------------------------------------

	bool terminateQ (const Problem& pbm, const Solver& solver,
			 const SetUpParams& setup)
	{
		StopCondition_1 stop;
		return stop.EvaluateCondition(pbm,solver,setup);
	}
}
#endif

