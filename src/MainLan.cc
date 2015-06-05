#include "newGA.hh"

int main (int argc, char** argv)
{
	using skeleton newGA;

	char path[MAX_BUFFER];
	int len;
	int longitud;

	system("clear");

	get_path(argv[0],path);
	len = strlen(path);
	longitud = MAX_BUFFER - len;

	strcat(path,"Config.cfg");
	ifstream f(path);
	if(!f) show_message(10);

	f.getline(&(path[len]),longitud,'\n');
	ifstream f1(path);
	if(!f1)	show_message(11);

	f.getline(&(path[len]),longitud,'\n');
	ifstream f2(path);
	if(!f2) show_message(12);

	Problem pbm;
	f2 >> pbm;

	Operator_Pool pool(pbm);
	SetUpParams cfg(pool);
	f1 >> cfg;


	Solver_Lan solver(pbm,cfg,argc,argv);
	solver.run();

	if (solver.pid()==0)
	{
		solver.show_state();
		cout << "Solution: " << solver.global_best_solution() << " Fitness: " << solver.global_best_solution().fitness();

		f.getline(&(path[len]),longitud,'\n');
	  	ofstream fexit(path);
	  	if(!fexit) show_message(13);
	  	fexit << solver.userstatistics();

		cout << endl << endl << " :( ---------------------- THE END --------------- :) " << endl;
	}
	return(0);
}
