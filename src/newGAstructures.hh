#ifndef INC_newGA_mallba_hh
#define INC_newGA_mallba_hh


#include <iostream>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <Rlist.h>
#include <Rarray.h>
#include <Messages.h>
#include <mallba.hh>
#include <States.hh>
#include <random.hh>
#include <time.hh>
#include <netstream.hh>
#include <assert.h>

using namespace std;

#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_

struct individual // index of a individual in the population and its fitness
{
	int    index;
	double fitness;
	double sel_parameter;
	bool   change;
};

/*int lessF(const struct individual &i1,const  struct individual &i2)
{
	return i1.fitness < i2.fitness;
}

int lessS(const struct individual &i1,const  struct individual &i2)
{
	return i1.sel_parameter < i2.sel_parameter;
}

int greaterF(const struct individual &i1,const  struct individual &i2)
{
	return i1.fitness > i2.fitness;
}

int greaterS(const struct individual &i1,const  struct individual &i2)
{
	return i1.sel_parameter > i2.sel_parameter;
}*/
#endif
#endif
