#include "poisson.h"
#include "solution.h"
#include "preconditioner.h"
#include "GMRES.h"
#include "fstream"
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include "systemSolver.h"
#include "DataBase.h"
#include "UCom.h"
#include <UINsInvterm.h>

SolveMRhs bgx;
SolveMRhs::SolveMRhs()
{
	;
}

SolveMRhs::~SolveMRhs()
{
	;
}

SolveMRhs Rank;
void SolveMRhs::Init()
{
	TempA = ArrayUtils<double>::onetensor(Rank.NUMBER);
	TempIA = ArrayUtils<int>::onetensor(Rank.RANKNUMBER+1);
	TempJA = ArrayUtils<int>::onetensor(Rank.NUMBER);
	TempB = ArrayUtils<double>::twotensor(Rank.RANKNUMBER,Rank.COLNUMBER);
	TempX = ArrayUtils<double>::twotensor(Rank.RANKNUMBER,Rank.COLNUMBER);
}
void SolveMRhs::Deallocate()
{
	ArrayUtils<double>::delonetensor(TempA);
	ArrayUtils<int>::delonetensor(TempIA);
	ArrayUtils<int>::delonetensor(TempJA);
	ArrayUtils<double>::deltwotensor(TempB);
	ArrayUtils<double>::deltwotensor(TempX);
	TempA = NULL;
	TempIA = NULL;
	TempJA = NULL;
	TempB = NULL;
	TempX = NULL;
}
void SolveMRhs::BGMRES()
{
	clock_t start, finish;
	double time;
	start = clock();
	Poisson* A = new Poisson;   // The operator to invert.
	Solution* x = new Solution(Rank.RANKNUMBER);  // The approximation to calculate.
	Solution* b = new Solution(Rank.RANKNUMBER);  // The forcing function for the r.h.s.
	Solution* residual = new Solution(Rank.RANKNUMBER);
	Preconditioner* pre =
		new Preconditioner(Rank.RANKNUMBER);      // The preconditioner for the system.
	int restart = 0;                    // Number of restarts to allow
	int maxIt = ONEFLOW::GetDataValue< int >("GMRESIterStep");
	double tol = ONEFLOW::GetDataValue< double >("GMRESTol");

	/**
	   produce the right-hand sides
	*/
	int i, j;
	for (i = 0; i < Rank.RANKNUMBER; i++)
	{
		for (j = 0; j < Rank.COLNUMBER; j++)
		{
			{
				(*b)(i, j) = Rank.TempB[i][j];
			}
		}
	}

	int result = GMRES(A, x, b, residual, pre, maxIt, restart, tol);

	// Output the solution
	for (int lupe = 0; lupe < Rank.COLNUMBER; lupe++)
	{
		for (int innerlupe = 0; innerlupe < Rank.RANKNUMBER; innerlupe++)
		{
			Rank.TempX[innerlupe][lupe] = (*x)(innerlupe, lupe);
		}
	}

	delete A;
	delete x;
	delete b;
	delete residual;
	delete pre;
	A = NULL;
	x = NULL;
	b = NULL;
	residual = NULL;
	pre = NULL;


#define SOLUTION
#ifdef SOLUTION
	
#endif
}
