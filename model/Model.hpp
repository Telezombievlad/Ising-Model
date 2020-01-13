// No Copyright. Vladislav Aleinik 2019
#ifndef POTTS_MODEL_MODEL_HPP_INCLUDED
#define POTTS_MODEL_MODEL_HPP_INCLUDED

#include "StateGraph.hpp"

#include <random>
#include <cstdlib>
#include <cmath>
#include <assert.h>

struct LatticePoint
{
	int stateX;
	int stateY;

	long double interactivityL; // x-1
	long double interactivityR; // x+1
	long double interactivityU; // y-1
	long double interactivityD;	// y+1
};

struct Lattice
{
	int sizeX, sizeY;
	FibonacciSpinStateGraph stateGraph;
	LatticePoint* points;

	Lattice(int latticeSizeX, int latticeSizeY,
			long double (*getInteractivity)(int, int),
	        int         (*getStateX)       (int, int),
            int         (*getStateY)       (int, int));

	~Lattice();

	inline LatticePoint& get(int x, int y);

	inline void metropolisStep();
	inline void metropolisSweep(char thread);
	void metropolisStepAt(int alteredX, int alteredY, int randomNum);

	double calculateMagnetization();
	double calculateEnergy();
};

Lattice::Lattice(int latticeSizeX, int latticeSizeY,
                 long double (*getInteractivity)(int, int),
                 int         (*getStateX)       (int, int),
                 int         (*getStateY)       (int, int)) : 
	sizeX (latticeSizeX),
	sizeY (latticeSizeY),
	stateGraph (),
	points (new LatticePoint[latticeSizeX * latticeSizeY])
{
	for (int x = 0; x < sizeX; ++x)
	{
		for (int y = 0; y < sizeY; ++y)
		{
			points[x * sizeY + y].stateX         = getStateX(x, y);
			points[x * sizeY + y].stateY         = getStateY(x, y);
			points[x * sizeY + y].interactivityR = getInteractivity(x, y);
			points[x * sizeY + y].interactivityD = getInteractivity(x, y);

			int neighborX = (x + 1) % sizeX; 
			int neighborY = (y + 1) % sizeY;

			points[neighborX * sizeY +         y].interactivityL = 
			points[        x * sizeY +         y].interactivityR;
			points[        x * sizeY + neighborY].interactivityU = 
			points[        x * sizeY +         y].interactivityD;
		}		
	}
}

Lattice::~Lattice()
{
	delete[] points;
}

inline LatticePoint& Lattice::get(int x, int y)
{
	return points[x * sizeY + y];
}

inline void Lattice::metropolisStep()
{
	static std::random_device rd;
	static std::uniform_int_distribution<int> dist{0, 100 * sizeX * sizeY};

	int randomNum = dist(rd);

	int alteredX = randomNum % sizeX;
	randomNum /= sizeX;
	int alteredY = randomNum % sizeY;
	randomNum /= sizeY;

	metropolisStepAt(alteredX, alteredY, randomNum);
}

inline void Lattice::metropolisSweep(char thread)
{
	static std::random_device rd;
	static std::uniform_int_distribution<int> dist{0, 100 * sizeX * sizeY};

	for (size_t iter = 0; iter < mc_iters_per_sample; ++iter)
	{
		int randomNum = dist(rd);

		int alteredX = randomNum % (sizeX/2);
		randomNum /= sizeX;
		int alteredY = randomNum % (sizeY/2);
		randomNum /= sizeY;

		if (thread & 0x01) alteredX += sizeX/2;
		if (thread & 0x02) alteredY += sizeY/2;

		metropolisStepAt(alteredX, alteredY, randomNum);
	}
}
	
void Lattice::metropolisStepAt(int alteredX, int alteredY, int randomNum)
{
	static std::random_device rd;
	static std::mt19937 gen{rd()};
	static std::uniform_real_distribution<double> distribution(0.0, 1.0);

	LatticePoint& alteredPoint = get(alteredX                                      , alteredY                                     );
	LatticePoint  neighbourL   = get((alteredX == 0)? (sizeX - 1) : (alteredX - 1) , alteredY                                     );
	LatticePoint  neighbourR   = get((alteredX + 1) % sizeX                        , alteredY                                     );
	LatticePoint  neighbourU   = get(alteredX                                      , (alteredY == 0)? (sizeY - 1) : (alteredY - 1));
	LatticePoint  neighbourD   = get(alteredX                                      , (alteredY + 1) % sizeY                       );

	Vector spinL = stateGraph.get(neighbourL.stateX, neighbourL.stateY) * alteredPoint.interactivityL;
	Vector spinR = stateGraph.get(neighbourR.stateX, neighbourR.stateY) * alteredPoint.interactivityR;
	Vector spinU = stateGraph.get(neighbourU.stateX, neighbourU.stateY) * alteredPoint.interactivityU;
	Vector spinD = stateGraph.get(neighbourD.stateX, neighbourD.stateY) * alteredPoint.interactivityD;
	Vector interactionVector = spinL + spinR + spinU + spinD + externalField;

	int newStateX = alteredPoint.stateX;
	if (randomNum & 0x01)    newStateX = (newStateX + 1) % STATE_GRAPH_SIZE_X;
	else if (newStateX == 0) newStateX = STATE_GRAPH_SIZE_X - 1;
	else                     newStateX = newStateX - 1;

	int newStateY = alteredPoint.stateY;
	if (randomNum & 0x02)    newStateY = (newStateY + 1) % STATE_GRAPH_SIZE_Y;
	else if (newStateY == 0) newStateY = STATE_GRAPH_SIZE_Y - 1;
	else                     newStateY = newStateY - 1;

	Vector curSpin = stateGraph.get(alteredPoint.stateX, alteredPoint.stateY);
	Vector nxtSpin = stateGraph.get(          newStateX,           newStateY);

	double curEnergy = -curSpin.scalar(interactionVector);
	double nxtEnergy = -nxtSpin.scalar(interactionVector);

	if (curEnergy > nxtEnergy)
	{
		alteredPoint.stateX = newStateX;
		alteredPoint.stateY = newStateY;
		return;
	}

	double acceptanceRatio = exp((curEnergy - nxtEnergy) / temperature);
	double toss = distribution(gen);

	if (toss < acceptanceRatio)
	{
		alteredPoint.stateX = newStateX;
		alteredPoint.stateY = newStateY;
	}
}

double Lattice::calculateMagnetization()
{
	double magnetization = 0.0;

	for (int x = 0; x < sizeX; ++x) {
	for (int y = 0; y < sizeY; ++y) {
		LatticePoint curPoint = get(x, y);
		Vector spin = stateGraph.get(curPoint.stateX, curPoint.stateY);

		magnetization += spin.z;
	}}
	
	magnetization /= sizeX*sizeY;

	return magnetization;
}

double Lattice::calculateEnergy()
{
	double energy = 0.0;

	for (int x = 0; x < sizeX; ++x) {
	for (int y = 0; y < sizeY; ++y) {
		LatticePoint& cur        = get(              x,               y);
		LatticePoint& neighbourR = get((x + 1) % sizeX,               y);
		LatticePoint& neighbourD = get(              x, (y + 1) % sizeY);

		Vector spin  = stateGraph.get(       cur.stateX,        cur.stateY);
		Vector spinR = stateGraph.get(neighbourR.stateX, neighbourR.stateY);
		Vector spinD = stateGraph.get(neighbourD.stateX, neighbourD.stateY);		

		energy -= spin.scalar(externalField) +
		          spin.scalar(spinR) * cur.interactivityR +
		          spin.scalar(spinD) * cur.interactivityD;
	}}

	return energy;
}

#endif  // POTTS_MODEL_MODEL_HPP_INCLUDED