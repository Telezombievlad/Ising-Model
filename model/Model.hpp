// No Copyright. Vladislav Aleinik 2019
#ifndef POTTS_MODEL_MODEL_HPP_INCLUDED
#define POTTS_MODEL_MODEL_HPP_INCLUDED

// MAXIM LOH
// VLADIK SUPER MOLODEC

#include "StateGraph.hpp"

#include <random>
#include <cstdlib>
#include <cmath>
#include <assert.h>

struct LatticePoint
{
	int stateX;
	int stateY;
};

struct Lattice
{
	int sizeX, sizeY, sizeZ;
	FibonacciSpinStateGraph stateGraph;
	LatticePoint* points;

	Lattice(int latticeSizeX, int latticeSizeY, int latticeSizeZ,
	        int (*getStateX) (int, int, int),
            int (*getStateY) (int, int, int));

	~Lattice();

	inline LatticePoint& get(int x, int y, int z);

	inline void metropolisStep();
	void metropolisStepAt(int alteredX, int alteredY, int alteredZ, int randomNum);

	double calculateMagnetization();
	double calculateEnergy();
};

Lattice::Lattice(int latticeSizeX, int latticeSizeY, int latticeSizeZ,
	             int (*getStateX) (int, int, int),
                 int (*getStateY) (int, int, int)) : 
	sizeX (latticeSizeX),
	sizeY (latticeSizeY),
	sizeZ (latticeSizeZ),
	stateGraph (),
	points (new LatticePoint[latticeSizeX * latticeSizeY * latticeSizeZ])
{
	for (int x = 0; x < sizeX; ++x) {
	for (int y = 0; y < sizeY; ++y) {
	for (int z = 0; z < sizeZ; ++z) {
		points[(x * sizeY + y) * sizeZ + z].stateX = getStateX(x, y, z);
		points[(x * sizeY + y) * sizeZ + z].stateY = getStateY(x, y, z);
	}}}
}

Lattice::~Lattice()
{
	delete[] points;
}

inline LatticePoint& Lattice::get(int x, int y, int z)
{
	return points[(x * sizeY + y) * sizeZ + z];
}

inline void Lattice::metropolisStep()
{
	static std::random_device rd;
	static std::uniform_int_distribution<int> dist{0, 100 * sizeX * sizeY * sizeZ};

	int randomNum = dist(rd);

	int alteredX = randomNum % sizeX;
	randomNum /= sizeX;
	int alteredY = randomNum % sizeY;
	randomNum /= sizeY;
	int alteredZ = randomNum % sizeZ;
	randomNum /= sizeZ;

	metropolisStepAt(alteredX, alteredY, alteredZ, randomNum);
}
	
void Lattice::metropolisStepAt(int alteredX, int alteredY, int alteredZ, int randomNum)
{
	static std::random_device rd;
	static std::mt19937 gen{rd()};
	static std::uniform_real_distribution<double> distribution(0.0, 1.0);

	int lX = (alteredX == 0)? (sizeX - 1) : (alteredX - 1);
	int uY = (alteredY == 0)? (sizeY - 1) : (alteredY - 1);
	int tZ = (alteredZ == 0)? (sizeZ - 1) : (alteredZ - 1);
	
	int rX = (alteredX + 1) % sizeX;
	int dY = (alteredY + 1) % sizeY;
	int bZ = (alteredZ + 1) % sizeZ;

	LatticePoint& alteredPoint = get(alteredX, alteredY, alteredZ);
	LatticePoint  neighbourL   = get(      lX, alteredY, alteredZ);
	LatticePoint  neighbourR   = get(      rX, alteredY, alteredZ);
	LatticePoint  neighbourU   = get(alteredX,       uY, alteredZ);
	LatticePoint  neighbourD   = get(alteredX,       dY, alteredZ);
	LatticePoint  neighbourT   = get(alteredX, alteredY,       tZ);
	LatticePoint  neighbourB   = get(alteredX, alteredY,       bZ);

	Vector spinL = stateGraph.get(neighbourL.stateX, neighbourL.stateY);
	Vector spinR = stateGraph.get(neighbourR.stateX, neighbourR.stateY);
	Vector spinU = stateGraph.get(neighbourU.stateX, neighbourU.stateY);
	Vector spinD = stateGraph.get(neighbourD.stateX, neighbourD.stateY);
	Vector spinT = stateGraph.get(neighbourT.stateX, neighbourT.stateY);
	Vector spinB = stateGraph.get(neighbourB.stateX, neighbourB.stateY);
	Vector interactionVector = (spinL + spinR + spinU + spinD + spinT + spinB) * interactivity +
	                           externalField;

	int newStateX = alteredPoint.stateX;
	int newStateY = alteredPoint.stateY;
	switch (randomNum % 4)
	{
		case 0:
		{
			newStateX = (newStateX + 1) % STATE_GRAPH_SIZE_X;
			break;
		}
		case 1:
		{
			if (newStateX == 0) newStateX = STATE_GRAPH_SIZE_X - 1;
			else                newStateX = newStateX - 1;
			break;
		}
		case 2:
		{
			newStateY = (newStateY + 1) % STATE_GRAPH_SIZE_Y;
			break;
		}
		case 3:
		{
			if (newStateY == 0) newStateY = STATE_GRAPH_SIZE_Y - 1;
			else                newStateY = newStateY - 1;
			break;
		}
	}

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
	for (int z = 0; z < sizeZ; ++z) {
		LatticePoint curPoint = get(x, y, z);
		Vector spin = stateGraph.get(curPoint.stateX, curPoint.stateY);

		magnetization += spin.z;
	}}}
	
	magnetization /= sizeX*sizeY*sizeZ;

	return magnetization;
}

double Lattice::calculateEnergy()
{
	double energy = 0.0;

	for (int x = 0; x < sizeX; ++x) {
	for (int y = 0; y < sizeY; ++y) {
	for (int z = 0; z < sizeY; ++z) {
		LatticePoint& cur        = get(              x,               y,               z);
		LatticePoint& neighbourR = get((x + 1) % sizeX,               y,               z);
		LatticePoint& neighbourD = get(              x, (y + 1) % sizeY,               z);
		LatticePoint& neighbourB = get(              x,               y, (z + 1) % sizeZ);

		Vector spin  = stateGraph.get(       cur.stateX,        cur.stateY);
		Vector spinR = stateGraph.get(neighbourR.stateX, neighbourR.stateY);
		Vector spinD = stateGraph.get(neighbourD.stateX, neighbourD.stateY);
		Vector spinB = stateGraph.get(neighbourB.stateX, neighbourB.stateY);	

		energy -= spin.scalar((spinR + spinD + spinB) * interactivity + externalField);
	}}}

	return energy;
}

#endif  // POTTS_MODEL_MODEL_HPP_INCLUDED