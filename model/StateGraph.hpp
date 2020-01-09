// No Copyright. Vladislav Aleinik 2019
#ifndef POTTS_MODEL_STATE_GRAPH_HPP_INCLUDED
#define POTTS_MODEL_STATE_GRAPH_HPP_INCLUDED

#include "Vector.hpp"
#include <cmath>
#include <cstdio>
#include <assert.h>

#define STATE_GRAPH_SIZE_X 1
#define STATE_GRAPH_SIZE_Y 2

struct FibonacciSpinStateGraph
{
	Vector spins[STATE_GRAPH_SIZE_X * STATE_GRAPH_SIZE_Y];

	FibonacciSpinStateGraph();
	inline Vector get(int x, int y) const;
};

FibonacciSpinStateGraph::FibonacciSpinStateGraph() : 
	spins ({})
{
	for (size_t x = 0; x < STATE_GRAPH_SIZE_X; ++x)
	{
		for (size_t y = 0; y < STATE_GRAPH_SIZE_Y; ++y)
		{
			double longtitude = 2.0 * M_PI * x/(STATE_GRAPH_SIZE_X);
			double latitude = acos(2.0 * y/(STATE_GRAPH_SIZE_Y-1) - 1) - M_PI/2;

			spins[x * STATE_GRAPH_SIZE_Y + y].x = cos(latitude) * cos(longtitude);
			spins[x * STATE_GRAPH_SIZE_Y + y].y = cos(latitude) * sin(longtitude);
			spins[x * STATE_GRAPH_SIZE_Y + y].z = sin(latitude);

			// printf("[%ld, %ld] {%f, %f, %f}\n", x, y, 
			//         spins[x * STATE_GRAPH_SIZE_Y + y].x,
			//         spins[x * STATE_GRAPH_SIZE_Y + y].y,
			//         spins[x * STATE_GRAPH_SIZE_Y + y].z);
		}
	}
}

inline Vector FibonacciSpinStateGraph::get(int x, int y) const
{	
	return spins[x * STATE_GRAPH_SIZE_Y + y];
}

#endif  // POTTS_MODEL_STATE_GRAPH_HPP_INCLUDED