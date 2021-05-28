//	-------------------------------------------------------------------------
//	Smallest circumscribing k-gon. https://github.com/atyuwen/smallest-k-gon.
//	Contact : atyuwen@gmail.com
//	Author : Yuwen Wu (https://atyuwen.github.io/)
//	License : Distributed under the MIT License.
//	-------------------------------------------------------------------------

#include <cassert>
#include <random>
#include <iostream>

#include "../include/convex.h"

static float random01()
{
	return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

int main()
{
	// Generate M points randomly.
	const int M = 20;
	std::vector<ayw::float2> points(M);
	for (int i = 0; i < M; ++i)
	{
		points.emplace_back(random01(), random01());
	}

	// Build convex hull from those points.
	ayw::convex n_gon;
	n_gon.build(points.begin(), points.end());
	std::cout << "n-gon. n = " << n_gon.vertices.size() << ", area = " << n_gon.area();

	for (int k = 8; k >=4; --k)
	{
		// Simplify the n-gon to k-gon.
		ayw::convex k_gon = n_gon;
		k_gon.simplify(k);

		// Check whether the simplified k-gon contains all the points.
		for (int i = 0; i < M; ++i)
		{
			if (!k_gon.contains(points[i]))
			{
				// Should not happens.
				std::cout << "error occurred.";
			}
		}

		std::cout << "k-gon. k = " << k << ", area = " << k_gon.area();

		// Clip the k-gon to [0, 1] rect.
		k_gon.clip();
		std::cout << "k-gon clipped. v = " << k_gon.vertices.size() << ", area = " << k_gon.area();
	}
}
