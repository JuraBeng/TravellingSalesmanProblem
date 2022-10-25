#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>
#include "Node.h"
#define N 12

using namespace std;

vector<double> distances;

void print_distances(float (&distMatrix)[N][N])
{
	for (int row = 0; row < N; row++)
	{
		for (int col = 0; col < N; col++)
		{
			cout << "City: " << row << " || Second City: " << col << " || Overall distance: " << distMatrix[row][col] << endl;
		}
	}
}

float**create_distance_matrix(vector<double> xs, vector<double> ys)
{
	float **distanceMatrix = 0;

	distanceMatrix = new float* [N];
	for (int row = 0; row < N; row++)
	{
		distanceMatrix[row] = new float[N];
		for(int col = 0; col < N; col++)
		{
			if (row == col)
			{
				distanceMatrix[row][col] = 0.0f;
			}
			else
			{
				float xDist = (xs[col] - xs[row]) * (xs[col] - xs[row]);
				float yDist = (ys[col] - ys[row]) * (ys[col] - ys[row]);
				float distance = sqrtf(xDist + yDist);
				distanceMatrix[row][col] = distance;

			}
		}
	}

	return distanceMatrix;
}


void compute_distance(float distanceMatrix[N][N], int pStartingPoint, vector<double> xs, vector<double> ys)
{
	// Global variables
	vector<int> permutations;
	float minDistances[N - 1];

	for(int i = 0; i < N - 1; i++)
	{
		minDistances[i] = FLT_MAX;
	}

	for (int i = 0; i < N; i++)
	{
		if (i == pStartingPoint)
			continue;
		permutations.push_back(i);
	}
	

		float localMinDistance = FLT_MAX;
		int localStartingPoint = pStartingPoint;
#pragma omp parallel for num_threads(16)
		for (int cycle = 0; cycle < N - 1; cycle++)
		{
			int localCurrentCity = localStartingPoint;
			vector<int> localPermutation = permutations;
			std::rotate(localPermutation.begin(), localPermutation.begin() + cycle, localPermutation.begin() + cycle + 1);

			do
			{
				float totalDistance = 0.0f;

				for (int i = 0; i < localPermutation.size(); i++)
				{
					int neighborCity = localPermutation[i];
					totalDistance += distanceMatrix[localCurrentCity][neighborCity];
					localCurrentCity = neighborCity;
				}

				totalDistance += distanceMatrix[localCurrentCity][localStartingPoint];

				if (totalDistance < localMinDistance)
				{
					localMinDistance = totalDistance;
				}

			} while (next_permutation(localPermutation.begin() + 1, localPermutation.end()));
			minDistances[cycle] = localMinDistance;
		}
	float globalMinDistance = FLT_MAX;
	for(int i = 0; i < N - 1; i++)
	{
		globalMinDistance = min(globalMinDistance, minDistances[i]);
	}

	cout << "Final min distance: " << globalMinDistance << endl;
}

float reduce_matrix(float (&distanceMatrix)[N][N])
{
	return 0.0f;
}

void branch_and_bound(float distanceMatrix[N][N], int pStartingPoint)
{

}

void read_tsp_file(const char* fname)
{
	std::ifstream file(fname);
	float distanceMatrix[N][N];
	vector<double> xs, ys;

	if (file.is_open())
	{
		std::string line;

		std::getline(file, line);
		std::getline(file, line);
		std::getline(file, line);
		std::getline(file, line);
		std::getline(file, line);
		std::getline(file, line);
		std::getline(file, line);

		while (std::getline(file, line)) {
			if (line[0] == 'E')
				break;

			stringstream sin(line);
			int id;
			double x, y;
			sin >> id >> x >> y;

			xs.push_back(x);
			ys.push_back(y);
		}

		unsigned int n = xs.size();

		distances.resize(n * n);

		//TODO: calculate distance matrix
		for(int row = 0; row < N; row++)
		{
			for(int col = 0; col < N; col++)
			{
				if(row == col)
				{
					distanceMatrix[row][col] = 0.0f;
				}
				else
				{
					float xDist = (xs[col] - xs[row]) * (xs[col] - xs[row]);
					float yDist = (ys[col] - ys[row]) * (ys[col] - ys[row]);
					float distance = sqrtf(xDist + yDist);
					distanceMatrix[row][col] = distance;

				}
			}
		}

		file.close();
	}
	else
	{
		cout << fname << " file not open" << endl;
		return;
	}
	//print_distances(distanceMatrix);
	//TODO: calculate shortest possible route
	for (int i = 0; i < 10; i++)
	{
		auto start = chrono::high_resolution_clock::now();
		compute_distance(distanceMatrix, 0, xs, ys);
		auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::seconds>(end - start);
		cout << duration.count() << endl;
	}

}


float reduce_row1(float(&dist_matrix)[5][5], int row)
{
	float min_val = FLT_MAX;
	for (int i = 0; i < 5; i++)
	{
		if (dist_matrix[row][i] == FLT_MAX)
			continue;
		if (dist_matrix[row][i] <= min_val)
			min_val = dist_matrix[row][i];
	}
	if (min_val != 0.0f)
	{
		for (int i = 0; i < 5; i++)
		{
			if (dist_matrix[row][i] == FLT_MAX)
				continue;
			dist_matrix[row][i] -= min_val;
		}
	}
	return min_val;
}

float reduce_col1(float(&dist_matrix)[5][5], int col)
{
	float min_val = FLT_MAX;
	for (int i = 0; i < 5; i++)
	{
		if (dist_matrix[i][col] == FLT_MAX)
			continue;
		if (dist_matrix[i][col] <= min_val)
			min_val = dist_matrix[i][col];
	}
	if (min_val != 0.0f)
	{
		for (int i = 0; i < 5; i++)
		{
			if (dist_matrix[i][col] == FLT_MAX)
				continue;
			dist_matrix[i][col] -= min_val;
		}
	}
	return min_val;
}

float reduce_matrix1(float(&dist_matrix)[5][5])
{
	float row_sum = 0.0f;
	float col_sum = 0.0f;
	for (int i = 0; i < 5; i++)
	{
		row_sum += reduce_row1(dist_matrix, i);
	}
	for (int j = 0; j < 5; j++)
	{
		col_sum += reduce_col1(dist_matrix, j);
	}
	return col_sum + row_sum;
}

void set_mat_inf1(float(&dist_matrix)[5][5], int row, int col)
{
	dist_matrix[row][col] = FLT_MAX;
}

void set_inf_row_col1(float(&dist_matrix)[5][5], int row, int col)
{
	for (int i = 0; i < 5; i++)
	{
		set_mat_inf1(dist_matrix, row, i);
		set_mat_inf1(dist_matrix, i, col);
	}
	return;
}

void set_inf_col1(float(&dist_matrix)[5][5], int col)
{
	for (int i = 0; i < 5; i++)
	{
		set_mat_inf1(dist_matrix, i, col);
	}
}

void set_inf_row1(float(&dist_matrix)[5][5], int row)
{
	for (int i = 0; i < 5; i++)
	{
		set_mat_inf1(dist_matrix, row, i);
	}
}

void set_mat_value1(float(&dist_matrix)[5][5], int row, int col, float val)
{
	dist_matrix[row][col] = val;
}

int main()
{
	float matrix[5][5] = {
		{FLT_MAX, 20.0f, 30.0f, 10.0f, 11.0f}, //03142 10 + 6 + 4 + 7 + 3
		{15.0f, FLT_MAX, 16.0f, 4.0f, 2.0f},
		{3.0f, 5.0f, FLT_MAX, 2.0f, 4.0f},
		{19.0f, 6.0f, 18.0f, FLT_MAX, 3.0f},
		{16.0f, 4.0f, 7.0f, 16.0f, FLT_MAX}
	};
	float tmp[5][5];
	memcpy(tmp, matrix, 25 * sizeof(float));
	float bound = reduce_matrix1(matrix);
	Node* root = new Node(0, 0.0f, matrix, std::vector<int>{}, bound);
	Node* curr_node = root;
	std::vector<int> vc = { 1,2,3, 4 };

	while (true)
	{
		Node *node = curr_node->solveChildren(vc);
		vc.erase(std::remove(vc.begin(), vc.end(), node->idx), vc.end());
		curr_node = node;
		if(vc.empty())
			break;
	}
	std::vector<int> path = root->getPath();
	path.push_back(0);
	float cost = 0.0f;
	for(int i = 0; i < path.size() - 1; i++)
	{
		cost += tmp[path[i]][path[i + 1]];
	}
	std::cout << cost << std::endl;
	return 0;
}
