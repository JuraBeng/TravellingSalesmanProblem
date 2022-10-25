#include "Node.h"


	float reduce_row(float(&dist_matrix)[5][5], int row)
	{
		float min_val = FLT_MAX;
		for (int i = 0; i < 5; i++)
		{
			if (dist_matrix[row][i] == FLT_MAX)
				continue;
			if (dist_matrix[row][i] <= min_val)
				min_val = dist_matrix[row][i];
		}
		if (min_val != 0.0f && min_val != FLT_MAX)
		{
			for (int i = 0; i < 5; i++)
			{
				if (dist_matrix[row][i] == FLT_MAX)
					continue;
				dist_matrix[row][i] -= min_val;
			}

		}
		if (min_val == FLT_MAX)
			min_val = 0.0f;
		return min_val;
	}

	float reduce_col(float(&dist_matrix)[5][5], int col)
	{
		float min_val = FLT_MAX;
		for (int i = 0; i < 5; i++)
		{
			if (dist_matrix[i][col] == FLT_MAX)
				continue;
			if (dist_matrix[i][col] <= min_val)
				min_val = dist_matrix[i][col];
		}
		if (min_val != 0.0f && min_val != FLT_MAX)
		{
			for (int i = 0; i < 5; i++)
			{
				if (dist_matrix[i][col] == FLT_MAX)
					continue;
				dist_matrix[i][col] -= min_val;
			}
		}
		if (min_val == FLT_MAX)
			min_val = 0.0f;
		return min_val;
	}

	float reduce_matrix(float(&dist_matrix)[5][5])
	{
		float row_sum = 0.0f;
		float col_sum = 0.0f;
		for (int i = 0; i < 5; i++)
		{
			row_sum += reduce_row(dist_matrix, i);
		}
		for (int j = 0; j < 5; j++)
		{
			col_sum += reduce_col(dist_matrix, j);
		}
		return col_sum + row_sum;
	}

	 void set_mat_inf(float(&dist_matrix)[5][5], int row, int col)
	{
		dist_matrix[row][col] = FLT_MAX;
	}

	 void set_inf_row_col(float(&dist_matrix)[5][5], int row, int col)
	{
		for (int i = 0; i < 5; i++)
		{
			set_mat_inf(dist_matrix, row, i);
			set_mat_inf(dist_matrix, i, col);
		}
		return;
	}

	 void set_inf_col(float(&dist_matrix)[5][5], int col)
	{
		for (int i = 0; i < 5; i++)
		{
			set_mat_inf(dist_matrix, i, col);
		}
	}

	 void set_inf_row(float(&dist_matrix)[5][5], int row)
	{
		for (int i = 0; i < 5; i++)
		{
			set_mat_inf(dist_matrix, row, i);
		}
	}

	 void set_mat_value(float(&dist_matrix)[5][5], int row, int col, float val)
	{
		dist_matrix[row][col] = val;
	}



Node::Node(int t_idx, float t_cost, float t_matrix[5][5], std::vector<int> t_explorable, float t_bound)
{
	bound = t_bound;
	explorable = t_explorable;
	idx = t_idx;
	cost = t_cost;
	memcpy(matrix, t_matrix, 25 * sizeof(float));
}

Node::Node(int t_idx, std::vector<int> t_explorable)
{
	explorable = t_explorable;
	idx = t_idx;
}

void Node::addChild(Node* node)
{
	children.push_back(node);
}

Node* Node::getMinChild()
{
	float min = FLT_MAX;
	int idx = 0;
	if (children.empty())
		return nullptr;
	for (int i = 0; i < children.size(); i++)
	{
		if(children[i]->cost <= min)
		{
			min = children[i]->cost;
			idx = i;
		}
	}
	return children[idx];
}

Node* Node::solveChildren(std::vector<int> data)
{
	std::vector<int> explore = explorable;
	explore.push_back(idx);
	for(int i = 0; i < data.size(); i++)
	{
		Node* node = new Node(data[i], explore);

		float tmp_mat[5][5];
		memcpy(tmp_mat, matrix, 25 * sizeof(float));

		//set rows of all predecessors to inf
		for(int i = 0; i < explore.size(); i++)
		{
			set_inf_row(tmp_mat, explore[i]);
		}
		set_inf_col(tmp_mat, node->idx);
		// set mat[curr][root] to inf
		float cache = tmp_mat[node->idx][explore[0]];
		set_mat_value(tmp_mat, node->idx, explore[0], FLT_MAX);
		float reduction = reduce_matrix(tmp_mat);

		// compute cost and set all other things to node
		float cost = reduction + matrix[idx][node->idx] + bound;

		node->setCost(cost);
		node->bound = cost;
		tmp_mat[node->idx][explore[0]] = cache;
		node->setMatrix(tmp_mat);
		children.push_back(node);
	}

	return getMinChild();
}

void Node::setCost(float t_cost)
{
	cost = t_cost;
}

void Node::setMatrix(float t_matrix[5][5])
{
	memcpy(matrix, t_matrix, 25 * sizeof(float));
}

std::vector<int> Node::getPath()
{
	std::vector<int> path_order;
	path_order.push_back(idx);
	Node* curr = getMinChild();
	while (curr != NULL)
	{
		path_order.push_back(curr->idx);
		curr = curr->getMinChild();
	}
	return path_order;
}

