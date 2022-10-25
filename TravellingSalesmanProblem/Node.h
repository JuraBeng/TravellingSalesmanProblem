#ifndef NODE_H_
#define NODE_H_
#include <iostream>
#include <vector>
class Node
{
public:
	Node(int t_idx, float t_cost, float t_matrix[5][5], std::vector<int> t_explorable, float t_bound);
	Node(int t_idx, std::vector<int> t_explorable);
	float cost;
	int idx;
	float bound;
	float matrix[5][5];
	void addChild(Node* node);
	Node* getMinChild();
	Node* solveChildren(std::vector<int> data);
	void setCost(float t_cost);
	void setMatrix(float t_matrix[5][5]);
	std::vector<int> getPath();

private:
	std::vector<int> explorable = {};
	std::vector<Node*> children = {};
};
#endif

