#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

#include <cstdlib>
#include <vector>

// To represent a data point corresponding to x and y = f(x)
class Node {
    private:
        double x_;
        double y_;
    public:
        // Constructor
        Node(double x, double y) : x_(x), y_(y) {};

        double x() const { return x_; }
        double y() const { return y_; }
};

class Nodes {
    private:
        std::vector<Node> nodes_;

    public:
        // Constructor
        Nodes(std::vector<Node> nodes) : nodes_(nodes) {}

        // Get the number of nodes
        int size() const { return nodes_.size(); }
        // Get the node at index i
        Node get(int i) const { return nodes_[i]; }
        // Operator to get node at index with braces
        Node operator[](int i) const { return nodes_[i]; }
};

// function to interpolate the given data points using Lagrange's formula
// xi corresponds to the new data point whose value is to be obtained
double *lagrangeInterpolation(Nodes, int, double[]);

double *chebyshevInterpolation(double, double, int, int, double[], double (double));

#endif