#ifndef GENETICS_HPP
#define GENETICS_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <random>
#include <cstdlib>
#include <ctime>
#include <algorithm>
using namespace std;

class Organism {
public:
    Organism(map<char, pair<string, string>> genes_alleles_in);

    Organism(Organism& parent1, Organism& parent2);

    Organism * reproduce();

    void die();

    void print_organism();

private:
    int age; // -1 means deceased
    double reproduction_probability;
    double death_probability;

    map<char, pair<string, string>> genes_alleles;
    pair<Organism *, Organism *> parents;
    Organism * mate;
    vector<Organism *> offspring;
};

class World {
public:
    void read_file(string filename);

    void make_world(int n);

    void time_step();

    void print_organisms();

    ~World();
private:
    size_t current_year;

    map<char, map<string, double>> genes_alleles_freqs;
    map<char, map<pair<string, string>, double>> genotypes_freqs;
    vector<Organism *> organisms;
};

#endif