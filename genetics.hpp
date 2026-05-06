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
#include<bits/stdc++.h>
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

    void find_mate_for(vector<Organism *> * organisms_in);

    void set_mate(Organism * mate_in);

    int get_age();

    string get_mating_status();

    map<char, pair<string, string>> * get_genes_alleles();

    // Returns a pointer to a newly created organism OR nullptr
    Organism * organism_time_step(vector<Organism *> * organisms_in);

private:
    int age; // -1 means deceased
    double mating_probability;
    double reproduction_probability;
    double death_probability;
    string mating_status; // "LOOKING", "NOT LOOKING", "HAS MATE"

    map<char, pair<string, string>> genes_alleles;
    pair<Organism *, Organism *> parents;
    Organism * mate;
    vector<Organism *> offspring;
};

class World {
public:
    void read_file(string filename);

    void make_world();

    void time_step();

    void print_organisms();

    void print_allele_frequencies();

    void add_new_organism(Organism * organism_in);

    Organism * get_organism(int i); // Indexing starting from i = 1.

    size_t get_current_year();

    size_t get_year_limit();

    ~World();
private:
    size_t current_year;
    size_t year_limit;
    size_t initial_organism_count;

    map<char, map<string, double>> genes_alleles_freqs;
    map<char, map<pair<string, string>, double>> genotypes_freqs;
    vector<Organism *> organisms;
};

#endif