#include "genetics.hpp"

Organism::Organism(map<char, pair<string, string>> genes_alleles_in) 
 : age(0), reproduction_probability(0), death_probability(0), 
   genes_alleles(genes_alleles_in), parents(nullptr, nullptr), 
   mate(nullptr), offspring({ }) {}

Organism::Organism(Organism& parent1, Organism& parent2)
: age(0), reproduction_probability(0), death_probability(0),
  parents(&parent1, &parent2), mate(nullptr), offspring({ }) {
    vector<pair<string, string>> genotype_options = { };
    auto it2 = (parent2.genes_alleles).begin();

    for (
        auto it1 = (parent1.genes_alleles).begin(); 
        it1 != (parent1.genes_alleles).end(); 
        it1++
    ) {
        genotype_options.push_back({it1->second.first, it2->second.first});
        genotype_options.push_back({it1->second.first, it2->second.second});
        genotype_options.push_back({it1->second.second, it2->second.first});
        genotype_options.push_back({it1->second.second, it2->second.second});

        srand(time(0));
        int rand_num = rand() % 3;

        genes_alleles.insert({it1->first, genotype_options.at(rand_num)});

        it2++;
    }

    parent1.offspring.push_back(this);
    parent2.offspring.push_back(this);
}

Organism * Organism::reproduce() {
    return new Organism(*this, *mate);
}

void World::read_file(string filename) {
    ifstream reader(filename);

    if (!reader.is_open()) {
        return;
    }

    string line;
    // Read one entire line from the file at a time
    while (getline(reader, line)) {
        
        // Skip empty lines just in case
        if (line.empty()) {
            continue;
        }

        // Turn the line into a stream so we can parse it just like 'reader'
        istringstream iss(line);
        
        char gene;
        // Try to read the first character as the gene name
        if (iss >> gene) {
            
            string allele;
            double freq;
            
            // Loop through the rest of the line, reading pairs until the line ends
            while (iss >> allele >> freq) {
                genes_alleles_freqs[gene][allele] = freq;
            }
        }
    }
}

void World::make_world(int n) {
    if(genes_alleles_freqs.empty()) {
        return;
    }

    current_year = 0;
    
    for (const auto& gene_pair : genes_alleles_freqs) {
        char gene = gene_pair.first;
        const auto& alleles = gene_pair.second;

        // Outer loop: pick the first allele
        for (auto it1 = alleles.begin(); it1 != alleles.end(); ++it1) {
            
            // Inner loop: pick the second allele 
            // Note: it2 starts at it1 so we don't repeat combinations like (Allele2, Allele1)
            for (auto it2 = it1; it2 != alleles.end(); ++it2) {
                
                string a1_name = it1->first;
                double a1_freq = it1->second;
                
                string a2_name = it2->first;
                double a2_freq = it2->second;

                // Calculate the frequency based on Hardy-Weinberg
                if (it1 == it2) {
                    // Homozygous (same allele): p^2
                    genotypes_freqs[gene][{a1_name, a2_name}] = a1_freq * a1_freq;
                } else {
                    // Heterozygous (different alleles): 2pq
                    genotypes_freqs[gene][{a1_name, a2_name}] = 2.0 * a1_freq * a2_freq;
                }
            }
        }
    }
    
    using MultiGeneGenotype = map<char, pair<string, string>>;
    
    // A list to hold all final combinations and their total calculated frequencies
    // Each item is a pair: { The Combination, The Total Frequency }
    vector<pair<MultiGeneGenotype, double>> combined_freqs;

    // Initialize with a single "empty" combination and a starting probability of 1.0. 
    // We need this base state to start multiplying against!
    combined_freqs.push_back({MultiGeneGenotype(), 1.0});

    // --- 2. Build the Combinations Dynamically ---

    // Go through one gene at a time
    for (const auto& gene_pair : genotypes_freqs) {
        char current_gene = gene_pair.first;
        const auto& single_gene_genotypes = gene_pair.second;

        // Create a temporary list to hold the new, expanded combinations for this step
        vector<pair<MultiGeneGenotype, double>> next_combined_freqs;

        // Take every combination we've built so far...
        for (const auto& existing_combo : combined_freqs) {
            
            // ...and combine it with every possible genotype of the CURRENT gene
            for (const auto& geno_pair : single_gene_genotypes) {
                
                // Make a copy of the existing combination map
                MultiGeneGenotype new_combo_map = existing_combo.first;
                
                // Add the current gene's genotype to it
                new_combo_map[current_gene] = geno_pair.first; 

                // Multiply the running frequency by the new genotype's frequency
                double new_combo_freq = existing_combo.second * geno_pair.second;

                // Save this new, expanded combination to our temporary list
                next_combined_freqs.push_back({new_combo_map, new_combo_freq});
            }
        }
        
        // Overwrite our main list with the newly expanded list before moving to the next gene
        combined_freqs = next_combined_freqs;
    }


    struct Allocation {
        MultiGeneGenotype genotype;
        int count;
        double remainder;
    };

    vector<Allocation> allocations;
    int current_total = 0;

    // 1. Calculate base counts and remainders
    for (const auto& combo : combined_freqs) {
        double exact_expected = combo.second * n;
        int base_count = static_cast<int>(exact_expected); // Rounds down (e.g., 33.8 becomes 33)
        double decimal_remainder = exact_expected - base_count; // The leftover (e.g., 0.8)

        allocations.push_back({combo.first, base_count, decimal_remainder});
        current_total += base_count;
    }

    // 2. Figure out how many organisms we are short due to rounding down
    int missing_organisms = n - current_total;

    // 3. Sort the allocations by their remainder, highest to lowest
    // (A lambda function tells std::sort how to compare our custom structs)
    sort(allocations.begin(), allocations.end(), [](const Allocation& a, const Allocation& b) {
        return a.remainder > b.remainder;
    });

    // 4. Hand out the missing organisms to those with the highest remainders
    for (int i = 0; i < missing_organisms; ++i) {
        allocations[i].count++; 
    }

    // 5. Actually build the Organism objects using the finalized exact counts
    organisms.reserve(n);

    for (const auto& alloc : allocations) {
        for (int i = 0; i < alloc.count; ++i) {
            // Using the requested constructor
            organisms.push_back(new Organism(alloc.genotype));
        }
    }
}

void Organism::print_organism() {
    cout << "  Age: " << age << endl;
    cout << "  Reproduction Probability: " << reproduction_probability << endl;
    cout << "  Death Probability: " << death_probability << endl;
    cout << "  Genes:" << endl;
    for (auto it = genes_alleles.begin(); it != genes_alleles.end(); it++) {
        cout << "    " << it->first << ": " 
            << it->second.first << it->second.second << endl;
    }
    cout << "  # of Offspring: " << offspring.size() << endl;
}

void World::print_organisms() {
    for (int i = 0; i < organisms.size(); i++) {
        cout << "Organism " << i + 1 << ":" << endl;
        organisms[i]->print_organism();
        cout << endl;
    }
}

World::~World() {
    for (int i = 0; i < organisms.size(); i++) {
        delete organisms[i];
    }
}

int main() {
    World newWorld;

    newWorld.read_file("input01.txt");
    newWorld.make_world(25);
    newWorld.print_organisms();

    return 0;
}