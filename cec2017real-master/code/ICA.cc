extern "C" {
#include "cec17.h"
}
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

using namespace std;

enum CountryType { COUNTRY, IMPERIALIST, COLONY };

class Country {
protected:
    std::vector<double> solution;
    double fitness;
    CountryType type;

public:
    Country(const std::vector<double>& solution, double fitness, CountryType type)
        : solution(solution), fitness(fitness), type(type) {}

    double getFitness() const {
        return fitness;
    }

    void setFitness(double newFitness) {
        fitness = newFitness;
    }

    std::vector<double> getSolution() const {
        return solution;
    }

    void setSolution(const std::vector<double>& newSolution) {
        solution = newSolution;
    }

    CountryType getType() const {
        return type;
    }

    void setType(CountryType newType) {
        type = newType;
    }
};

class Imperialist : public Country {
private:
    std::vector<Country*> colonies; // List of colonies

public:
    Imperialist(const std::vector<double>& solution, double fitness)
        : Country(solution, fitness, IMPERIALIST) {}

    // Method to add a colony to the list
    void addColony(Country* colony) {
        colonies.push_back(colony);
    }

    void removeColony(Country* colony) {
        colonies.erase(std::remove(colonies.begin(), colonies.end(), colony), colonies.end());
    }

    // Method to get the list of colonies
    std::vector<Country*> getColonies() const {
        return colonies;
    }
};

class Colony : public Country {
public:
    Colony(const std::vector<double>& solution, double fitness)
        : Country(solution, fitness, COLONY) {}
};

class ICA{
private:
    vector<Country> CountryList;
    vector<Imperialist*> Imperialists;
    int dim;
    int max_evals;
    double ImperialistPercentage;
    int numEmpires;
    int numColonies;
    double phi;
    double beta;
    double xi;
    double revolution_rate;
    int evals;
    int seed;
    mt19937 gen;

public:
    ICA(int dim, int max_evals, double ImperialistPercentage, vector<Country> CountryList, double phi, double beta, double xi, double revolution_rate, int seed){
        this->dim = dim;
        this->max_evals = max_evals;
        this->ImperialistPercentage = ImperialistPercentage;
        this->numEmpires = round(ImperialistPercentage * CountryList.size());
        this->numColonies = CountryList.size() - numEmpires;
        this->phi = phi;
        this->beta = beta;
        this->xi = xi;
        this->revolution_rate = revolution_rate;
        this->evals = CountryList.size();
        this->seed = seed;
        this->CountryList = CountryList;
        initializeEmpires();
        gen = mt19937(seed);
    }

    vector<double> ImperialistPower() {
        vector<double> costs;
        for (auto* imp : Imperialists) {
            costs.push_back(imp->getFitness());
            //cout << imp->getFitness() << endl;
        }

        double max_cost = *max_element(costs.begin(), costs.end());
        vector<double> normalized_costs;
        for (auto cost : costs) {
            normalized_costs.push_back(max_cost - cost);
        }

        double total_normalized_cost = accumulate(normalized_costs.begin(), normalized_costs.end(), 0.0);
        vector<double> powers;
        for (auto normalized_cost : normalized_costs) {
            powers.push_back(normalized_cost / total_normalized_cost);
        }

        return powers;
    }

    void initializeEmpires(){
    vector<Country> countries = CountryList;
    // Ordena los países por aptitud en orden ascendente (menor fitness primero)
    sort(countries.begin(), countries.end(), [](Country a, Country b) { return a.getFitness() < b.getFitness(); });

    // Inicializa los imperialistas con los países de menor fitness
    for (int i = 0; i < numEmpires; i++){
        Imperialist *imp = new Imperialist(countries.front().getSolution(), countries.front().getFitness());
        Imperialists.push_back(imp);
        countries.erase(countries.begin());
    }

    // Calcula el poder de los imperialistas
    vector<double> power = ImperialistPower();
    vector<int> num_colonies_imperialist;
    for (int i = 0; i < numEmpires; i++){
        num_colonies_imperialist.push_back(round(power[i] * numColonies));
    }

    // Distribuye las colonias
    for (int i = 0; i < numEmpires; i++){
        for (int j = 0; j < num_colonies_imperialist[i]; j++){
            if (!countries.empty()){
                int index = rand() % countries.size();
                Colony *col = new Colony(countries[index].getSolution(), countries[index].getFitness());
                Imperialists[i]->addColony(col);
                countries.erase(countries.begin() + index);
            }
        }
    }

    // Maneja imperios sin colonias
    for (int i = 0; i < Imperialists.size(); i++){
        if (Imperialists[i]->getColonies().size() == 0){
            Colony *col = new Colony(Imperialists[i]->getSolution(), Imperialists[i]->getFitness());
            Imperialists.erase(Imperialists.begin() + i);
            countries.push_back(Country(col->getSolution(), col->getFitness(), COUNTRY));
            i--; // Ajusta el índice ya que hemos eliminado un elemento
        }
    }

    // Distribuye las colonias restantes
    while(!countries.empty()){
        int index = rand() % Imperialists.size();
        Colony *col = new Colony(countries.back().getSolution(), countries.back().getFitness());
        Imperialists[index]->addColony(col);
        countries.pop_back();
    }
}



    double EuclideanDistance(vector<double> a, vector<double> b){
        double sum = 0;
        for (int i = 0; i < a.size(); i++){
            sum += pow(a[i] - b[i], 2);
        }
        return sqrt(sum);
    }

    double EuclideanDistance(Country a, Country b){
        return EuclideanDistance(a.getSolution(), b.getSolution());
    }


    void advance(vector<double>& a, vector<double> b, double step){
        for (int i = 0; i < a.size(); i++){
            if (a[i] < b[i]){
                a[i] = a[i] + step;
            } else {
                a[i] = a[i] - step;
            }
        }
    }

    void deviation(vector<double>& a, double theta, double step){
        double p = abs(theta)/2*3.1416;
        int changes = round(p * a.size());
        for (int i=0; i<changes; i++){
            int index = rand() % a.size();
            a[index] = a[index] + step;
        }
    }
    

    void Assimilation(){
        uniform_real_distribution<> dis2(-phi, phi);
        for (int i = 0; i < Imperialists.size(); i++){
            Imperialist *imp = Imperialists[i];
            vector<Country*> colonies = imp->getColonies();
            for (int j = 0; j < colonies.size(); j++){
                if (evals >= max_evals){
                    break;
                }
                Country *col = colonies[j];
                double d = EuclideanDistance(*imp, *col);
                uniform_real_distribution<> dis1(0, beta*d);
                vector<double> new_sol = col->getSolution();
                double x = dis1(gen);
                double theta = dis2(gen);
                int steps = d/10;
                int step = 0;
                while (step < steps){
                    advance(new_sol, imp->getSolution(), 10);
                    deviation(new_sol, theta, 1);
                    step++;
                }
                double new_fitness = cec17_fitness(&new_sol[0]);
                evals++;
                col->setSolution(new_sol);
                col->setFitness(new_fitness);
            }
        }
    }

    void Revolution() {
    std::uniform_int_distribution<> dis_index(0, dim - 1); // Distribución para seleccionar índices
    std::normal_distribution<> dis_mutation(0, 1); // Distribución normal para las mutaciones
    std::uniform_real_distribution<> dis_probability(0.0, 1.0); // Distribución para la probabilidad de revolución

    for (auto& imp : Imperialists) {
        auto colonies = imp->getColonies();
        for (auto& col : colonies) {
            if (dis_probability(gen) < revolution_rate) { // Condición para aplicar la revolución
                if (evals >= max_evals) {
                    break; // Salir si se alcanza el máximo de evaluaciones
                }

                auto new_sol = col->getSolution(); // Obtener la solución actual
                int changes = std::max(1, static_cast<int>(dim * 0.05)); // Número de cambios basado en el 5% de las dimensiones

                for (int i = 0; i < changes; i++) {
                    int index = dis_index(gen); // Seleccionar índice aleatorio
                    new_sol[index] += dis_mutation(gen); // Aplicar mutación en el índice seleccionado
                }

                double new_fitness = cec17_fitness(&new_sol[0]); // Evaluar la nueva solución
                evals++;
                col->setSolution(new_sol); // Actualizar la solución de la colonia
                col->setFitness(new_fitness); // Actualizar el fitness de la colonia
            }
        }
    }
}


    void Relocation() {
        for (int i = 0; i < Imperialists.size(); i++) {
            Imperialist *imp = Imperialists[i];
            vector<Country*> colonies = imp->getColonies();
            for (int j = 0; j < colonies.size(); j++) {
                Country* col = colonies[j];

                if (col->getFitness() < imp->getFitness()) {
                    // Swap the solutions and fitness values of the imperialist and the colony
                    vector<double> tempSolution = imp->getSolution();
                    double tempFitness = imp->getFitness();

                    imp->setSolution(col->getSolution());
                    imp->setFitness(col->getFitness());

                    col->setSolution(tempSolution);
                    col->setFitness(tempFitness);
                }
            }
        }
    }

    vector<double> EmpiresTotalPower() {
        vector<double> powers;
        for (int i = 0; i < Imperialists.size(); i++) {
            Imperialist* imp = Imperialists[i];
            double imp_power = imp->getFitness();
            vector<Country*> colonies = imp->getColonies();
            double cols_power = 0;
            for (int j = 0; j < colonies.size(); j++) {
                cols_power += colonies[j]->getFitness();
            }
            powers.push_back(imp_power + xi * cols_power);
        }

        double max_power = *min_element(powers.begin(), powers.end());
        for (int i = 0; i < powers.size(); i++) {
            powers[i] -= max_power;
        }
        for (int i = 0; i < powers.size(); i++) {
            powers[i] = powers[i] / accumulate(powers.begin(), powers.end(), 0.0);
        }
        return powers;
    }

    void WeakestColonyTransfer() {
    Country* weakestColony = nullptr;
    Imperialist* weakestImperialist = nullptr;
    double maxFitnessImp = numeric_limits<double>::lowest();
    double maxFitnessCol = numeric_limits<double>::lowest();

    for (const auto& imp : Imperialists) {
        if (imp->getFitness() > maxFitnessImp) {
            maxFitnessImp = imp->getFitness();
            weakestImperialist = imp;
        }
    }

    for (const auto& col : weakestImperialist->getColonies()) {
        if (col->getFitness() > maxFitnessCol) {
            maxFitnessCol = col->getFitness();
            weakestColony = col;
        }
    }

    weakestImperialist->removeColony(weakestColony);
    vector<double> p = EmpiresTotalPower();
    vector<double> r;
    uniform_real_distribution<> dis(0, 1);

    for (size_t i = 0; i < Imperialists.size(); ++i) {
        r.push_back(dis(gen));
    }

    vector<double> d;
    for (size_t i = 0; i < Imperialists.size(); ++i) {
        d.push_back(p[i] - r[i]);
    }

    int index = distance(d.begin(), max_element(d.begin(), d.end()));
    Imperialists[index]->addColony(weakestColony);
    }


    void EmpireElimination() {
        vector<Colony*> new_colonies;

    for (auto it = Imperialists.begin(); it != Imperialists.end(); ) {
        if ((*it)->getColonies().empty()) {
            Colony *col = new Colony((*it)->getSolution(), (*it)->getFitness());
            new_colonies.push_back(col);
            it = Imperialists.erase(it); // Eliminar y avanzar el iterador
        } else {
            ++it;
        }
    }

    for (auto& new_col : new_colonies) {
        Imperialist* imp = Imperialists[rand() % Imperialists.size()];
        imp->addColony(new_col);
    }
    }


    void Competition(vector<double> &sol, double &fitness){
        while (evals<max_evals){
            Assimilation();
            Revolution();
            Relocation();
            WeakestColonyTransfer();
            EmpireElimination();
            sort(Imperialists.begin(), Imperialists.end(), [](Imperialist* a, Imperialist* b) { return a->getFitness() > b->getFitness(); });
        }
        //cout << "Evals: " << evals << endl;
        sort(Imperialists.begin(), Imperialists.end(), [](Imperialist* a, Imperialist* b) { return a->getFitness() > b->getFitness(); });
        sol = Imperialists[0]->getSolution();
        fitness = Imperialists[0]->getFitness();
        //cout << "Evals: " << evals << " Num empires " << Imperialists.size() << endl;
    }
    
};


int main() {
  vector<double> sol;
  int dim = 10;
  int seed = 42;
  std::uniform_real_distribution<> dis(-100.0, 100.0);
  
  for (int funcid = 1; funcid <= 30; funcid++) {
    vector<double> sol(dim);
    vector<double> bestsol(dim);
    const int num_countries = dim*20;
    double fitness, bestfitness = -1;

    cec17_init("ICA", funcid, dim);

    cerr <<"Warning: output by console, if you want to create the output file you have to comment cec17_print_output()" <<endl;
    //cec17_print_output(); // Comment to generate the output file

    std::mt19937 gen(seed); // Inicio semilla
    vector<Country> CountryList;
    for (int c = 0; c < num_countries; c++) {
      for (int i = 0; i < dim; i++) {
        sol[i] = dis(gen);
      }
      Country country = Country(sol, cec17_fitness(&sol[0]), COUNTRY);
      CountryList.push_back(country);
    }
    ICA ica(dim, 10000*dim, 0.3, CountryList, 3.1416/3, 10, 0.3, 0.0, seed);
    ica.Competition(sol, fitness);
    bestfitness = fitness;
    bestsol = sol;


    cout <<"Best Random[F" <<funcid <<"]: " << scientific <<cec17_error(bestfitness) <<endl;
  }

  return 0;
}
