extern "C" {
#include "cec17.h"
}
#include <iostream>
#include <vector>
#include <random>

using namespace std;

int main() {
  vector<double> sol;
  int dim = 10;
  int seed = 42;
  std::uniform_real_distribution<> dis(-100.0, 100.0);

  for (int funcid = 1; funcid <= 30; funcid++) {
    vector<double> sol(dim);
    vector<double> bestsol(dim);
    double fitness;
    double best = -1;

    // Set the function to use in fitness
    cec17_init("random", funcid, dim);
    // If it is commented the output is print in console, instead of external files.
    // cec17_print_output();

    std::mt19937 gen(seed); // Start seed
    int evals = 0;

    while (evals < 10000*dim) {
      // Generate random solution
      for (int i = 0; i < dim; i++) {
        sol[i] = dis(gen);
      }

      // Evaluate the solution
      fitness = cec17_fitness(&sol[0]);
      // Increase count
      evals += 1;

      // Calculate the best one
      if (evals == 1 || fitness < best) {
        best = fitness;
        bestsol = sol;
      }
    }

    // Show the error of the best solution
    cout <<"Best Random[F" <<funcid <<"]: " << scientific <<cec17_error(best) <<endl;
  }

}
