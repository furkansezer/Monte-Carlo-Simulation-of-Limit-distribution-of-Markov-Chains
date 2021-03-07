#include <iostream>
#include <random>
#include <cmath>
#include <chrono>
#include <fstream>
using namespace std;
#define EPSI 0.0005


class MATRIX {
public:
	MATRIX(int s) {
		int i, j, k;
		double row_sum, control, random_prob;
		size = s;

		int m, random_state;
		steady_dist = new float[size];
		for (int i = 0; i < size; i++) {
			steady_dist[i] = 0;
		}

		t = new float* [size];
		for (i = 0; i < size; i++) {
			t[i] = new float[size];
		}


		for (i = 0; i < size; i++) {
			row_sum = 1.0;

			for (j = 0; j < size; j++) {

				unsigned seed = chrono::system_clock::now().time_since_epoch().count();
				default_random_engine generator(seed);
				uniform_real_distribution<float> prob(0.01, 0.02);
				prob.reset();
				control = prob(generator);

				t[i][j] = control;
				row_sum -= control;

			}
			for (k = 0; k < size; k++) {

				unsigned seed_new = chrono::system_clock::now().time_since_epoch().count();
				default_random_engine generator_2(seed_new);

				uniform_real_distribution<float> prob_2(0, row_sum);

				random_state = rand() % size;
				random_prob = prob_2(generator_2);
				t[i][random_state] += random_prob;
				row_sum -= random_prob;

			}
			m = rand() % size;
			t[i][m] += row_sum;

		}
	}

	void monte_carlo() {
		int current_state = rand() % size;
		steady_dist[current_state]++;
		int i, j;
		float current_random;
		int minimum_of_greater_states, maximum_of_smaller_states;

		for (i = 0; i < 200000; i++) {

			unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed);
			uniform_real_distribution<float> our_r(0.0, 1.0);
			current_random = our_r(generator);
			minimum_of_greater_states = -1;
			maximum_of_smaller_states = -1;
			for (j = 0; j < size; j++) {
				if (isgreater(current_random, t[current_state][j])) { //isgreater is a macro in cmath library.
					if (maximum_of_smaller_states == -1) {
						maximum_of_smaller_states = j;
					}
					else {
						if (isgreater(t[current_state][j], t[current_state][maximum_of_smaller_states])) {
							maximum_of_smaller_states = j;
						}
					}
				}
				else {
					if (minimum_of_greater_states == -1) {
						minimum_of_greater_states = j;
					}
					else {
						if (isgreater(t[current_state][minimum_of_greater_states], t[current_state][j])) {
							minimum_of_greater_states = j;
						}
					}

				}
			}
			if (minimum_of_greater_states == -1) {
				current_state = maximum_of_smaller_states;
				steady_dist[maximum_of_smaller_states]++;
			}
			else {
				current_state = minimum_of_greater_states;
				steady_dist[minimum_of_greater_states]++;
			}


		}
		for (i = 0; i < size; i++) {
			steady_dist[i] /= (float)200000;
		}

	}
	void matrix_multiplication() {
		float** current_matrix, ** multiplication;
		float* column_average = new float[size];
		current_matrix = new float* [size];
		multiplication = new float* [size];
		double current_column_sum, total_squared_distance;
		int i, j, k;
		bool convergence_is_reached = false;
		for (i = 0; i < size; i++) {
			current_matrix[i] = new float[size];
			multiplication[i] = new float[size];
		}
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				current_matrix[i][j] = t[i][j];
			}
		}
		while (convergence_is_reached == false) {

			convergence_is_reached = true;
			total_squared_distance = 0.0;

			for (i = 0; i < size; i++) {
				for (j = 0; j < size; j++) {
					multiplication[i][j] = 0.0;
				}
			}

			for (i = 0; i < size; i++) {
				for (j = 0; j < size; j++) {
					for (k = 0; k < size; k++)
					{

						multiplication[i][j] += current_matrix[i][k] * current_matrix[k][j];
					}
				}
			}


			for (i = 0; i < size; i++) {

				for (j = 0; j < size; j++) {

					current_matrix[i][j] = multiplication[i][j];
				}
			}

			for (i = 0; i < size; i++) {
				current_column_sum = 0;
				for (j = 0; j < size; j++) {
					current_column_sum += current_matrix[j][i];
				}
				column_average[i] = current_column_sum / (double)size;
			}
			int r = rand() % size;
			for (i = 0; i < size; i++) {
				total_squared_distance = pow(column_average[i] - current_matrix[r][i], 2.0);
			}
			if (isgreater(sqrt(total_squared_distance), EPSI)) {
				convergence_is_reached = false;

			}

		}



		for (i = 0; i < size; i++) {

			steady_dist[i] = column_average[i];
		}
		for (int i = 0; i < size; i++) {
			delete[] * (current_matrix + i);
			delete[] * (multiplication + i);
		}
		delete[] current_matrix;
		delete[] multiplication;
		delete[] column_average;


	}
	void print_steady(ofstream& ofs) {

		for (int i = 0; i < size; i++) {
			ofs << steady_dist[i] << " ";
		}
		ofs << endl;
	}
	void print(ofstream& ofs) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				ofs << t[i][j] << " ";
			}
			ofs << endl;
		}
	}
	void absorbtion() {
		t[0][0] = 1;
		for (int i = 1; i < size; i++) {
			t[0][i] = 0;
		}
	}
	~MATRIX() {
		for (int i = 0; i < size; i++) {
			delete[] * (t + i);
		}
		delete[] t;
		delete[] steady_dist;

	}


private:
	float** t;
	int size;
	float* steady_dist;


};

int main() {
	ofstream fs1("C:\\Users\\Desktop\\nfive.txt", ios::out);
	ofstream ofs1("C:\\Users\\Desktop\\monte_carlo.txt", ios::out);
	ofstream tfs1("C:\\Users\\Desktop\\matrix.txt", ios::out);
	ofstream fs2("C:\\Users\\Desktop\\n25.txt", ios::out);
	ofstream ofs2("C:\\Users\\Desktop\\monte_carlo_25.txt", ios::out);
	ofstream tfs2("C:\\Users\\Desktop\\matrix_25.txt", ios::out);
	ofstream fs3("C:\\Users\\Desktop\\n50.txt", ios::out);
	ofstream ofs3("C:\\Users\\Desktop\\monte_carlo_50.txt", ios::out);
	ofstream tfs3("C:\\Users\\Desktop\\matrix_50.txt", ios::out);

	MATRIX N_FIVE(5);
	N_FIVE.print(fs1);
	N_FIVE.monte_carlo();
	N_FIVE.print_steady(ofs1);
	N_FIVE.matrix_multiplication();
	N_FIVE.print_steady(tfs1);
	MATRIX N_25(25);
	N_25.print(fs2);
	N_25.monte_carlo();
	N_25.print_steady(ofs2);
	N_25.matrix_multiplication();
	N_25.print_steady(tfs2);
	MATRIX N_50(50);
	N_50.print(fs3);
	N_50.monte_carlo();
	N_50.print_steady(ofs3);
	N_50.matrix_multiplication();
	N_50.print_steady(tfs3);

	N_FIVE.absorbtion(); //Making one of states of the matrices absorbing
	N_25.absorbtion();
	N_50.absorbtion();

	N_FIVE.monte_carlo();
	N_FIVE.print_steady(ofs1);
	N_FIVE.matrix_multiplication();
	N_FIVE.print_steady(tfs1);

	N_25.monte_carlo();
	N_25.print_steady(ofs2);
	N_25.matrix_multiplication();
	N_25.print_steady(tfs2);

	N_50.monte_carlo();
	N_50.print_steady(ofs3);
	N_50.matrix_multiplication();
	N_50.print_steady(tfs3);

	return 0;

}