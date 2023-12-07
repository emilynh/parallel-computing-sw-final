#include "mat.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#define M 3      // match
#define MM -3    // mismatch
#define W -2     // gap score
#define max(a, b) (((a) > (b)) ? (a) : (b)) // return maximum of two values
#define min(a, b) (((a) < (b)) ? (a) : (b)) // return minimum of two values



void read_sequence_from_file(const std::string& filename, std::vector<char>& seq, int line_number) {
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        for (int i = 0; i < line_number; ++i) {
            if (!std::getline(file, line)) {
                std::cerr << "Error reading line " << line_number << " from file: " << filename << std::endl;
                file.close();
                return;
            }
        }
        seq.assign(line.begin(), line.end());
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

std::pair<int, int> fill_cpu(Matrix h, Matrix d, char seqA[], char seqB[]) {

  int full_max_id = 0;
  int full_max_val = 0;

  for (int i = 1; i < h.height; i++) {
    for (int j = 1; j < h.width; j++) {

      // scores
      int max_score = 0;
      int direction = 0;
      int tmp_score;
      int sim_score;

      // comparison positions
      int id = i * h.width + j;                  // current cell
      int abov_id = (i - 1) * h.width + j;       // above cell, 1
      int left_id = i * h.width + (j - 1);       // left cell, 2
      int diag_id = (i - 1) * h.width + (j - 1); // upper-left diagonal cell, 3

      // above cell
      tmp_score = h.elements[abov_id] + W;
      if (tmp_score > max_score) {
        max_score = tmp_score;
        direction = 1;
      }

      // left cell
      tmp_score = h.elements[left_id] + W;
      if (tmp_score > max_score) {
        max_score = tmp_score;
        direction = 2;
      }

      // diagonal cell (preferred)
      char baseA = seqA[j - 1];
      char baseB = seqB[i - 1];
      if (baseA == baseB) {
        sim_score = M;
      } else {
        sim_score = MM;
      }

      tmp_score = h.elements[diag_id] + sim_score;
      if (tmp_score >= max_score) {
        max_score = tmp_score;
        direction = 3;
      }

      // assign scores and direction
      h.elements[id] = max_score;
      d.elements[id] = direction;

      if (max_score > full_max_val) {
        full_max_id = id;
        full_max_val = max_score;
      }
    }
  }
  return std::make_pair(full_max_id, full_max_val);
}


void traceback(Matrix d, int max_id, char seqA[], char seqB[],
               std::vector<char> &seqA_aligned,
               std::vector<char> &seqB_aligned) {

  int max_i = max_id / d.width;
  int max_j = max_id % d.width;

  // traceback algorithm from maximum score to 0
  while (max_i > 0 && max_j > 0) {

    int id = max_i * d.width + max_j;
    int dir = d.elements[id];

    switch (dir) {
    case 1:
      --max_i;
      seqA_aligned.push_back('-');
      seqB_aligned.push_back(seqB[max_i]);
      break;
    case 2:
      --max_j;
      seqA_aligned.push_back(seqA[max_j]);
      seqB_aligned.push_back('-');
      break;
    case 3:
      --max_i;
      --max_j;
      seqA_aligned.push_back(seqA[max_j]);
      seqB_aligned.push_back(seqB[max_i]);
      break;
    case 0:
      max_i = -1;
      max_j = -1;
      break;
    }
  }
}

// print aligned sequnces
void io_seq(std::vector<char> &seqA_aligned, std::vector<char> &seqB_aligned) {

  std::cout << "Aligned sub-sequences of A and B: " << std::endl;
  int align_len = seqA_aligned.size();
  std::cout << "   ";
  for (int i = 0; i < align_len + 1; ++i) {
    std::cout << seqA_aligned[align_len - i];
  }
  std::cout << std::endl;

  std::cout << "   ";
  for (int i = 0; i < align_len + 1; ++i) {
    std::cout << seqB_aligned[align_len - i];
  }
  std::cout << std::endl << std::endl;
}



void smith_water_cpu(Matrix h, Matrix d, char seqA[], char seqB[]) {

  // populate scoring and direction matrix and find id of max score
  std::pair<int, int> result = fill_cpu(h, d, seqA, seqB);
  int max_id = result.first;
  int score = result.second;
  // traceback
  std::vector<char> seqA_aligned;
  std::vector<char> seqB_aligned;
  traceback(d, max_id, seqA, seqB, seqA_aligned, seqB_aligned);

  // print aligned sequences
  io_seq(seqA_aligned, seqB_aligned);

  std::cout << std::endl;
  std::cout << "CPU result: " << std::endl;
  std::cout << "Max score:" << score << std::endl;
}

char* vectorToCharArray(const std::vector<char>& vec) {
    char* arr = new char[vec.size() + 1]; // +1 để thêm ký tự null
    std::copy(vec.begin(), vec.end(), arr);
    arr[vec.size()] = '\0'; // thêm ký tự null vào cuối mảng
    return arr;
}

int main() {
  std::vector<char> seqA;
  read_sequence_from_file("D:\\Gpu-SW\\src\\a.txt", seqA, 1);
  std::vector<char> seqB;
  read_sequence_from_file("D:\\Gpu-SW\\src\\b.txt", seqB, 1);
  std::cout << "Seq A with length " << seqA.size() << " is: ";
  for (int i = 0; i < seqA.size(); i++)
    std::cout << seqA[i];
  std::cout << std::endl;
  std::cout << "Seq B with length " << seqB.size() << " is: ";
  for (int i = 0; i < seqB.size(); i++)
    std::cout << seqB[i];
  std::cout << std::endl;
  char* arrA = vectorToCharArray(seqA);
  char* arrB = vectorToCharArray(seqB);
  // initialize scoring and direction matrices
  Matrix scr_cpu(seqA.size() + 1, seqB.size() + 1); // cpu score matrix
  Matrix dir_cpu(seqA.size() + 1, seqB.size() + 1); // cpu direction
  Matrix scr_gpu(seqA.size() + 1, seqB.size() + 1); // gpu score matrix
  Matrix dir_gpu(seqA.size() + 1, seqB.size() + 1); // gpu direction matrix

  // apply initial condition of 0
  for (int i = 0; i < scr_cpu.height; i++) {
    for (int j = 0; j < scr_cpu.width; j++) {
      int id = i * scr_cpu.width + j;
      scr_cpu.elements[id] = 0;
      dir_cpu.elements[id] = 0;
      scr_gpu.elements[id] = 0;
      dir_gpu.elements[id] = 0;
    }
  }

  // visualize initial scoring matrix
  // io_score(std::string("init.dat"), scr_cpu, seqA, seqB);

  // CPU
  auto start_cpu = std::chrono::steady_clock::now();
  smith_water_cpu(scr_cpu, dir_cpu, arrA, arrB); // call CPU smith water
  auto end_cpu = std::chrono::steady_clock::now();
  auto diff = end_cpu - start_cpu;
  std::cout << "   CPU time = "
            << std::chrono::duration<double, std::milli>(diff).count() << " ms"
            << std::endl;
  std::cout << std::endl;

  // GPU
  // smith_water_gpu(scr_gpu, dir_gpu, seqA, seqB); // call GPU smith water
  // deallocate memory
  scr_cpu.cpu_deallocate();
  dir_cpu.cpu_deallocate();
  scr_gpu.cpu_deallocate();
  dir_gpu.cpu_deallocate();

  return 0;
}

