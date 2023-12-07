#include "mat.h"
#include <chrono>
#include <fstream>
#include <ios>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>
#include <thread>
#include <mutex>
#include <algorithm>


#define M 3      // match
#define MM -3    // mismatch
#define W -2     // gap score
#define A_LEN 16 // 16, 32, 64, 256, 1024, 2048, 8192  len of sequence A
#define B_LEN 16 // 16, 32, 64, 256, 1024, 2048, 8192 len of sequence B
#define max(a, b) (((a) > (b)) ? (a) : (b)) // return maximum of two values
#define min(a, b) (((a) < (b)) ? (a) : (b)) // return minimum of two values

// Forward declarations of scoring kernel
__global__ void fill_gpu(Matrix h, Matrix d, char seqA[], char seqB[],
                         const int *k);

// generate random sequence of length n
// Thay thế hàm read_sequence_from_file
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

__global__ void fill_gpu(Matrix h, Matrix d, char seqA[], char seqB[],const int k, int max_id_val[]) {

  // scores
  int max_score = 0;
  int direction = 0;
  int tmp_score;
  int sim_score;

  // row and column index depending on anti-diagonal
  int i = threadIdx.x + 1 + blockDim.x * blockIdx.x;
  if (k > A_LEN + 1) {
    i += (k - A_LEN);
  }
  int j = ((k) - i) + 1;

  // comparison positions
  int id = i * h.width + j;
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

  // similarity score for diagonal cell
  char baseA = seqA[j - 1];
  char baseB = seqB[i - 1];
  if (baseA == baseB) {
    sim_score = M;
  } else {
    sim_score = MM;
  }

  // diagonal cell (preferred)
  tmp_score = h.elements[diag_id] + sim_score;
  if (tmp_score >= max_score) {
    max_score = tmp_score;
    direction = 3;
  }

  // assign scores and direction
  h.elements[id] = max_score;
  d.elements[id] = direction;

  // save max score and position
  if (max_score > max_id_val[1]) {
    max_id_val[0] = id;
    max_id_val[1] = max_score;
  }
}

// traceback: starting at the highest score and ending at a 0 score
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

// input output function to visualize matrix
void io_score(std::string file, Matrix h, char seqA[], char seqB[]) {
  std::ofstream myfile_tsN;
  myfile_tsN.open(file);

  // print seqA
  myfile_tsN << '\t' << '\t';
  for (int i = 0; i < A_LEN; i++)
    myfile_tsN << seqA[i] << '\t';
  myfile_tsN << std::endl;

  // print vertical seqB on left of matrix
  for (int i = 0; i < h.height; i++) {
    if (i == 0) {
      myfile_tsN << '\t';
    } else {
      myfile_tsN << seqB[i - 1] << '\t';
    }
    for (int j = 0; j < h.width; j++) {
      myfile_tsN << h.elements[i * h.width + j] << '\t';
    }
    myfile_tsN << std::endl;
  }
  myfile_tsN.close();
}

void smith_water_cpu(Matrix h, Matrix d, char seqA[], char seqB[]) {

  // populate scoring and direction matrix and find id of max score
  std::pair<int, int> result = fill_cpu(h, d, seqA, seqB);
  int max_id = result.first;
  // traceback
  std::vector<char> seqA_aligned;
  std::vector<char> seqB_aligned;

  // print aligned sequences
  io_seq(seqA_aligned, seqB_aligned);

  std::cout << std::endl;
  std::cout << "CPU result: " << std::endl;

  // print cpu populated direction and scoring matrix
  io_score(std::string("score.dat"), h, seqA, seqB);
  io_score(std::string("direction.dat"), d, seqA, seqB);
}

void smith_water_gpu(Matrix h, Matrix d, char seqA[], char seqB[]) {
    // allocate and transfer sequence data to device
    char *d_seqA, *d_seqB;
    cudaMalloc(&d_seqA, strlen(seqA) * sizeof(char));
    cudaMalloc(&d_seqB, strlen(seqB) * sizeof(char));
    cudaMemcpy(d_seqA, seqA, strlen(seqA) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_seqB, seqB, strlen(seqB) * sizeof(char), cudaMemcpyHostToDevice);

    // initialize matrices for gpu
    int Gpu = 1;
    Matrix d_h(strlen(seqA) + 1, strlen(seqB) + 1, Gpu);
    Matrix d_d(strlen(seqA) + 1, strlen(seqB) + 1, Gpu);
    d_h.load(h, Gpu);
    d_d.load(d, Gpu);

    // max id and value
    int *d_max_id_val;
    std::vector<int> h_max_id_val(2, 0);
    cudaMalloc(&d_max_id_val, 2 * sizeof(int));
    cudaMemcpy(d_max_id_val, h_max_id_val.data(), 2 * sizeof(int), cudaMemcpyHostToDevice);

    cudaEvent_t start, stop;
    float time;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    // loop over diagonals of the matrix
    for (int i = 1; i <= (strlen(seqA) + strlen(seqB) - 1); i++) {
      // i là chỉ số đường chéo, ma trận đang có max là (strlen(seqA) + strlen(seqB) - 1) đường chéo
        int col_idx = max(0, (i - strlen(seqB)));
        // biến này là chỉ số cột đầu tiên của một đường chéo
        int diag_len = min(i, (strlen(seqA) - col_idx));
        // biến này là chiều dài của đường chéo
        // launch the kernel: one block by length of diagonal
        int blks = 1; //1 block
        dim3 dimBlock(diag_len / blks); //kích thước của mỗi block, ở đây mỗi block có số thread bằng số phần tử có trong đường chéo ngược
        dim3 dimGrid(blks);
        fill_gpu<<<dimGrid, dimBlock>>>(d_h, d_d, d_seqA, d_seqB, i, d_max_id_val);
        cudaDeviceSynchronize();
    }

    // copy data back
    size_t size = (strlen(seqA) + 1) * (strlen(seqB) + 1) * sizeof(float);
    cudaMemcpy(h.elements, d_h.elements, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(d.elements, d_d.elements, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_max_id_val.data(), d_max_id_val, 2 * sizeof(int), cudaMemcpyDeviceToHost);

    // std::cout << "   Max score of " << h_max_id_val[1] << " at " << h_max_id_val[0] << std::endl;
    int max_id = h_max_id_val[0];
    std::vector<char> seqA_aligned;
    std::vector<char> seqB_aligned;
    traceback(d, max_id, seqA, seqB, seqA_aligned, seqB_aligned);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);

    // visualize output
    // io_seq(seqA_aligned, seqB_aligned);
    // io_score(std::string("score_gpu.dat"), h, seqA, seqB);
    // io_score(std::string("direction_gpu.dat"), d, seqA, seqB);

    // deallocate memory
    d_h.gpu_deallocate();
    d_d.gpu_deallocate();
    cudaFree(d_seqA);
    cudaFree(d_seqB);
    cudaFree(d_max_id_val);
}


char* vectorToCharArray(const std::vector<char>& vec) {
    char* arr = new char[vec.size() + 1]; // +1 để thêm ký tự null
    std::copy(vec.begin(), vec.end(), arr);
    arr[vec.size()] = '\0'; // thêm ký tự null vào cuối mảng
    return arr;
}

int countLines(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
      std::cerr << "Unable to open file: " << filename << std::endl;
      return -1; // Return -1 to indicate an error
  }

  int lineCount = 0;
  std::string line;
  while (std::getline(file, line)) {
      ++lineCount;
  }

  file.close();
  return lineCount;
}
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream tokenStream(s);
    std::string token;
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::mutex mutex;

// ... (Không thay đổi các hàm khác)

void process_and_compare_line(const std::string& line, std::vector<char>& seqA, std::vector<std::pair<int, std::string>>& scores) {
    // Split the line into gene sequence and number
    std::vector<std::string> parts = split(line, '\t');

    if (parts.size() == 2) {
        std::vector<char> seqB(parts[0].begin(), parts[0].end());

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

        // CPU
        std::pair<int, int> result = fill_cpu(scr_cpu, dir_cpu, arrA, arrB);
        int score = result.second;

        std::pair<int, std::string> element(score, parts[1]);
        scores.push_back(element);

        // GPU
        smith_water_gpu(scr_gpu, dir_gpu, arrA, arrB);

        // deallocate memory
        scr_cpu.cpu_deallocate();
        dir_cpu.cpu_deallocate();
        scr_gpu.cpu_deallocate();
        dir_gpu.cpu_deallocate();
    } else {
        std::cerr << "Invalid line format: " << line << std::endl;
    }
}
void read_and_compare_sequences_from_file(const std::string& filename, std::vector<char>& seqA, std::vector<std::pair<int, std::string>>& scores, int num_threads) {
  std::ifstream file(filename);
  if (!file.is_open()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      exit(1);
  }

  std::string line;
  std::vector<std::thread> threads;

  // Process each line along with seqA
  while (std::getline(file, line)) {
      threads.emplace_back(process_and_compare_line, line, std::ref(seqA), std::ref(scores));
  }

  // Wait for all threads to complete
  for (std::thread& t : threads) {
      t.join();
  }

  file.close();
}

int main() {
  std::vector<char> seqA;
  int length = countLines("D:\\Gpu-SW\\src\\dog.txt");
  std::cout << "Số chuỗi trong file tham chiếu:" << length << std::endl;
  //Read seqA from line 1
  read_sequence_from_file("D:\\Gpu-SW\\src\\a.txt", seqA, 1);

  std::vector<std::pair<int, std::string>> scores;
  const int num_threads = 2; // Get the number of available threads
  auto start_time = std::chrono::high_resolution_clock::now();

  read_and_compare_sequences_from_file("D:\\Gpu-SW\\src\\dog.txt", seqA, scores, num_threads);
  std::cout << "Số điểm tương đồng được tính:" << scores.size() << std::endl;
  auto end_time = std::chrono::high_resolution_clock::now();

  // Sort the scores list
  std::sort(scores.begin(), scores.end(), [](const auto& a, const auto& b) {
    return a.first < b.first;  // Compare based on the score (first element of the pair)
  });
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
  std::cout << "Thời gian chạy: " << duration.count() << " microseconds" << std::endl;
  // std::cout << "Scores list: ";

  int maxScore = -1; // Initialize maxScore to a value that is guaranteed to be less than any actual score
  std::string maxNumber;
  for (const auto& element : scores) {
    int score = element.first;
    std::string number = element.second;
    // std::cout << "(" << score << "," << number << ")" << std::endl;
    if (score > maxScore){
      maxScore = score;
      maxNumber = number;
    }
  }
  std::cout << maxNumber << std::endl;

}
