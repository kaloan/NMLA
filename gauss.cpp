#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace std;

template <typename T>
using matrix = vector<vector<T>>;

template <typename T>
int maxOfColumnUnder(const matrix<T>& A, const int& k)
{
  int maxIndex = k;
  for (int j=k+1; j<A.size(); ++j)
  {
    if(abs(A[maxIndex][k]) < abs(A[j][k])) maxIndex = j;
  }
  return maxIndex;
}

template<typename T>
vector<T> GaussStable(matrix<T> A, vector<T> b) {
  if(A.size() == 0 || A.size() != A[0].size() || A.size() != b.size()){
    throw invalid_argument("Matrix is not square or doesn't have as many rows as right hand side");
  }
  for(int k = 0; k < A.size(); ++k) {
    auto maxPivotRowIndex = maxOfColumnUnder(A, k);
    if(maxPivotRowIndex > k) {
      swap(A[k], A[maxPivotRowIndex]);
      swap(b[k], b[maxPivotRowIndex]);
    }

    b[k] = b[k] / A[k][k];
    for (int s = A.size() - 1; s >= k; --s){
      A[k][s] = A[k][s] / A[k][k];
    }

    for(int j = k + 1; j < A.size(); ++j) {
      b[j] = b[j] - A[j][k] * b[k];
      for (int s = A.size() - 1; s >= k; --s){
        A[j][s] = A[j][s] - A[j][k] * A[k][s];
      }
    }
  }

  for(int k = A.size() - 1; k > 0; --k) {
    for(int j = k - 1; j >= 0; --j) {
      b[j] = b[j] - A[j][k] * b[k];
      A[j][k] = A[j][k] - A[j][k] * A[k][k];
    }
  }

  return b;
}

int main ()
{
  // Solution of Ax=b is {8/3, -16/3, 3}
  matrix<double> A = {
    {1, 2, 3},
    {4, 5, 6},
    {0, 0, 1}
  };
  vector<double> b = {1, 2, 3};
  auto x = GaussStable(A, b);
  for (const auto& xi : x)
  {
    cout << xi << '\n';
  }
  return 0;
}
