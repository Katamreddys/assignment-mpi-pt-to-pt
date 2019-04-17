#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void generateReduceData (int* arr, size_t n) {
  for (size_t i=0; i<n; ++i) {
    if (i%2 ==0 )
      arr[i] = i/2;
    else
      arr[i] = -(i/2);
  }
}


void generateMergeSortData (int* arr, size_t n) {
  for (size_t i=0; i<n; ++i) {
    arr[i] = i;
  }

  for (int i = n-1; i>0; --i) {
    int swapwith = rand()%i;
    int temp = arr[i];
    arr[i] = arr[swapwith];
    arr[swapwith] = temp;
  }  
}


void generatePrefixSumData (int* arr, int n) {
  for (int i=0; i<n; ++i) {
    if (i%2 == 0 )
      arr[i] = i/2;
    else
      arr[i] = -(i/2);
  }
}


void checkPrefixSumResult (const int* arr, size_t n) {
  int correct = 1;

  int *data = malloc(sizeof(int)*n);
  generatePrefixSumData (data, n);

  correct = correct && (arr[0] == 0);

  int sum = 0;
  for (int i=0; i<n; ++i) {
    sum += data[i];
    if (arr[i+1] != sum) {
      correct = 0;
    }
  }
  
  if (correct) {  
    printf ("checked\n");
  }
  else {
    printf ("incorrect\n");
  }
  free(data);
}


//this is based on how MS data is generated
void checkMergeSortResult (const int* arr, size_t n) {
  int correct = 1;
  
  for (int i=0; i< n-1; ++i) {
    correct = correct && (arr[i] == i);
  }

  if (correct) {  
    printf ("checked\n");
  }
  else {
    printf ("incorrect\n");
  }
}


// assumes n and m are powers of 10, 
//
void generateLCS(char* X, int m, char* Y, int n) {

  const char A = 'A'; 
  const char B = 'B';
  const char C = 'C'; 
  
  int x = (n>m) ? m : n;

  // X ~ AABBAAAAAA
  for (int i=0; i<(2*x/10); ++i) {
    X[i] = A;
    X[i+(2*x/10)] = B;
  }
  for (int i=(4*x/10); i<m; ++i) {
    X[i] = A;
  }

  // Y ~ CCCCCCCCBB
  for (int i=0; i<(n-(2*x/10)); ++i) {
    Y[i] = C;
  }
  for (int i=n-(2*x/10); i<n; ++i) {
    Y[i] = B;
  }
}


void checkLCS(char* X, int m, char* Y, int n, int result) {
  int x = (n>m) ? m : n;

  if ((x*2/10) == result) {
    printf ("checked\n");
  }
  else {
    printf ("incorrect\n");
  }
}


double generate2DHeat(long n, long i, long j) {
  if (j == n/2)
    return 1.0;
  else
    return 0.0;
}

void generate2DHeatOld(double** H, long n, long rank, long P) {
  
  long sqrtP = sqrt(P);
  long sub_n = n/sqrtP;

  long mid = (P%2) ? sub_n/2 : 0;

  for (long i=0; i<sub_n; ++i) {
    for (long j=0; j<sub_n; ++j) {
      H[i][j] = 0.0;
    }

    if (rank%sqrtP == sqrtP/2)
      H[i][mid] = 1.0;
  }
}

int not_equal(double a, double b) {
  return fabs(a-b) > 0.0000001;
}

//this function return 1 on correct. 0 on incorrect. Note that it may return 1 on incorrect. But a return of 0 means it is definitely incorrect
int check2DHeat(long n, long i, long j, double v, long k) {
  long mid = n/2;

  //don't deal with the wave bouncing on the border
  //so return true in these cases
  if (k > n/2 -1)
    return 1;

  long dist = mid-j;
  if (dist < 0)
    dist = -dist;

  //only deal with the front of the wave
  //so return true in all other cases
  if (k != dist)
    return 1;

  double correct_value = pow(0.2, k);
  return !(not_equal (v, correct_value));
}


int check2DHeatOld(double** H, long n, long rank, long P, long k) {
  long len = (k+1)*2 + 1;
  double** sequence = malloc (2*sizeof(double*));
  sequence[0] = malloc (sizeof(double)*len);
  sequence[1] = malloc (sizeof(double)*len);
  
  long mid = len/2;

  
  long curr = 0;
  long next = 1;

  // init sequence
  for (long i=0; i<len; ++i) { 
    sequence[0][i] = (i==mid);
    sequence[1][i] = (i==mid);
  }

  // calculate sequence up to point of symmetry
  for (long i=0; i<(mid-1); ++i) {
    sequence[next][mid] = ( sequence[curr][mid-1]*2 + sequence[curr][mid]*3 ) / 5.0;
    for (long j=mid-1; j>0; --j) {
      sequence[next][j] = ( sequence[curr][j-1]+ sequence[curr][j]*3+sequence[curr][j+1] ) / 5.0;
    }
    curr = (curr+1)%2;
    next = (next+1)%2;
  }

  // copy front into back 
  for (long i=1; i<mid; ++i) {
    sequence[curr][len-i-1] = sequence[curr][i];
  }
  
  long ret = 0;
  
  // where am i?
  if (rank%(long)sqrt(P) == (long)sqrt(P)/2 || rank%(long)sqrt(P) == (long)(sqrt(P)/2)-1 ) { 
    long sub_n = n/sqrt(P);

    long start = (P%2) ? sub_n/2 : sub_n;

    long foreward = rank%(long)sqrt(P) == sqrt(P)/2;
    long backward = P%2 || (rank%(long)sqrt(P) == (long)(sqrt(P)/2)-1);

    // check
    for (long i=0; i<mid-1; ++i) {
      if (foreward && not_equal(H[0][(start+i)%sub_n], sequence[curr][mid+i])) {
        ret = -1; goto end;
      } 

      if (i != 0 && backward && not_equal(H[0][start-i], sequence[curr][mid-i])) {
	ret = -1; goto end;
      }
    }
  }

 end:
  free(sequence[0]);
  free(sequence[1]);
  free(sequence);

  return ret;

}
