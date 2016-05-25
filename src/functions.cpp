/******************************************************************************
  *  Name:    createM
  *  Autor:   Gregor Danler
  *  Datum:   25.05.2015
  *  Version: 1
  *****************************************************************************
  *
  * Kurzbeschreibung:
  * *****************
  *
  * Matrixoperationen für (n x n)- Matrizen
  * 
  * Konventionen: 
  *  - Zeilenindex i  mit i = 0,1, ... n-1
  *  - Spaltenindex j mit j = 0,1, ... n-1
  *  - A_ij :=A[i*n +j]
  * 
**/

/// Includefiles:
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <cmath>

using namespace std;

// Unterprogramme:
// ***************

// Erzeugen der Nullmatrix
void mkNULLM (double A[], int n){
  int k;
  for (int i=0; i<n; i++){
    k=i*n;
    for (int j=0; j<n; j++){
      A[k+j]=0;
    }
  }
}

// Erzeugen der Einsmatrix
void mkEinsM (double A[], int n){
  mkNULLM(A, n);
  for (int i=0; i<n; i++) A[i*n+i]=1;
}

//Erzeugen einer Zufallsmatrix
void mkRandM (double A[], int n){
  srand (time(NULL));
  int k;
  for (int i=0; i<n; i++){
    k=i*n;
    for (int j=0; j<n ; j++) A[k+j]=(rand() % 10) - 5;
  }
}

// Ausgabe einer Matrix
void printM (char X, double A[], int n){
  cout << X << "= " << endl;
    for(int i=0; i<n; i++){
      int k =i*n;
      for(int j=0; j<n; j++){
        cout <<"\t"<<  A[k+j];
      }
      cout << endl;
   }
}

// Betragssumme einer Matrixzeile i = 0, ... , n-1
double SumAbsRowM (double A[], int n, int i){
 double ret=0;
 int k=i*n;
 for (int j=0;j<n;j++) ret += abs(A[k+j]);
 return ret;
}

// Vergleich des Betrag eines Diagonalelement mit dem Zeilensummenbetrag
int ChckDiagM (double A[], int n){
  int ret=0;
  double delm;
  for(int i=1; i<n ; i++){
    int k=i*n;
    delm = abs(A[k+i]);
    if (delm  > (SumAbsRowM(A, n, i) - delm )){
      ret=1;
    } else {
      ret=0;
      i=n;
    }
  }
  return ret;
}

// 
//
// Programmstart:
// **************
int main(int argc, char *argv[]) {
  int n=100;
  int l=n*n;
  double A[l], B[l];
  mkNULLM (A, n);
  mkEinsM (B, n);
  cout << "Nullmatrix: " << ChckDiagM(A,n) << endl
       << "Einsmatrix: " << ChckDiagM(B,n) << endl;
}
