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

// Dimension der Matrix
//int dimof(double A[]){
//  return sqrt(sizeof(A)/sizeof(A[0]));
//}

// Absolutbetrag 
double MyABS(double v){
	if (v>=0) 
		return v;
	else
		return -v;
}


// Erzeugen der Nullmatrix
void mkNULLM (double A[],int n){
  int k;
  for (int i=0; i<n; i++){ k=i*n;
    for (int j=0; j<n; j++){
      A[k+j]=0;
    }
  }
}

// Erzeugen der Einsmatrix
void mkEinsM (double A[],int n){
  mkNULLM(A,n);
  for (int i=0; i<n; i++) A[i*n+i]=1;
}

// Erzeugen einer Matrix mit fortlaufend steigenden Wert
void mkSteigM (double A[], int n){
	int k;
	int eintrag=1;
  for (int i=0; i<n; i++){ k=i*n;
    for (int j=0; j<n; j++){
      A[k+j]=eintrag;
			eintrag++;
    }
  }
}

//Erzeugen einer Zufallsmatrix
void mkRandM (double A[],int n){
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
       cout << "\t" << A[k+j];
    }
    cout << endl;
  }
}

// Betragssumme einer Matrixzeile i = 0, ... , n-1
double SumAbsRowM (double A[], int n, int i){
 double ret=0;
 int k=i*n;
 for (int j=0;j<n;j++) ret += MyABS(A[k+j]);
 return ret;
}


// Auf schwache diagonale Dominanz kontrollieren
int ChckDiagM (double A[], int n){
  int ret=0;
  double delm;
  for(int i=1; i<n ; i++){
    int k=i*n;
    delm = MyABS(A[k+i]);
    if (delm  >= (SumAbsRowM(A, n, i) - delm )){
      ret=1;
    } else {
      ret=0;
      i=n;
    }
  }
  return ret;
}

// Zeilentausch: r1 mit r2
void ChangeRow (double A[], int n, int r1, int r2){
	double Z;
	int k = r1*n;
	int l = r2*n;
	for (int j=0; j<n ; j++) {
		Z=A[k+j];
		A[k+j]=A[l+j];
		A[l+j]=Z;
	}	
}

// Pivotmaximierungssuche
int ChckMaxPivot (double A[], int n){
	int ret=0;
	double piv=MyABS(A[0]);
	for (int i=1; i<n; i++){
		if(piv < A[n*i]) ret=1;
		//piv < A[n+i] ? ret=1;
	}
	return ret;
}

// SetMaxPivot
void SetMaxPivot(double A[], int n){
	int indPind=0;
	int maxPind=0;
	
	double indPval=MyABS(A[0]);

	for (int i=1; i<n; i++){
		if (MyABS(A[n*i]) > indPval) maxPind = i;
	}
	if (indPind != maxPind) ChangeRow(A, n, indPind, maxPind);
}



//
// Programmstart:
// **************
int main(int argc, char *argv[]) {
  int n=5;
  int l=n*n;
  double A[l], B[l];
  mkNULLM (A, n);
  mkSteigM (B, n);
	cout << MyABS(-3.1) << endl;
	SetMaxPivot(B,n);
	
	ChangeRow(B, n, 0, 4);
	cout << ChckMaxPivot(B, n) << endl;
	printM('B', B, n);
  return 0;
}
