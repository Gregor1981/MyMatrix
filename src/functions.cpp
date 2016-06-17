/******************************************************************************
* Name: createM
* Autor: Gregor Danler
* Datum: 10.06.2015
* Version: 1
*****************************************************************************
*
* Kurzbeschreibung:
* *****************
*
*
* Hinweise:
* - Zeilenindex i mit i = 0,1, ... n-1
* - Spaltenindex j mit j = 0,1, ... n-1
* - A_ij :=A[i*n +j]
* - Noch zu testen: max n (Derzeit n=1022)
* 
* Matrixoperationen für (n x n)- Matrizen:
* ****************************************
*
* 1. Erzeugung von Matrizen:
*
* - void mkNULLM (double A[],int n)
* 	Erzeugt eine Nullmatrix
*
* - void mkEinsM (double A[],int n)
* 	Erzeugt die Einheitsmatrix
*
* -	void mkSteigM (double A[], int n)
*	Erzeugen einer Matrix mit fortlaufend steigenden Wert
*
* -	void mkLogM (double A[], int n)
*	Erzeugen einer Matrix mit fortlaufend logarithmisch steigenden Wert
*
* -	void mkRandM (double A[],int n)
*	Erzeugen einer Zufallsmatrix Werte -5 bis 5
*
* 
* 2. Ausgabe:
* 
* -	void printM (char X, double A[], int n, int m)
*	Ausgabe einer (n x m)-Matrix
*
* -	double SumAbsRowM (double A[], int n, int i)
*	Betragssumme einer Matrixzeile i = 0, ... , n-1
*
*
*
* 3. Manipulatinen:
*
* -	void ChangeRow (double A[], int n, int m,int r1, int r2)
*	Zeilentausch einer (nxm)-Matrix: Zeile r1  mit r2
*
* -	int SetMaxPivot(double A[], int n, int ps) 
*	Pivotelement mit maximalen Eintrag erzeugen
*	Ausgabewert: 1 bei Zeilentausch sonst 0 (Wichtig für Determinante)
*	Parameter ps: Pivotschritt
*
* -	void PivotSet(double A[], int n, int ps)
*	Pivotschritt: 
*		Definiert l_ij = a_ij/a_jj 			mit i = ps+1, ps+2, ... , n und j = ps+1
*		Berechnet a_ik = a_ik -l_ij*a_jk 	mit k = ps+1, ps+2, ... , n
*		Setzt 	  a_ij = l_ij	
*
* -	double Gauss(double A[], int n)
*	Gaussalgorithmus: LR-Zerlegung
*
*
*
* 	
**/


#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iomanip> //std::setprecision

using namespace std;

// Unterprogramme:
// ***************


// Funktionen R -> R:

double MyABS(double v);


// Erzeugung von speziellen Matrizen
void mkNULLM (double A[],int n);
void mkEinsM (double A[],int n);
void mkDiag1M (double A[],int n);
void mkSteigM (double A[], int n);
void mkLogM (double A[], int n);
void mkRandM (double A[],int n);


// Ausgabefunktionen
void printM (char X, double A[], int n, int m);
double SumAbsRowM (double A[], int n, int i);

// Matrixmanipulationen
void ChangeRow (double A[], int n, int m, int r1, int r2);
int SetMaxPivot(double A[], double b[], int n, int ps);
void PivotSet(double A[], double b[], int n, int ps);
double Gauss(double A[], double b[], int n);
void EvalX(double A[], double b[], int n , double x[]);

// Auf schwache diagonale Dominanz kontrollieren
int ChckDiagM (double A[], int n);



//
// Programmstart:
// **************
int main(int argc, char *argv[]) {
  int n=1022;
  int laenge=n*n;
  double A[laenge];
  double b[n];
  double x[n];
  for (int i=0; i<n; i++) x[i]=1;  
  for (int i=0; i<n; i++) b[i]=i;
  //mkRandM (A, n);
  mkDiag1M (A, n);
  printM('x',x,n,1);
 // printM('A',A,n,n);
  cout << "nu=" << Gauss(A,b,n) << endl;
  EvalX(A,b,n,x);
 // printM('A',A,n,n);
 // printM('b',b,n,1);
  printM('x',x,n,1);
 return 0;
}



//Funktionsdeklarationen
// Auswerten von x
void EvalX(double A[], double b[], int n , double x[]){
 	for (int i=(n-1); i>=0; i--){
 		double sum = 0;
 		for(int j=(i+1); j<n;j++){
 			sum=sum+A[i*n+j]*x[j]; 
 		}
 		x[i]=(b[i]-sum)/A[i*n+i];
 	}
}


//Algorithmus
double Gauss(double A[], double b[], int n){
	int nu=0;
	for(int l=0; l<(n-1); l++){
		nu=nu + SetMaxPivot(A,b,n,l);
		PivotSet(A,b,n,l);
	}
	return nu;
}

// Pivotschritt
void PivotSet(double A[], double b[], int n, int ps){
	double l_ips;
	int k;
	for(int i=ps+1; i<n; i++){
		k    = i*n;
		l_ips = A[k+ps]/A[n*ps+ps]; 
		b[i] = b[i]-l_ips*b[ps];
		A[k+ps]=l_ips;
		int k2=n*ps;
		for(int j=ps+1; j<n; j++){
			A[k+j] = A[k+j] - l_ips*A[k2+j];
		}
	}
}

// Pivotelement mit maximalen Eintrag erzeugen
// Liefert Wert 1 bei Zeilentausch sonst 0 (Wichtig für Determinante)
int SetMaxPivot(double A[], double b[], int n, int ps){
	int indPind=ps;
	int maxPind=ps;
	double indPval=MyABS(A[ps*n+ps]);

    for (int i=ps+1; i<n; i++){
		if (MyABS(A[n*i+ps]) > indPval)
		{ 
			maxPind = i;
			indPval = MyABS(A[n*i+ps]);
		}
	}
    if (indPind != maxPind){
		ChangeRow(A, n, n, indPind, maxPind);
		ChangeRow(b, n, 1, indPind, maxPind);
		return 1;
	}
	return 0;		
}

// Zeilentausch: r1 mit r2
void ChangeRow (double A[], int n, int m, int r1, int r2){
        double Z;
        int k = r1*n;
        int l = r2*m;
        for (int j=0; j<m ; j++) {
                Z=A[k+j];
                A[k+j]=A[l+j];
                A[l+j]=Z;
        }
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

// Betragssumme einer Matrixzeile i = 0, ... , n-1
double SumAbsRowM (double A[], int n, int i){
 double ret=0;
 int k=i*n;
 for (int j=0;j<n;j++) ret += MyABS(A[k+j]);
 return ret;
}

// Ausgabe einer Matrix
void printM (char X, double A[], int n, int m){
  cout << X << "= " << endl;
  for(int i=0; i<n; i++){
    int k =i*m;
    for(int j=0; j<m; j++){
       cout << "\t" << setprecision(4) << A[k+j];
    }
    cout << endl;
  }
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

// Erzeugen einer Matrix mit fortlaufend logarithmisch steigenden Wert
void mkLogM (double A[], int n){
        int k;
        int eintrag=1;
  for (int i=0; i<n; i++){ k=i*n;
    for (int j=0; j<n; j++){
      A[k+j]=log1p(eintrag);
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

// Erzeugen der Einsmatrix
void mkEinsM (double A[],int n){
  mkNULLM(A,n);
  for (int i=0; i<n; i++) A[i*n+i]=1;
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
// Erzeugen einer Diagonalmatrix mit einsen in der Oberen Diagonale
void mkDiag1M (double A[],int n){
	mkNULLM(A,n);
	for(int i=0; i<n; i++){
		for(int j=i; j<n; j++)
			A[i*n+j]=1;
	}
}
// Absolutbetrag
double MyABS(double v){
        if (v>=0)
                return v;
        else
                return -v;
}
