#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace std;
using namespace boost::numeric::ublas;

template<class T>
bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
{
	typedef permutation_matrix<std::size_t> pmatrix;

	matrix<T> A(input);
	pmatrix pm(A.size1());

	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	inverse.assign(identity_matrix<T> (A.size1()));
	lu_substitute(A, pm, inverse);

	return true;
}

void printMenu() {
	cout << "MENU" << endl << endl;
	cout << "1 Start simplex" << endl;
	cout << "2 test 1" << endl;
	cout << "4 test 2" << endl;
	cout << "3 wyjscie" << endl << endl;
	cout << "Opcja: ";
}

int _tmain(int argc, _TCHAR* argv[])
{
	bool menu = true;
	while(menu) {
		printMenu();
		int zmienna;
		cin >> zmienna;

		switch(zmienna) {
			case 1: {
				system("cls");
				cout << "Podaj liczbe dostawocow: ";
				int providersNumber;
				cin >> providersNumber;
				int providersAndDeliver = providersNumber+1;
				cout << endl << "Podaj liczbe odbiorcow: ";
				int receiversNumber;
				cin >> receiversNumber;

				matrix<double> INFO_MATRIX (receiversNumber, providersAndDeliver);

				/*cout << endl << "Podaj zapotrzebowanie: " << endl;
				for (unsigned i = 0; i < INFO_MATRIX.size1 (); ++ i)
					for (unsigned j = 0; j < INFO_MATRIX.size2 (); ++ j) {
						if(j == providersNumber) {
							cout << "Podaj oferowana miesieczna wielkosc dostawy: ";
							cin >> INFO_MATRIX (i, j);
						} else {
						cout << "Odbiorcy " << i;
						cout << ", dostawcy " << j;
						cout << ": ";
						cin >> INFO_MATRIX (i, j);
						}
					}*/

				INFO_MATRIX(0,0) = 2;
				INFO_MATRIX(0,1) = 1;
				INFO_MATRIX(0,2) = 1000;
				INFO_MATRIX(1,0) = 4;
				INFO_MATRIX(1,1) = 4;
				INFO_MATRIX(1,2) = 3000;
				INFO_MATRIX(2,0) = 1.5;
				INFO_MATRIX(2,1) = 0;
				INFO_MATRIX(2,2) = 300;

				/*INFO_MATRIX(0,0) = 50;
				INFO_MATRIX(0,1) = 40;
				INFO_MATRIX(0,2) = 50;
				INFO_MATRIX(0,3) = 20;
				INFO_MATRIX(0,4) = 70;
				INFO_MATRIX(1,0) = 40;
				INFO_MATRIX(1,1) = 80;
				INFO_MATRIX(1,2) = 70;
				INFO_MATRIX(1,3) = 30;
				INFO_MATRIX(1,4) = 50;
				INFO_MATRIX(2,0) = 60;
				INFO_MATRIX(2,1) = 40;
				INFO_MATRIX(2,2) = 70;
				INFO_MATRIX(2,3) = 80;
				INFO_MATRIX(2,4) = 80;*/

				//dodanie zapotrzebowania
				std::vector < int > requestion;
				for (unsigned i = 0; i < providersNumber; ++ i) {
					cout << "Podaj miesieczne zapotrzebowanie dostawcy " << i << ": ";
					double temp;
					cin >> temp;
					requestion.push_back(temp);
				}

				for (unsigned i = providersNumber; i < providersNumber+INFO_MATRIX.size2(); ++ i) {
					requestion.push_back(0);
				}
				
				//tworzenie pierwszej tablicy simplexowej
				matrix<double> SIMPLEX_TABLE(receiversNumber, providersAndDeliver+receiversNumber);
				identity_matrix<double> IDENTITY_MATRIX (receiversNumber);

				for (unsigned i = 0; i < providersAndDeliver-1; ++ i) {
					column(SIMPLEX_TABLE, i) = column(INFO_MATRIX, i);
				}
				for (unsigned i = providersAndDeliver-1; i < SIMPLEX_TABLE.size2()-1; ++ i) {
					column(SIMPLEX_TABLE, i) = column(IDENTITY_MATRIX, i-providersAndDeliver+1);
				}
				column(SIMPLEX_TABLE, SIMPLEX_TABLE.size2()-1) = column(INFO_MATRIX, INFO_MATRIX.size2()-1);

				cout << endl << "Pierwsza macierz simplex" << SIMPLEX_TABLE;

				//1stworzenie cb
				matrix<double> CB(1, receiversNumber);
				for (unsigned j = 0; j < receiversNumber; ++ j) {
						CB (0, j) = 0;
				}

				//stworzenie zj
				matrix<double> ZJ;
				ZJ = prod(CB, SIMPLEX_TABLE);

				//stworzenie cj-zj
				std::vector < int > cjzj;
				for (unsigned i = 0; i < ZJ.size2()-1; ++ i) {
					cjzj.push_back( requestion[i] - ZJ(0, i) );
				}

				/*std::vector < int > zmienne_bazowe_x;
				for (unsigned i = 0; i < providersNumber+receiversNumber; ++ i) {
					zmienne_bazowe_x.push_back( i );
				}
				for( size_t i = 0; i < zmienne_bazowe_x.size(); i++ )
					printf( "%d, ", zmienne_bazowe_x[ i ] );*/

				std::vector < int > zmienne_bazowe_y;
				for (unsigned i = providersNumber; i < providersNumber+receiversNumber; ++ i) {
					zmienne_bazowe_y.push_back( i );
				}


				// stworzenie bi/p o odpowiedniej dlugosci
				std::vector < int > bip;
					for (unsigned i = 0; i < receiversNumber; ++ i) {
						bip.push_back( 0 );
					}






				matrix<double> BASEA(receiversNumber, receiversNumber);
				matrix<double> INVERTA (receiversNumber, receiversNumber);
				matrix<double> TEMP_PROD;

				//stworzenie A=bazowej
				for (unsigned i = 0; i < receiversNumber; ++ i) {
					column(BASEA, i) = column(SIMPLEX_TABLE, providersNumber+i);
				}

				for (unsigned iteration = 0; iteration < 2; ++ iteration) {	//	ZMIENNA
					cout << endl << "---------------" << iteration << " ITERACJA -----------------" << endl << endl;

					// Pozycja najwiekszego cjzj
					auto biggestCJZJ = std::max_element(std::begin(cjzj), std::end(cjzj));
					auto positionBiggestCJZJ = std::distance(std::begin(cjzj), biggestCJZJ);

					// obliczenie bi/p
					for (unsigned i = 0; i < receiversNumber; ++ i) {
						if (SIMPLEX_TABLE(i, positionBiggestCJZJ) != 0) {
							bip[ i ] = SIMPLEX_TABLE(i, SIMPLEX_TABLE.size2()-1) / SIMPLEX_TABLE(i, positionBiggestCJZJ);
						}
						else {
							bip[ i ] = 1000000;
						}
					}
					for( size_t i = 0; i < bip.size(); i++ )
						printf( "%d, ", bip[ i ] );

					//dostosowanie A bazowaj
					auto lowestBip = std::min_element(std::begin(bip), std::end(bip));
					auto positionLowestBip = std::distance(std::begin(bip), lowestBip);
					
					column(BASEA, positionLowestBip) = column(SIMPLEX_TABLE, positionBiggestCJZJ);	//	ZMIENNA

					//if(iteration != 0) {
					//	column(BASEA, BASEA.size2()-iteration) = column(SIMPLEX_TABLE, iteration-1);
					//}
				cout << endl << endl << "BASEA" << BASEA;

					//stworzenie A-bazowej odwrotnej
					matrix<double> INVERTA (receiversNumber, receiversNumber);
					InvertMatrix(BASEA, INVERTA);
				cout << endl << endl << "INVERTA" << INVERTA;

					SIMPLEX_TABLE = prod(INVERTA, SIMPLEX_TABLE);
					//Nadpisanie pierwszej tablicy simplexowej
					//column(SIMPLEX_TABLE, iteration) = column(IDENTITY_MATRIX, IDENTITY_MATRIX.size2()-1-iteration);
					//column(SIMPLEX_TABLE, SIMPLEX_TABLE.size2()-2-iteration) = column(INVERTA, INVERTA.size2()-1-iteration);
					//column(SIMPLEX_TABLE, SIMPLEX_TABLE.size2()-1) = column(TEMP_PROD, TEMP_PROD.size2()-1);
				cout << endl << endl << "Druga tablica simplex" << SIMPLEX_TABLE;

					// Dostosowanie CB
					CB(0, CB.size2()-1-iteration) = requestion[iteration];
				cout << endl << endl << "CB" << CB;

					ZJ = prod(CB, SIMPLEX_TABLE);
				cout << endl << endl << "ZJ" << ZJ;

					for (unsigned i = 0; i < cjzj.size(); ++ i) {
						cjzj[i] = ( requestion[i] - ZJ(0, i) );
					}

				for( size_t i = 0; i < cjzj.size(); i++ )
					printf( "%d, ", cjzj[ i ] );



				// podmienienie starwej wartosci w A bazowej
				column(BASEA, positionLowestBip) = column(SIMPLEX_TABLE, positionBiggestCJZJ);


				}














					cout << endl << endl;
			}
				break;
			case 2: {
				system("cls");
				matrix<double> m (3, 3);
				for (unsigned i = 0; i < m.size1 (); ++ i)
					for (unsigned j = 0; j < m.size2 (); ++ j)
						m (i, j) = 3 * i + j;
				cout << "Macierz m" << m << std::endl;

				m(0,0) = 1;

				cout << "Macierz m(0,0) = 1: " << m << std::endl;

				matrix<double> mm (3, 3);
				for (unsigned i = 0; i < mm.size1 (); ++ i)
					for (unsigned j = 0; j < mm.size2 (); ++ j)
						mm (i, j) = 3 * i + j;
				cout << "Macierz mm: " << mm << endl;


				matrix<double> C;
				C = prod(m,mm);

				cout << "Macierz C" << C << endl;

				matrix<double> z (3, 3);
				InvertMatrix(m, z);

				cout << "A=" << m << endl << "Z=" << z << endl;

				cout << endl << endl;
			}
				break;
			case 3: {
				system("cls");
				menu = false;
				cout << "Mozesz teraz wyjsc z programu naciskajac ENTER";
			}
				break;
			case 4: {
				system("cls");
				matrix<double> m (3, 3);
				for (unsigned i = 0; i < m.size1 (); ++ i)
					for (unsigned j = 0; j < m.size2 (); ++ j)
						m (i, j) = 3 * i + j;
				cout << "Macierz m" << m << endl;
				cout << "Ilosc kolumn m" << endl;

				matrix<double> mm (3, 3);
				for (unsigned i = 0; i < mm.size1 (); ++ i)
					for (unsigned j = 0; j < mm.size2 (); ++ j)
						mm (i, j) = 3 * i + j;
				cout << "Macierz mm" << mm << endl;
				cout << "Kolumna mm" << column(mm, 0) << endl;

				m.resize(m.size1(), m.size2()+1, true);
				cout << "Macierz m" << m << endl;
				cout << "Ilosc kolumn m" << m.size2() << endl;

				column(m, m.size2()-1) = column(mm, 0);
				cout << "Macierz m" << m << endl;

				cout << endl << endl;
			}
				break;
		}
	}

	cin.sync(); 
	cin.get();
	return 0;
}
