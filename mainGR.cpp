#include "funciones.h"
#include <fstream>
#include <time.h>
#include <vector>
#include <string>
#include <map>

using namespace std;


int main(int argc, char *argv[])
{
	// para ejecutar ->  ./a.out -i 100-300 -t 80 -th 0.75 -pi 10 (poblacion inicial greedy)
    srand(time(NULL));
    string s1=argv[1];
    string s2=argv[2];
    string s3=argv[3];
    string s4=argv[4];
	string s5=argv[5];
    string s6=argv[6];
	string s7=argv[7];
	string s8=argv[8];
    char str[20];
    strcpy(str,s4.c_str());
    double timeLimit = atof(str);
	TIMELIM = timeLimit;
    if(s1 != "-i") return 0;
    if(s2 !="100-300" && s2!="100-600"&&s2!="100-800"&&s2!="200-300"&&s2!="200-600"&&s2!="200-800") return 0;
    if(s3 != "-t") return 0;
    if(timeLimit < 0 && timeLimit > 100) return 0;
    if(s5 != "-th") return 0;
    if(s6 !="0.75" && s6 != "0.80" && s6 != "0.8" && s6!="0.85") return 0;
	if(s7 != "-pi") return 0;
	int pi = stoi(s8);
	if(pi > 15 ) return 0;


    string aux = s2.substr(0,3);
    int filas = stoi(aux);
    aux = s2.substr(4,3);
    int columnas = stoi(aux);
    float threshold = stof(s6);
	
	vector<string> texts;
	for (int i = 0; i < 1; i++){
		texts.push_back(s2);
	}
	texts.at(0) += "-001.txt"; //AHORA SOLO ESTAMOS USANDO LA INSTANCIA 1 DSP USAREMOS LAS OTRAS 
	
	vector<long double> promedios;
	vector<char> alfabeto{'A', 'C', 'G', 'T'};
	vector<vector<char>> vec1; 


	fstream file1(texts.at(0), fstream::in | fstream::out);
	
	fillMatrix(vec1, file1,filas,columnas);

	vector<char> v_j;	  // vector columna en el que se esta trabajando
	vector<char> s_g;	  // solucion greedy
	vector<vector<char>> sols_vec(pi);

	vector<int> scount_g; // cantidad de la letra elegida en una iteracion de columna VJ.
	 // vector de 100 pos cada una representa la fila i de cada instancia y muestra en
	// cuantas casillas difiere de la solucion
	vector<cromosoma> pob; //poblacion incial
	cromosoma crom_aux;
	vector<int> fitness;
	vector<int> rep_j; // vector con repeticiones tamaño 4
	//---------------------------------------------------------PARTE AG -----------------------------------------------------------------------
	
	
	//inicializacion variables de tiempo
	auto start = chrono::high_resolution_clock::now();
	auto end = chrono::high_resolution_clock::now();
	chrono::duration<double> diff = end - start;
	vector<cromosoma> mutRepVec;
	//Se inicia el algoritmo genetico con la población inicial ya creada y evaluada
	for ( int i = 0; i < pi ; i++){ // se genera la población inicial y se evalua el fitness respectivo.
		crom_aux = greedy(threshold,vec1,v_j,rep_j,alfabeto,s_g,scount_g,filas,columnas);
			
		pob.push_back(crom_aux);
	}


	vector<cromosoma> padres;
	vector<cromosoma> nuevos_hijos;
	vector<cromosoma> cruza;
	int sizeof_column = vec1.at(0).size();
	float lim1 = static_cast<float>(sizeof_column)*threshold;
	int limite = static_cast<int>(lim1);
	cout<<"limite: "<< limite<< endl;//threshold de float a int
	int limit_padres = 10;
	int i = 1;
	cout<<"generacion: 0 "<<endl;
	printPobFit(pob);
	printPobStats(pob); 
	while(diff.count() < timeLimit){ //se detiene hasta que se acabe el tiempo limite de ejecución 
		
		//torneo se escogen los dos mejores padres de un conjunto 
		cromosoma padre1;
		cromosoma padre2;
		padre1=(torneo(pob)); 
		padre2=(torneo(pob));
		//Se recombina con los padres escogidos mediante una tecnica de cruzamiento (uniform crossover)
		cruza = uniformCrossover(padre1,padre2); 		 
		mutacion(cruza); 
		vector<int>solRep1(vec1.size(),0);
		vector<int>solRep2(vec1.size(),0);
		setNewFitness(cruza.at(0),vec1,threshold,solRep1);
		setNewFitness(cruza.at(1),vec1,threshold,solRep2);
		padres.push_back(padre1);
		padres.push_back(padre2);
		nuevos_hijos.push_back(cruza.at(0));
		nuevos_hijos.push_back(cruza.at(1));
		local_search(cruza.at(0),vec1,limite);
		//cout<<"check local search 1"<<endl;
		local_search(cruza.at(1),vec1,limite);
		//cout<<"check local search 2"<<endl;

		cruza.clear();
		//se genera una nueva población y se evalua ()
		//10 padres x 2 hijos por par de padres =  10 hijos nuevos  + 15 poblacion inicial(incluye padres) = 25 total
		// -> steady state() pob nueva 15 -> escojo 5 hijos de los mejores (nueva gen) y 5 padres de los mejores (gen vieja)
		if(padres.size() >= limit_padres){
			mutacion(pob);
			//localsearch
			cout<<endl;
			cout<<"generacion: "<<i<<" "<<endl;
			i++; //numero de iteración
			printPobFit(pob);
			printPobStats(pob);
			cout<<"tiempo actual:" <<diff.count() <<endl;
			steady_state(pob,5,nuevos_hijos,padres);
		}

		//ahora nuestra nueva población es nuestra población inicial de la siguiente iteración
		//en caso de terminar antes break;
		bool termino = verificar(pob,filas,threshold); 	
		if(termino) {
			cout<<endl;
			cout<<"generacion: "<<i<<" "<<endl;
			i++; //numero de iteración
			printPobFit(pob);
			printPobStats(pob);
			end = chrono::high_resolution_clock::now();
			diff = end - start;
			break; //si todos cumplen el fitness se termina
		}
		end = chrono::high_resolution_clock::now();
		diff = end - start;

	}
	cout<<"tiempo termino:" <<diff.count() <<endl;
	return 0;
}
