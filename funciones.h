#include"estructuras.h"

using namespace std;
double TIMELIM;
double CURRENTTIME=0;
void fillMatrix(vector<vector<char>> &vec, fstream &file,int fila,int columna)
{
	vector<char> aux;
	int flag = 0;
	string trash;
	for (int i = 0; i < fila; ++i)
	{
		vec.push_back(aux);
		for (int j = 0; j < columna+2; ++j)
		{
			if (flag < columna)
			{
				vec.at(i).push_back(file.get());
				flag++;
			}
			else
			{
				trash = file.get();
				flag++;
				if (flag == columna+2)
					flag = 0;
			}
		}
	}
}
void printMatrix(vector<vector<char>> &vec)
{
	for (int i = 0; i < vec.size(); ++i)
	{
		for (int j = 0; j < vec.at(i).size(); ++j)
		{
			cout << vec.at(i).at(j);
		}
		cout << endl;
		// cout<<vec.at(i).size();
	}
}

void record_pos_dif(vector<char>& inst, int pos, vector<char>& solucion, vector<int>& char_dif_x_fila)
{
	// cout<<" "<<inst.size()<<" "<< solucion.size()<<" ";
	for (int i = 0; i < inst.size(); i++){
		if (inst.at(i) != solucion.at(pos)) char_dif_x_fila.at(i) += 1;
	}
}

int notificar(vector<int> &sol_posRecord, int limite,int columnas, int filas){ // sol_posrecord size es 100
	int total = 0; //cuantas filas cumplen con el threshold
	// cout<< "compare_th: "<<compare_th<<endl;
	for (int i = 0; i < sol_posRecord.size(); i++){
		//cout<<"FILA: "<<i<< " valor: "<<sol_posRecord.at(i)<<" ";
		if (sol_posRecord[i] >= limite ) total += 1;
	}
	cout<<total<<" de "<< filas <<"strings cumplen con el threshold"<<endl;
	return total;
}

void firstchoice(vector<int> &rep, vector<char> &sol, vector<int> &n_char_sol){//repeticiones por letra de alfabeto en una columna, letra escogida , cantidad de esa escogida
	int save = 0;
	for (int i = 0; i < rep.size() - 1; i++){
		if (rep.at(i) > rep.at(i + 1))
			save = i + 1;
	}
	n_char_sol.push_back(rep.at(save));
	switch (save)
	{
	case 0:
		sol.push_back('A');
		break;
	case 1:
		sol.push_back('C');
		break;
	case 2:
		sol.push_back('G');
		break;
	case 3:
		sol.push_back('T');
		break;
	}

}



//resuelve un archivo de texto. una instancia.
cromosoma greedy (float threshold,vector<vector<char>>& mat, vector<char>& vj,vector<int>& reps_j,vector<char>& alf,vector<char>& sol,vector <int>& n_char_sol,int filas,int columnas){
    
	random_device random; // genera un dado que se tira aleatoriamente
	mt19937 engine{random()}; // se planta una semilla https://cpppatterns.com/patterns/choose-random-element.html 
	uniform_real_distribution<> prob(0,0.1);// que genera una distribucion real desde 0 a 0.1
	
	float e = 0.1;//
	sol.clear(); //se limpia el vector solucion
    n_char_sol.clear(); 
	vector<int> solRep(filas, 0);// cuantos chars difieren de la solucion en una fila 
	float lim1 = static_cast<float>(mat.at(0).size())*threshold;
	int lim = static_cast<int>(lim1); //threshold traducido a cantidad de caracteres 
	map<char,int> repeticiones;// mapa para saber cuanto se repiten los caracteres del alfabeto, se actualiza x columna
    for (int j = 0; j < lim; j++){ // x columna
		vj.clear();
		reps_j.clear();
		repeticiones['A'] = 0;
		repeticiones['C'] = 0;
		repeticiones['G'] = 0;
		repeticiones['T'] = 0;
		for (int i = 0; i < mat.size(); i++){//
			vj.push_back(mat.at(i).at(j)); //lleno la columna actual 
		}
		for(int k = 0; k < vj.size();k++){
				repeticiones[vj.at(k)]++;
		}
		reps_j.push_back(repeticiones['A']);
		reps_j.push_back(repeticiones['C']);
		reps_j.push_back(repeticiones['G']);
		reps_j.push_back(repeticiones['T']);
		
		firstchoice(reps_j,sol,n_char_sol); //se mete a sol el caracterer menos repetido.	
		record_pos_dif(vj,j,sol,solRep); //por columna se guarda en solRep cuantas filas difieren del caracter escogido 
	}
	//hasta aqui se arma una sol con las menos repetidas en sol hasta el threshold

	for (int j = lim; j < mat.at(0).size(); j++){ //desde donde se construyo se continua 
		repeticiones['A'] = 0; //reciclamos el mapa repeticiones a un mapa de puntajes por 
		repeticiones['C'] = 0;// letra a medida de que se revise una columna
		repeticiones['G'] = 0;
		repeticiones['T'] = 0;
		for (int i = 0; i < filas; i++){
			for (map<char,int>::iterator it= repeticiones.begin();  it != repeticiones.end(); it++){//taamaño del alfabeto
				int aux = 0;
				if (it->first != mat.at(i).at(j)){
					it->second++;
					aux = solRep[i];
					aux++;
				}
				if(aux >= lim) it->second++;
			}//termina for de 4 pasadas	
		}//termino for de filas
		//------------------------------- se genera gen random ------------------------------------
		float p = prob(engine);
		uniform_int_distribution<int> dist(0,alf.size()-1); //se usa otra funcion con el engine para escoger una posicion aleatoria
		
		if(p<e){
			sol.push_back(alf[dist(engine)]); // se llama dist(engine) para agregar la posicion aleatoria escogida
		}else{
			auto pr = max_element(repeticiones.begin(),repeticiones.end(),[](const auto &x, const auto &y){
				return x.second < y.second; // funcion para comṕarar valores del mapa
			});
			//cout<<"letra escogida para iteracion :"<< j <<" "<<pr->first<<" puntaje: "<<pr->second<<endl;
			sol.push_back(pr->first);
		}
	}
	int fitness = notificar(solRep,lim,columnas,filas); //calcula fitness 
	cromosoma aux = {sol,fitness,solRep};
    //notificar(solRep,lim,columnas,filas);
	return aux;

}

void printPobFit(vector<cromosoma>& pob){
    for (int i = 0; i < pob.size(); ++i){
        cout<<"fitness "<<i<<" : "<<pob.at(i).fitness<<endl;
    }
}
void printPobStats(vector<cromosoma>& pob){
    long double prom=0;
    cout.precision(9);
    for (int i = 0; i < pob.size(); i++){
        prom=prom + pob.at(i).fitness;
    }
    prom = prom / pob.size();
    cout<<"fitness promedio:  "<<prom<<endl;
    long double sumatoria=0;
    long double varianza=0;
    for (int i = 0; i < pob.size(); i++){
        sumatoria=sumatoria + pow(pob.at(i).fitness-prom,2);
    }
    varianza = sumatoria/(pob.size());
    cout<<"varianza:  "<<varianza<<endl;
    long double desviacion=sqrt(varianza);
    cout<<"Desviacion "<<desviacion<<endl;
    return;
}

cromosoma greedyfalso (float threshold, vector<vector<char>>& mat, vector<char>& vj,vector<int>& reps_j,vector<char>& alf,vector<char>& sol,vector <int>& n_char_sol,int filas,int columnas){
    
	random_device random; // genera un dado que se tira aleatoriamente
	mt19937 engine{random()}; // se planta una semilla https://cpppatterns.com/patterns/choose-random-element.html 
	uniform_real_distribution<> prob(0,0.1);// que genera una distribucion real desde 0 a 0.1
	
	float e = 0.1;//
	//sol.clear(); //se limpia el vector solucion
    //n_char_sol.clear(); 
    float lim1 = static_cast<float>(mat.at(0).size())*threshold;
	int lim = static_cast<int>(lim1);
	vector<int> solRep(filas, 0);// cuantos chars difieren de la solucion en una fila 
	map<char,int> repeticiones;// mapa para saber cuanto se repiten los caracteres del alfabeto, se actualiza x columna
    for (int j = 0; j < mat[0].size(); j++){ // x columna
		vj.clear();
		reps_j.clear();
		repeticiones['A'] = 0;
		repeticiones['C'] = 0;
		repeticiones['G'] = 0;
		repeticiones['T'] = 0;
		for (int i = 0; i < mat.size(); i++){//
			vj.push_back(mat.at(i).at(j)); //lleno la columna actual 
		}
		for(int k = 0; k < vj.size();k++){
				repeticiones[vj.at(k)]++;
		}
		reps_j.push_back(repeticiones['A']);
		reps_j.push_back(repeticiones['C']);
		reps_j.push_back(repeticiones['G']);
		reps_j.push_back(repeticiones['T']);
		
		firstchoice(reps_j,sol,n_char_sol); //se mete a sol el caracterer menos repetido.	
		record_pos_dif(vj,j,sol,solRep); //por columna se guarda en solRep cuantas filas difieren del caracter escogido 
	}
    int fitness = notificar(solRep,lim,columnas,filas);
    cromosoma aux = {sol,fitness,solRep};
    
	return aux;

}

void mutacion2(vector<cromosoma>& vec,cromosoma mutvec){
  
    int h = rand()%vec.size();
    for (int i = 0; i < vec[h].gen.size(); i++){
        if( (rand()%1) == 1 ){
            vec[h].gen[i] = mutvec.gen.at(i);
        }
    }    
}
void mutacion(vector<cromosoma>& vec){
    vector<char> alf{'A', 'C', 'G', 'T'};
    random_device random;
    mt19937 engine{random()};
    uniform_int_distribution <int> alfa_pos(0,3);
    int h = rand()%vec.size();
    for (int i = 0; i < vec[h].gen.size(); i++){
        if( (rand()%1) == 1 ){
            int pos_aux = alfa_pos(engine);
            while(vec[h].gen[i] == alf[pos_aux]){ //asegurarse que no sea el mismo char que el gen del hijo
                int pos_aux = alfa_pos(engine);
            }
            vec[h].gen[i] = alf[pos_aux];
        }
    }    
}

bool verificar(vector<cromosoma>& pob,int filas,float threshold){ //verificar termino.
    float lim = (float)filas*threshold;
    long double prom=0;
	cout.precision(9);
	for (int i = 0; i < pob.size(); i++){
		prom = prom + pob.at(i).fitness;
	}
	prom = prom / pob.size();
	long double sumatoria = 0;
	long double varianza = 0;
	for (int i = 0; i < pob.size(); i++){
		sumatoria=sumatoria + pow(pob.at(i).fitness-prom,2);
	}
	varianza = sumatoria/(pob.size());
	long double desviacion=sqrt(varianza);
    if(prom > lim) return true;
    else return false;
}


char random_gen(){ // "GEN" random para la cruza de padres
    vector<char> alf{'A', 'C', 'G', 'T'};
    random_device random;
    mt19937 engine{random()};
    uniform_int_distribution <int> alfa_pos(0,3); // rand pos 0 a 3 tam alf
    int pos_gen = alfa_pos(engine); // posicion random
    char gen = alf.at(pos_gen); // sacar el char de la posicion aleatoria de vector alfabeto
    return gen; // retorna char aleatorio
}

void setNewFitness(cromosoma& crom, vector<vector<char>>& mat ,float threshold, vector<int> solRep_aux){ 
    int fitness = 0;
    int filas = mat.size();
    int columnas = crom.gen.size();
    float lim = (float)columnas*threshold; 
    float aux = 0; //contador de posiciones diferentes

    for (int i = 0; i < filas; i++){
        for(int j = 0; j < columnas; j++){
            if(crom.gen.at(j) != mat.at(i).at(j)) aux++;
            if(aux >=  lim){
                fitness++;
                aux = 0;
                break;
            }
        }
        aux = 0;
    }
    for (int j = 0; j < columnas; j++){
        for (int i = 0; i < filas; i++){
            if(mat[i][j] != crom.gen[j]) solRep_aux[i]++;
        }
        
    }
    crom.solRep = solRep_aux;
    crom.fitness = fitness;
    return;
}

vector<cromosoma> uniformCrossover(cromosoma& padre1,cromosoma& padre2 ){  
    vector<cromosoma> cruza;
    cromosoma hijo1;
    cromosoma hijo2;
    int cromSize = padre1.gen.size();
    for(int i = 0; i<cromSize ;i++){
        int aux = rand() % 2; //da 0 o 1
        if(aux){
            hijo1.gen.push_back(padre1.gen.at(i));//padre1 le hereda su gen al hijo1
            hijo2.gen.push_back(random_gen());//hijo 2 genera un GEN RANDOM
        }else{
            hijo2.gen.push_back(padre2.gen.at(i));//padre2 le hereda su gen al hijo2
            hijo1.gen.push_back(random_gen());//hijo 1 genera un GEN RANDOM
        } 
    }
    hijo1.fitness = 0;
    hijo2.fitness = 0;
    cruza.push_back(hijo1);
    cruza.push_back(hijo2);
    return cruza; //devolver los dos nuevos genes cruzados
}

cromosoma torneo(vector <cromosoma>& poblacion){
    int s=2; //nro de seleccion de participantes de torneo, comunmente se usa 2
    random_device random;//ej: uniform_real_distribution<> prob(0,0.1); uniform_int_distribution<int> dist(0,alf.size()-1); 
    mt19937 engine{random()};
    uniform_int_distribution<int> pos(0,poblacion.size()-1); //se escogen dos cromosomas random para que compitan
    int pos1 = pos(engine);
    int pos2 = pos(engine);

    cromosoma ganador;
    
    while (pos1 == pos2) pos2 = pos(engine); //para evitar torneo entre iguales.
    
    cromosoma padre1 = poblacion.at(pos1);
    cromosoma padre2 = poblacion.at(pos2);
   
    if(padre1.fitness >= padre2.fitness) ganador = padre1;
    else ganador = padre2;

    return ganador;
}

void steady_state(vector<cromosoma>& poblacion_inicial,int elitismo, vector<cromosoma>& nuevos,vector<cromosoma>& padres){ 
        vector<cromosoma> new_pob; // antigua poblacion tomada 5
        vector<int> fitness_h; //posiciones tomadas
        //vector nuevos son los hijos nuevos
        cromosoma best;
        int best_fit;
        
        for (int i = 0; i < poblacion_inicial.size()-1; i++){
            for (int j = 0; j < poblacion_inicial.size() - i-1; j++){
                if(poblacion_inicial[j].fitness < poblacion_inicial[j + 1].fitness){
                    cromosoma aux = poblacion_inicial[j];
                    poblacion_inicial[j] = poblacion_inicial[j + 1];
                    poblacion_inicial[j + 1] = aux;
                }
            }
        }

        for (int i = 0; i < elitismo; i++){  //primeros 5 con mayor fitness
            new_pob.push_back(poblacion_inicial[i]); 
        }

        for (int i = 0; i < nuevos.size()-1; i++){   
            for (int j = 0; j < nuevos.size() - i-1; j++){
                if(nuevos[j].fitness < nuevos[j + 1].fitness){
                    cromosoma aux = nuevos[j];
                    nuevos[j] = nuevos[j + 1];
                    nuevos[j + 1] = aux;
                }
            }
        }
        for (int i = 0; i < elitismo; i++){
            new_pob.push_back(nuevos[i]);
        }
        
        poblacion_inicial.clear();
        padres.clear();
        poblacion_inicial = new_pob;
        return;
}
void local_search(cromosoma &crom,vector<vector<char>> &mat,int limite){
	vector<char> best;// mejor primera solucion
	vector<char> alfabeto {'A','C','G','T'};
    //vector<char> sol_aux = crom.gen;
    int filas = mat.size();
    float comp = static_cast<float> (filas);
	//vector<char> vj;
	int sol_position = crom.gen.size()-1;
	bool flag = true;
    int intento = 3; // tres intentos
    int fitness_aux = crom.fitness;
	while(sol_position != 0){// comienzo del ciclo
        //cout<<"entre al loop inicial de local search mi posicion en rvision es: "<<sol_position <<endl;
        random_device random;//tienes que ajustar el ciclo a que vaya en incremento de a 1 hasta que llegue al final
        mt19937 engine{random()};// avanza 1, vuelve  hacioa atra pero con un nuevo limite por donde buscar
        uniform_int_distribution<int> pos_alf(0,3);
        
        
        while(intento != 0){ // dura maximo 3 veces
            //cout<<"entre el loop para probar distintos chars para cambiar intento: "<< intento <<endl;
            char cambio = alfabeto[pos_alf(engine)];
            //int media_principal = 0;
            //int media_siguiente = 0;
            int fitness_aux = crom.fitness;
            vector <int>numRep_aux = crom.solRep;
            if(cambio == crom.gen[sol_position])cambio = alfabeto[pos_alf(engine)];// evita el mismo char que el escogido original
		    for (int i = 0; i < filas; i++){
                //cout<<"numrep en posicion "<<i<<":"<< numRep_aux[i]<<endl;
                
                //media_siguiente += numRep_aux[i];
                //media_principal += crom.solRep[i];

                if(numRep_aux[i] >= limite) continue; // evitar que entre antes de aumentarlo si ya cumple con el fitness
                
                if(mat[i][sol_position] != cambio){
                    numRep_aux[i]++;
                    //media_siguiente++;
                }
                if(mat[i][sol_position] == cambio){
                    numRep_aux[i]--;
                    //media_siguiente--;
                }
                

                if(numRep_aux[i] >= limite){//solo si logra mejorar fitness y no era parte del fitness antes
                    //cout<<"entre con el char: "<< cambio <<endl;
                    fitness_aux++;
                    flag = false;
                    //break;
                }
		    }

            //float media1 = static_cast<float>(media_principal)/comp;
            //float media2 = static_cast<float>(media_siguiente)/comp;

            //cout<<media_principal<<" "<<media_siguiente<<endl;
            if(flag == false /*|| media2 > media1*/){
                crom.gen[sol_position] = cambio;
                crom.fitness = fitness_aux;
                crom.solRep = numRep_aux;
            }
           
            intento--;
        }
        intento = 3;
        if(flag == false){
            sol_position = crom.gen.size()-1;
            flag = true;
        }else{
            sol_position--;
        }
        
	}
}
