#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double randP(); //Gera um valor aleatorio entre 0 e 10,000...
double rande(); //Gera o epsylon aleatorio
void prencx(); //Preenche os vetores com valores entre 10,000... e -10,000...
double regrasimp(); //Calcula a intregral pela regra de simpson
double funcgama(double x); //Calcula a funcao gama
double v_Levy (int k); //Calcula a equacao do voo de levy
void local (int i); //Polinizacao local
void global (int i); //Polinizacao global
void g_atual (); //Analiza o g*(os melhores valores pop anterior)
void r_prox_geracao (); //Avalia quem vai ser a proxima geracao
void aval_solucao (); //Avalia a melhor solucao
double func (double x, double y); //Implementa as funcoes
double funcrecu (double x, double y, int s); //Implementa uma função recursiva para calcular apolarização
double rev (double x, double y); //revisa o resultado no while

#define PROB 0.7
#define MAX_ITER 10000
#define lamb 1.5
#define D 2
#define POP 10*D
#define ERRO 1.e-5
#define nm 9
#define cax 10

int k=0,ss=0,s=0,tt=0;
double atual[POP][D], anterior[POP][D], g[D], r_temp_atual[POP], r_temp_prox[POP];
double mx[nm],my[nm];

int main()
{
	srand( (unsigned)time(NULL));
	int i, j, corrida, min;
	double tempo, r=0;
	clock_t Ticks[2];
	
	FILE *ar;
	

	//ar=fopen("testeB2.txt","a+");
	//ar=fopen("testeCB3.txt","a+");
	//ar=fopen("testeE.txt","a+");
	//ar=fopen("testeGP.txt","a+");
	//ar=fopen("testeMRP.txt","a+");
	//ar=fopen("testePRD.txt","a+");
	//ar=fopen("testeSF2.txt","a+");
    
    //ar=fopen("testeBL.txt","a+");
    //ar=fopen("testeBR.txt","a+");
	ar=fopen("testeH.txt","a+");

	for (corrida = 0; corrida < nm; corrida++){
		Ticks[0] = clock();
		prencx();
		g_atual();
		k=0;
		
		//while ((fabs(rev(g[0],g[1]))>ERRO)&&(k < MAX_ITER)){
		//while ((fabs(rev(g[0],g[1]))>ERRO)&&(k < MAX_ITER)){ 
		//while ((fabs(rev(g[0],g[1])+1.0)>ERRO)&&(k < MAX_ITER)){ 
		//while ((fabs(rev(g[0],g[1])-3.0)>ERRO)&&(k < MAX_ITER)){
		//while ((fabs(rev(g[0],g[1]))>ERRO)&&(k < MAX_ITER)){
		//while ((fabs(rev(g[0],g[1])-0.9)>ERRO)&&(k < MAX_ITER)){
		//while ((fabs(rev(g[0],g[1]))>ERRO)&&(k < MAX_ITER)){
            
        //while ((fabs(rev(g[0],g[1]))>ERRO)&&(k < MAX_ITER)){
        //while ((fabs(rev(g[0],g[1])-0.397887)>ERRO)&&(k < MAX_ITER)){ 
		while ((fabs(rev(g[0],g[1]))>ERRO)&&(k < MAX_ITER)){

			for (i = 0; i < POP; i++){
				r = rand()%10;
				r /= 10;
				if ( r < PROB){
					global(i);
				} else {
					local(i);
				}
				r_prox_geracao();
				aval_solucao();
			}
			g_atual();
			k++;     	    
		}
		ss += 1;
		if(k < MAX_ITER){
			mx[ss] = g[0];
			my[ss] = g[1];
		}
		Ticks[1] = clock();
		tempo = (Ticks[1] - Ticks[0]) * 1000.0 / CLOCKS_PER_SEC;
		
		tt = tt + tempo;
		
		fprintf(ar,"Corrida: %d, Interacao: %d, Minimos: ( %g , %g ), Funcao: %g\n",corrida,k,g[0],g[1],rev(g[0],g[1]));
		fprintf(ar,"Tempo: %g\n",tempo);
	}
	tempo = tt;
	fprintf(ar,"\nSucesso, Tempo: %g\n\n",tempo);
	fclose(ar);
    return 0;
}

double randP(){
    double x=1+rand()%(cax*1000);
    x/=11;
    x/=11;
    x/=11;
    return x;
}

double rande(){
    double x=rand()%1000;
    x/=11;
    x/=11;
    x/=11;
    return x;
}

void prencx(){
    int i,j;
    for(i=0;i<POP;i++){
        for(j=0;j<D;j++){
           anterior[i][j]=(randP()-randP());
            atual[i][j]=anterior[i][j];
        }
    }
}

double regrasimp(){
    double x0=0,x2=1000;
	double h,x1,s;

	h = (x2 - x0)/(double)2;

	x1 = (x0 + x2)/(double)2;

	s=(h/3)*(funcgama(x0)+4*(funcgama(x1))+funcgama(x2));
	
	return s;
}


double funcgama(double x){

    double f;
    f=pow(x,lamb-1)*exp(-x);
	return f;
}

double v_Levy (int k)
{
    double L, s0, s, in;
    s0 = (double) k;
    s = 1+1e-5*s0;
    in = regrasimp();
    L = ((lamb * in * sin((M_PI * lamb) / 2))) / (M_PI*(pow(s, 1+lamb)));
    return L;
}

void local (int i)
{
    int j;
    for (j = 0; j < D; j++){
        double eps=rande();
        atual[i][j] =anterior[i][j] + (eps * (anterior[rand() % POP][j] - anterior[rand() % POP][j]));
    }
}

void global (int i)
{
    int j;
    for (j = 0; j < D; j++){
        atual[i][j] =anterior[i][j] + (v_Levy(k) * (anterior[i][j] - g[j]));
    }
}

void g_atual ()
{
    int i, ind=0;
    double x[POP], y[POP], m_solution;

    for (i = 0; i < POP; i++){
        x[i] =anterior[i][0];
        y[i] =anterior[i][1];
    }

    for (i = 0; i < POP; i++){
        r_temp_atual[i] = func(x[i], y[i]);
    }

    m_solution = r_temp_atual[0];

    for (i = 1; i < POP; i++){
        if (r_temp_atual[i] <= m_solution){
            m_solution = r_temp_atual[i];
            ind = i;
        }
    }
    
    g[0] = x[ind]; g[1] = y[ind];

}

void r_prox_geracao ()
{
    int i;
    double x[POP], y[POP];

    for (i = 0; i < POP; i++){
        x[i] = atual[i][0];
        y[i] = atual[i][1];
    }

    for (i = 0; i < POP; i++){
        r_temp_prox[i] = func(x[i], y[i]);
    }
}

void aval_solucao ()
{
    int i;


    for (i = 0; i < POP; i++){
        if (r_temp_prox[i] < r_temp_atual[i]){
            r_temp_atual[i] = r_temp_prox[i];
            anterior[i][0] = atual[i][0];
            anterior[i][1] = atual[i][1];
        }
    }
}

double func (double x, double y)
{
    double resp;
    s = ss;
    if (s < 1){
        resp = rev(x,y) + 1.e-8;
        return resp;
    }else if (s >= 1){
        resp = (rev(x,y) + 1.e-8)/funcrecu(x,y,s);
        return resp;
    }
}

double funcrecu (double x, double y, int s)
{
    double resp;
    if (s == 1){
        resp = atan(sqrt(pow(x-mx[s],2)+pow(y-my[s],2)));
        return resp;
    }else if(s >= 2){
        resp = atan(sqrt(pow(x-mx[s],2)+pow(y-my[s],2)))*funcrecu(x,y,s-1);
        return resp;
    }
}

double rev (double x, double y)
{
    double resp, b = 5.1/(4*(M_PI*M_PI)), c = 5.0/M_PI, h = 1.0/(8*M_PI);						 
	
	//Bohachevsky 2 Problem x = (0,0) e f(x) = 0 [-50,50]
	//resp = pow(x,2) + 2*(pow(y,2)) - 0.3*(cos(3*M_PI*x))*(cos(4*M_PI*y)) + 0.3;	    
	
    //Camel Back–3 Three Hump Problem x = (0,0) e f(x) = 0 [-5,5]
	//resp = 2*(pow(x,2)) - 1.05*(pow(x,4)) + (1.0/6)*(pow(x,6)) + x*y + pow(y,2);     		
	
    //Easom Function x = (M_PI, M_PI) e f(x) = -1 [-10,10]
	//resp = -cos(x)*cos(y)*exp(-((x-M_PI)*(x-M_PI)) - (y-M_PI)*(y-M_PI));      													

	//Goldstein and Price x = (0,-1) e f(x) = 3 [-2,2]
	//resp = (1+(pow((x+y+1),2))*(19-14*x+3*(pow(x,2))-14*y+6*x*y+3*(pow(y,2))))*(30+(pow((2*x-3*y),2))*(18-32*x+12*(pow(x,2))+48*y-36*x*y+27*(pow(y,2))));	

	//Modified Rosenbrock Problem x = (0.3412,0.1164), (1,1) e f(x) = 0 [-5,5]
    //resp = 100*(pow(y-(pow(x,2)),2)) + (pow(6.4*(pow(y-0.5,2))-x-0.6,2));                     						

	//Periodic Problem x = (0,0) e f(x) = 0.9 [-10,10]
    //resp = 1 + (pow(sin(x),2)) + (pow(sin(y),2)) - 0.1*(exp(-(pow(x,2))-(pow(y,2))));   

	//Schaffer 2 Problem x = (0,0) e f(x) = 0 [-100,100]
    //resp = (pow((pow(x,2))+(pow(y,2)),0.25))*((pow(sin(50*(pow((pow(x,2))+(pow(y,2)),0.1))),2))+1);  
	 
      
	
	//Becker and Lago Problem x = (-5,5), (5,-5), (-5,-5), (5,5) e f(x) = 0 [-10,10]
	//resp = pow((fabs(x) - 5),2) + pow((fabs(y) - 5),2);
	
	//Branin Problem x = (-M_PI , 12.275), (M_PI , 2.275), (9.42478, 2.475) e f(x) = 0.397887 [-10,10]  
	//resp = pow((y - (b*(x*x)) + (c*x) - 6),2) + 10*(1-h)*cos(x) + 10;  
	
	//Himmelblau problem 9 raizes f(x) = 0 [-5,5] 
	resp = pow(4*pow(x,3)+4*x*y+2*pow(y,2)-42*x-14,2)+pow(4*pow(y,3)+2*pow(x,2)+4*x*y-26*y-22,2);
	
    return resp;
}    
