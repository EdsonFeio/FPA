#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>

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

#define PROB 0.7
#define MAX_ITER 10000
#define lamb 1.5
#define D 2
#define POP 10*D
#define ERRO 1.e-3

int k,ss;
double atual[POP][D], anterior[POP][D], g[D], r_temp_atual[POP], r_temp_prox[POP];

int main()
{
    srand( (unsigned)time(NULL));
    int i, j, corrida;
    double tempo,r=0;
    clock_t Ticks[2];
    
    char d;
    FILE *ar;
    do{
        //ar=fopen("testeR.txt","a+");
        //ar=fopen("testeE.txt","a+");
        ar=fopen("testeB.txt","a+");
        
        Ticks[0] = clock();
        for (corrida = 0; corrida < 100; corrida++){
    	    Ticks[2] = clock();
            prencx();
    	    g_atual();
            k=0;
            //while ((fabs(func(g[0],g[1]))>ERRO)&&(k < MAX_ITER)){
            //while ((fabs(func(g[0],g[1])+1.0)>ERRO)&&(k < MAX_ITER)){
            while ((fabs(func(g[0],g[1])-0.397887)>ERRO)&&(k < MAX_ITER)){
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
            if (k < MAX_ITER) ss++;
            Ticks[3] = clock();
            tempo = (Ticks[3] - Ticks[2]) * 1000.0 / CLOCKS_PER_SEC;
    	    fprintf(ar,"Corrida: %d, Interacao: %d, Minimos: ( %g , %g ), Funcao: %g\n",corrida,k,g[0],g[1],func(g[0],g[1]));
            fprintf(ar,"Tempo: %g\n",tempo);
            Sleep(1000);
        }
        Ticks[1] = clock();
        tempo = (Ticks[1] - Ticks[0]) * 1000.0 / CLOCKS_PER_SEC;
        fprintf(ar,"\nSucesso: %d%%, Tempo: %g\n",ss,tempo-100000);
        fclose(ar);
        printf("\acontinuar? s ou n?");
        fflush(stdin); d=getchar(); system("cls");
    }while(d=='s'||d=='S');
    
    return 0;
}

double randP(){
    double x=1+rand()%10000;
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
    double resp, b = 5.1/(4*(M_PI*M_PI)), c = 5.0/M_PI, h = 1.0/(8*M_PI), alf = pow(10,-2);
     
    //resp = pow((x-1), 2) + (100 * pow(y - pow(x, 2), 2));                     //Rosenbrock Function x = (1, …, 1) e f(x) = 0
    //resp = -cos(x)*cos(y)*exp(-((x-M_PI)*(x-M_PI)) - (y-M_PI)*(y-M_PI));      //Easom Function x = (M_PI, M_PI) e f(x) = -1
    resp = pow((y - (b*(x*x)) + (c*x) - 6),2) + 10*(1-h)*cos(x) + 10;         //Branin Function x = (-M_PI , 12.275), (M_PI , 2.275), (9.42478, 2.475) e f(x) = 0.397887                
    
    return resp + alf;
}
