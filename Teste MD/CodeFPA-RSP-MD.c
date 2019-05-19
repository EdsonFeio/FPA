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
double func (int i, int o); //Implementa as funcoes
double funcrecu (int i, int o, int s); //Implementa uma função recursiva para calcular apolarização
double rev (int i, int o); //revisa o resultado no while
double revF (); //revisa o resultado no while


#define PROB 0.7
#define MAX_ITER 50000
#define lamb 1.5
#define D 2
#define POP 10*D
#define ERRO 1.e-5
#define nm 5

int k=0,s=0;
double atual[POP][D], anterior[POP][D], g[D], r_temp_atual[POP], r_temp_prox[POP];
double mx[nm][D];

int main()
{
    srand( (unsigned)time(NULL));
    int i, j, corrida;
    double tempo,r=0;
    clock_t Ticks[2];
	
    FILE *ar;

    ar=fopen("testeBR.txt","a+");

    Ticks[0] = clock();
    for (corrida = 0; corrida < nm; corrida++){
	    Ticks[2] = clock();
        prencx();
	    g_atual();
        k=0;

        while ((fabs(rev(g[0],g[1])-0.397887)>ERRO)&&(k < MAX_ITER)){ 

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
        s += 1;
        if(k < MAX_ITER){
            for (i = 0; i < D; i++){
                mx[s][i] = g[i];
            }
        }
        Ticks[3] = clock();
        tempo = (Ticks[3] - Ticks[2]) * 1000.0 / CLOCKS_PER_SEC;
	    fprintf(ar,"Corrida: %d, Interacao: %d, Minimos: ( %g , %g ), Funcao: %g\n",corrida,k,g[0],g[1],revF());
        fprintf(ar,"Tempo: %g\n",tempo);
        Sleep(1000);
    }
    Ticks[1] = clock();
    tempo = (Ticks[1] - Ticks[0]) * 1000.0 / CLOCKS_PER_SEC;
    fprintf(ar,"\nSucesso, Tempo: %g\n\n",tempo-3000);
    fclose(ar);
    
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
    double m_solution;

    for (i = 0; i < POP; i++){
        r_temp_atual[i] = func(i, 1);
    }

    m_solution = r_temp_atual[0];

    for (i = 1; i < POP; i++){
        if (r_temp_atual[i] <= m_solution){
            m_solution = r_temp_atual[i];
            ind = i;
        }
    }
    for (i = 0; i < D; i++){
		g[i] = anterior[ind][i];
	}
}

void r_prox_geracao ()
{
    int i;

    for (i = 0; i < POP; i++){
        r_temp_prox[i] = func(i, 2);
    }
}

void aval_solucao ()
{
    int i, j;

    for (i = 0; i < POP; i++){
		for (j = 0; j < D; j++){
			if (r_temp_prox[i] < r_temp_atual[i]){
				r_temp_atual[i] = r_temp_prox[i];
				anterior[i][j] = atual[i][j];
			}
		}
    }
}

double func (int i, int o)
{
    double resp;
    if (s < 1){
        resp = rev(i,o) + 1.e-2;
        return resp;
    }else if (s >= 1){
        resp = (rev(i,o) + 1.e-2)/funcrecu(i,o,s);
        return resp;
    }
}

double funcrecu (int i, int o, int s)
{
    double resp, t = 0;
    double x[D];
	int j;
	
	if (o == 1){
		for (j = 0; j < D; j++){
			x[j] = anterior[i][j];
		}
	}else if(o == 2){
		for (j = 0; j < D; j++){
			x[j] = atual[i][j];
		}
	}
	
	for (j = 0; j < D; j++){
        t = t + pow(x[j]-mx[s][j],2);
    }
	
    if (s == 1){
        resp = atan(sqrt(t));
        return resp;
    }else if(s >= 2){
        resp = atan(sqrt(t))*funcrecu(i,o,s-1);
        return resp;
    }
}

double rev (int i, int o)
{
	double x[D], resp, b = 5.1/(4*(M_PI*M_PI)), c = 5.0/M_PI, h = 1.0/(8*M_PI);
	int j;
	
	if (o == 1){
		for (j = 0; j < D; j++){
			x[j] = anterior[i][j];
		}
	}else if(o == 2){
		for (j = 0; j < D; j++){
			x[j] = atual[i][j];
		}
	}
	
    //Branin Problem x = (-M_PI , 12.275), (M_PI , 2.275), (9.42478, 2.475) e f(x) = 0.397887        
	resp = pow((x[1] - (b*(x[0]*x[0])) + (c*x[0]) - 6),2) + 10*(1-h)*cos(x[0]) + 10;      

    return resp;
}    

double revF ()
{
    double x[D], resp, b = 5.1/(4*(M_PI*M_PI)), c = 5.0/M_PI, h = 1.0/(8*M_PI);
	int j;
	
	for (j = 0; j < D; j++){
		x[j] = g[j];
	}
	
    //Branin Problem x = (-M_PI , 12.275), (M_PI , 2.275), (9.42478, 2.475) e f(x) = 0.397887        
	resp = pow((x[1] - (b*(x[0]*x[0])) + (c*x[0]) - 6),2) + 10*(1-h)*cos(x[0]) + 10;        

    return resp;
} 
