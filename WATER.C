#include<stdio.h>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include<malloc.h>

#define sigma 1.0
#define T 300.0
#define pasoDg 10
#define units 100
#define Rcut (2.5)*sigma
#define N 4
#define M_PI 3.1415
#define CONST ((9.035*9.035*100)/(16*M_PI*M_PI*6.02*6.02*1.38))

/*double const M_PI = 3.1415;*/

void init_celda( );
void monte_carlo( long int x );
void front_condt( int n );
int check_overlap( int n );
/*int suma_capas( );
int S_CAPAS( );*/
int paso_g,cuenta_g,Xunits;
long int actualiza_g,ACEPTANCY,COUNT;
double aleat( );
double set_energy(  );
double DELTA( int n );
double g_en_r( );
double X[ 1 ];
double ENERGY,L,A;
double XYZ[3*(N*N*N+1)];
double g[units];
double gbackup[units];
float B;

main(){
	long int pruebas,j,i,PASOS;
	int h;

	FILE *ENERGIA;
	FILE *RDF;

	printf( "CALCULO DE ENERGIA OPTIMIZADA\n\n" );
	
	if( XYZ != 0 ){
		printf( "se reservo memoria para XYZ\n" );
	}
	else{
		printf( "no se reservo memoria para XYZ\n" );
	}

	printf( "introduzca el numero de pasos monte-carlo" );
	scanf( "%ld" , &PASOS );

	/***************PASO DONDE COMIENZA EL CONTEO DE g(r)********************/

	printf( "introduzca el paso de actualizacion de g(r)" );
	scanf( "%ld" , &actualiza_g );

	/************************************************************************/


	/*************GROSOR DE LA CAPA DONDE SE CUENTAN LAS PARTICULAS**********/


	/************************************************************************/

	printf( "introduzca valor de CAMPO" );
	scanf( "%f" , &B );

	X[0] = 123457.0;

	ACEPTANCY = 0;
	COUNT = 0;

	L = 0.31*N/0.2791;
	A = L/N;

	/*cuenta_g = S_CAPAS();*/

	if(( g != 0 ))printf("se reservo memoria para g\n");
	else printf("no se reservo memoria para g") ;
	
	if(( gbackup != 0 ))printf("se reservo memoria para gbackup\n");
	else printf("no se reservo memoria para gbackup") ;

	init_celda( );

	/*for( h = 0; h < cuenta_g; h ++ )gbackup[h] = g[h] = 0.0;*/

	i = 0;

	ENERGY = set_energy(  );

	printf( "ENERGIA INICIAL = %1.10lf\n",(float)0.001*1.38*6.022*ENERGY/(N*N*N) );
	printf( "LONGITUD DE LA CELDA = %1.10lf\n",L);
	printf("cuenta_g=%d\n",cuenta_g);

	ENERGIA =  fopen("ENERGIO.txt","w");

	for( i = 0; i < PASOS; i ++)
	{
		monte_carlo( i );
		printf( "%ld\t %1.10lf\n",i,0.001*1.38*6.022*ENERGY/(N*N*N) );
		  fprintf(ENERGIA,"%ld\t%1.10lf\n",i,(0.001*(1.38*6.022*(ENERGY))/(N*N*N)));

	}
	fclose( ENERGIA );

	RDF = fopen("RDF_B000.txt","w");
	for( h = 0; h < units; h ++)
	{
		g[h] =   g[h]*pow(L,3)/((PASOS-actualiza_g)*4*pow(N,3)*pow(N,3)*M_PI*
			    pow((h+1)*pow(pasoDg,-1) ,2)*pow(pasoDg,-1))   ;

		fprintf(RDF,"%1.2f\t%1.2lf\n",h*pow(pasoDg,-1),g[h]);
		/*   printf("%1.2f\t%1.2lf\n",sigma + (h)*pow(paso_g,-1),g[h]);*/

	}
	fclose( RDF );
	
	printf("FINALIZADO");

	system("PAUSE");
	return(0);
}


/*****ESTA FUNCION GENERA IONES DISTRIBUIDOS  DENTRO DE  LA CELDA ELEMENTAL***/

void init_celda( ) {
	int i,j,k,l;
	l = 0;

	for(i = 0; i < 3*(N*N*N+1); i ++)XYZ[ i ] = 0;

	for(i = 0; i < N; i ++)
	{
		for(j = 0; j < N; j ++)
		{
			for(k = 0; k < N; k ++,l ++)
			{

				XYZ[ 3*l     ] = 0.5*A + k*A + 0.5*(A-sigma)*(aleat()-0.5);
				XYZ[ 3*l + 1 ] = 0.5*A + j*A + 0.5*(A-sigma)*(aleat()-0.5);
				XYZ[ 3*l + 2 ] = 0.5*A + i*A + 0.5*(A-sigma)*(aleat()-0.5);


			}
		}
	}
}

/***********************GENERA LOS PASOS MONTE - CARLO***********************/


void monte_carlo( long int x )
{
	double DELTA_E,RADIO_INTERACCION;
	double paso_unidades_sigma;
	int ch_prb,i,j,n,s,ch_pass,pass,k,l,m,w,M[ 3 ];
	long int h;

	M[ 0 ] =  0;
	M[ 1 ] =  1;
	M[ 2 ] = -1;

	/*paso_unidades_sigma = 0.051;*/
	paso_unidades_sigma = 0.30;

	for(w = 0; w <N*N*N; w ++)
	{
		n = (int)( N*N*N*aleat() );

		for( i = 0; i < 3; i ++ )
		{
			XYZ[3*(N*N*N+1)+i] = XYZ[3*n + i];
			XYZ[3*n+i] += paso_unidades_sigma*(aleat()-0.5);
		}

		ch_pass = check_overlap( n );

		if(ch_pass == 0)
        {
		  for(i = 0; i < 3; i ++)XYZ[3*n + i ]=XYZ[3*(N*N*N+1)+i];
		}

		if(ch_pass == 0)continue;

		DELTA_E = DELTA( n );

		if( DELTA_E <= 0.0 )
        {
				ENERGY += DELTA_E;
		}
		else
		{
			ch_prb = ( exp(-0.003*DELTA_E )>aleat() );

			if( ch_prb == 1 )
            {
			  ENERGY += DELTA_E;
			}
			else
			{
				for(i=0; i<3  ;i++)XYZ[ 3*n + i ] = XYZ[3*(N*N*N+1)+i];
			}

		}
		       front_condt( n );
	}
	
	
	if( x >= actualiza_g )
    {
			for( s=0; s<units; s++ )gbackup[s] = 0;


			for(k = 0; k < 3; k ++)
			{
			for(m = 0; m < 3; m ++)
			{
			for(l = 0; l < 3; l ++)
			{

			  for(i = 0; i < N*N*N; i ++)
			  {
			  for(j = 0; j < N*N*N; j ++)
			  {
				if( (i==j)&&( (M[k]==0)&&((M[m]==0)&&(M[l]==0)) ) )continue;

								RADIO_INTERACCION =
								    sqrt( pow( (XYZ[  3*i      ] - (XYZ[  3*j      ] + M[ m ]*L)), 2) +
									    pow( (XYZ[ (3*i) + 1 ] - (XYZ[ (3*j) + 1 ] + M[ k ]*L)), 2) +
									    pow( (XYZ[ (3*i) + 2 ] - (XYZ[ (3*j) + 2 ] + M[ l ]*L)), 2) );


								if(RADIO_INTERACCION < pasoDg+1)
                                {
								 for( s=0; s<units; s++ )
								 {
								  gbackup[s] += ((RADIO_INTERACCION >= s*pow(pasoDg,-1) )&&
												    (RADIO_INTERACCION < (s+1)*pow(pasoDg,-1) ));
								 }
                                }


				}
				}


			  }
			  }
			  }

			for(i=0; i<units; i++)g[i]+=gbackup[i];
	 }


}


/**********ESTA FUNCION CALCULA EL NUMERO DE CAPAS A CONTAR PARA g(r)********/

/*int S_CAPAS(){
	int j = 0;
	double CAPAS;

	CAPAS = pow(paso_g,-1);

	do{

		CAPAS += pow(paso_g,-1);
		j ++;

	}
	while(CAPAS < Xunits);

	g =(double*)malloc( j*sizeof(double) );

	return( j );
}*/

/*******CUENTA LA CANTIDAD DE CAPAS A TENER EN CUENTA PARA CALCULAR g(r)******/

/*int suma_capas( ){
	int      j=0;
	do{

		g =(double*)malloc( sizeof(double) +j*sizeof(double) );

		j++;

	}
	while(j/paso_g <= 3);

	return( j );
}*/


/*******************GENERADOR DE NUMEROS PSEUDO - ALEATORIOS*******************/

double aleat(){
	unsigned k;
	double m;

	m=pow(2,31)-1;
	k=16807U;

	X[0]=fmod((k*X[0]),m);

	return(X[0]/m);

}


/******************************************************************************/



/* ESTA FUNCION CALCULA LA ENERGIA DE LA CELDA ELEMENTAL Y
SUS REPLICAS */

double set_energy( ){
	int i,j,m,l,k,s,n,M[ 3 ];
	double RADIO_INTERACCION,E;

	E = 0;

	M[ 0 ] =  0;
	M[ 1 ] =  1;
	M[ 2 ] = -1;



	/****************  ENERGIA  CON LAS CELDAS IMAGENES    *******************/

	for( s=0; s<cuenta_g; s++ )gbackup[s] = 0;


	for(k = 0; k < 3; k ++)
	{
		for(m = 0; m < 3; m ++)
		{
			for(l = 0; l < 3; l ++)
			{

				for(i = 0; i < N*N*N-1; i ++)
				{
					for(j = i+1; j < N*N*N; j ++)
					{
						/*if( (i==j)&&( (M[k]==0)&&((M[m]==0)&&(M[l]==0)) ) )continue;*/

						RADIO_INTERACCION =
						    sqrt( pow( (XYZ[  3*i      ] - (XYZ[  3*j      ] + M[ m ]*L)), 2) +
							    pow( (XYZ[ (3*i) + 1 ] - (XYZ[ (3*j) + 1 ] + M[ k ]*L)), 2) +
							    pow( (XYZ[ (3*i) + 2 ] - (XYZ[ (3*j) + 2 ] + M[ l ]*L)), 2) );



						if((RADIO_INTERACCION) <= Rcut)E +=
						    (4*423.4)*(pow(( sigma/ RADIO_INTERACCION),12)-
							    pow((sigma/ RADIO_INTERACCION) ,6));

						E += ( CONST*B*B/pow(RADIO_INTERACCION, 3) )*
						    (1 - 3*pow(XYZ[ 3*i ] - (XYZ[ 3*j  ]+M[ m ]*L) ,2)/pow(RADIO_INTERACCION,2));

					}
				}

			}
		}
	}


	return( E );
}



/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/



double DELTA(int n){
	int i,m,l,k,M[ 3 ];
	double RADIO_INTERACCION;
	double E1 = 0;
	double E2 = 0;

	M[ 0 ] =  0;
	M[ 1 ] =  1;
	M[ 2 ] = -1;

	/*******************ENERGIA ANTES DE MOVER LA PARTICULA*********************/

	for(k = 0; k < 3; k ++)
	{
		for(m = 0; m < 3; m ++)
		{
			for(l = 0; l < 3; l ++)
			{
				/*ENERGIA DE LAS PARTICULAS DE LA CELDA BASE
				CON LA ESCOGIDA Y SUS IMAGENES           */

				for( i = 0; i < N*N*N; i ++ )
				{
					if(  (i==n)&&( (M[k]==0)&&( (M[m]==0)&&(M[l]==0) ) )  )continue;

					RADIO_INTERACCION =
					    sqrt( pow( XYZ[ 3*i     ] - (XYZ[3*(N*N*N+1)  ] + M[ k ]*L), 2) +
						    pow( XYZ[ 3*i+1   ] - (XYZ[3*(N*N*N+1)+1] + M[ m ]*L), 2) +
						    pow( XYZ[ 3*i+2   ] - (XYZ[3*(N*N*N+1)+2] + M[ l ]*L), 2) );


					if(RADIO_INTERACCION <= Rcut)E1 +=
					    (4*423.4)*(pow(( sigma/ RADIO_INTERACCION),12)-
						    pow(( sigma/ RADIO_INTERACCION) ,6));

					E1 += ( CONST*B*B/pow(RADIO_INTERACCION, 3) )*
					    (1 - 3*pow(XYZ[ 3*i ] - (XYZ[3*(N*N*N+1)  ]+M[ k ]*L) ,2)/	        pow(RADIO_INTERACCION,2));


				}
			}
		}
	}


	/*****************ENERGIA DESPUES DE MOVER LA PARTICULA*********************/


	for(k=0; k<3; k++)
	{
		for(m=0; m<3; m++)
		{
			for(l=0; l<3; l++)
			{

				for(i=0; i<N*N*N; i++)
				{
					if(  (i==n)&&( (M[k]==0)&&( (M[m]==0)&&(M[l]==0) ) )  )continue;


					/*ENERGIA DE LAS PARTICULAS DE LA CELDA
					BASE  CON LA MOVIDA Y SUS IMAGENES  */

					RADIO_INTERACCION =
					    sqrt( pow( XYZ[ 3*i     ] - (XYZ[ 3*n     ] + M[ k ]*L), 2)  +
						    pow( XYZ[ 3*i + 1 ] - (XYZ[ 3*n + 1 ] + M[ m ]*L), 2)  +
						    pow( XYZ[ 3*i + 2 ] - (XYZ[ 3*n + 2 ] + M[ l ]*L), 2)  );

					if(RADIO_INTERACCION <= Rcut)E2 +=
					    (4*423.4)*(pow(( sigma/ RADIO_INTERACCION),12)-
						    pow(( sigma/ RADIO_INTERACCION) ,6));

					E2 += ( CONST*B*B/pow(RADIO_INTERACCION, 3) )*
					    (1 - 3*pow(XYZ[ 3*i ] - (XYZ[ 3*n ]+M[ k ]*L) ,2)/	       pow(RADIO_INTERACCION,2));

				}
			}
		}
	}







	return((E2 - E1));
}



/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/*ESTA FUNCION REALIZA EL CONTROL DE FRONTERAS*/

void front_condt( int n ){
	int i;

	for( i = 0; i < 3; i++ )
	{
		if( XYZ[ 3*n + i ] > L )XYZ[ 3*n + i ] -= L;
		if( XYZ[ 3*n + i ] < 0 )XYZ[ 3*n + i ] += L;
	}

}



/******************************************************************************/



int check_overlap(int n){
	int j,k,l,m,pass,M[ 3 ];

	pass = 1;
	M[ 0 ] =  0;
	M[ 1 ] =  1;
	M[ 2 ] = -1;

	for( k = 0; k < 3; k ++ )
	{
		for( m = 0; m < 3; m ++ )
		{
			for( l = 0; l < 3; l ++ )
			{
				for( j = 0; j < N*N*N; j ++ )
				{
					if( (j==n) )continue;
					pass =
					    (sqrt(pow( XYZ[3*n    ] - (XYZ[ 3*j    ] + M[ k ]*L), 2)+
						    pow( XYZ[3*n + 1] - (XYZ[ 3*j + 1] + M[ m ]*L), 2)+
						    pow( XYZ[3*n + 2] - (XYZ[ 3*j + 2] + M[ l ]*L), 2)) !=0);

					if(pass == 0)break;
				}
				if(pass == 0)break;
			}
			if(pass == 0)break;
		}
		if(pass == 0)break;
	}
	return( pass );
}
/* #include"g_en_r.c"*/
