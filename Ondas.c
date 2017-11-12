#include <stdio.h>
#include <stdlib.h>
#include <math.h> // POR ESTA COSA TENGO QUE COMPILAR CON EL COMANDO cc Ondas.c -o Ondas.x -lm y ./Ondas.x

#define N_PUNTOS 129
#define N_PUNTOS_T 101
#define N_T 100+1 
#define c 250 /* m/s */ 
#define L 0.64 /* m */ 
#define T 1.0 /* s */ 
#define PI 3.141592

float sinf(float x);

int main ()
{
	/******** PRIMERA PARTE ********/

	/* LEE LOS DATOS */	

	FILE *datos_C; // datos es el archivo 
	datos_C = fopen("cond_ini_cuerda.dat","r");	
	
	/* GUARDA LOS DATOS EN UN ARRAY 129 x 2 */
		
	const int filas_datos_C = N_PUNTOS;
	const int cols_datos_C = 2;
	float Datos_C[filas_datos_C][2]; // Datos es la matriz de datos
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																											
	for(int j = 0; j < filas_datos_C; j++){
        	for (int i = 0; i < cols_datos_C; i++){
            		fscanf(datos_C," %f" , &Datos_C[j][i]);
			//printf("Datos_C[%d][%d]=%f", j, i, Datos_C[j][i]);
        		//printf("\n");
		}
	}

	/* ARRAYS DE CONDICIONES INICIALES */

	float x[N_PUNTOS];
	// CAMBIE DE filas_datos a N_PUNTOS-1 PARA GRAFICAR CON LA MISMA DIMENSION 
	for(int i = 0; i < N_PUNTOS; i++) {
		x[i] = Datos_C[i][0];
            	//printf("x[%d]=%f", i, x[i]);
        	//printf("\n");
    	}

	float phi_i[filas_datos_C];
	
	for(int i = 0; i < filas_datos_C; i++) {
		phi_i[i] = Datos_C[i][1];
            	//printf("phi_i[%d]=%f", i, phi_i[i]);
        	//printf("\n");
    	}

	/* VARIABLES GLOBALES */
	
	float Dx = x[1]-x[0];
	float Dt = 0.000005; 
	float alpha = c*Dt/Dx; /* condicion de estabilidad alpha < 1 */
	
	//printf("%E", alpha);

	/* DEFINO ARRAY 2D (slns de las iteraciones) */
	
	float (*phi)[N_PUNTOS] = malloc(sizeof(float[N_T][N_PUNTOS]));

	/* CONDICIONES DE FRONTERA */

	for( int j = 0; j < N_T; j++){ 
		for( int i = 0; i < N_PUNTOS; i++){
			if( i == 0  || i == N_PUNTOS-1){	
				phi[j][i] = 0.0; 
			}
		//printf("phi[%d][%d]=%f", j, i, phi[j][i]);
        	//printf("\n");
		}
	}

	/* ITERACIONES */

	for( int j = 0; j < N_T; j++){ 
		for( int i = 1; i < N_PUNTOS-1; i++){
			if( j == 0 ){	// EN t = 0 LA SLN SON LAS CONDICIONES INICIALES	
				phi[j][i] = phi_i[i]; 
			}else if( j == 1 ){  // LA PRIMERA ITERACION SON LAS CONDICIONES INICIALES
				phi[j][i] = phi[0][i] + (alpha*alpha/2.0) * (phi[0][i+1] - 2.0 * phi[0][i] + phi[0][i-1]);
			}else{
				phi[j][i] = (2.0 * (1.0-alpha*alpha)) * phi[j][i] - phi[j-1][i] + (alpha*alpha) * (phi[j][i+1] + phi[j][i-1]);
			}
	
		//printf("phi[%d][%d]=%f", j, i, phi[j][i]);
		//printf("%f", phi[1][i]-phi[0][i]);
        	//printf("\n");
		}
	}
	
	/* GENERA ARCHIVOS DE DATOS */ 

	FILE *x_x = fopen("x_x.dat", "w");
	
	for(int i = 0; i < N_PUNTOS; i++)
	{
            	fprintf(x_x, "%f ",  x[i]);
        	fprintf(x_x, "\n");
    	}
	fclose(x_x);

	FILE *phi_0 = fopen("phi_0.dat", "w");
	
	for(int i = 0; i < N_PUNTOS; i++)
	{
            	fprintf(phi_0, "%f ",  phi[0][i]);
        	fprintf(phi_0, "\n");
    	}
	fclose(phi_0);
	
	FILE *phi_12 = fopen("phi_12.dat", "w");
	
	for(int i = 0; i < N_PUNTOS; i++)
	{
            	fprintf(phi_12, "%f ",  phi[12][i]);
        	fprintf(phi_12, "\n");
    	}
	fclose(phi_12);

	FILE *phi_25 = fopen("phi_25.dat", "w");
	
	for(int i = 0; i < N_PUNTOS; i++)
	{
            	fprintf(phi_25, "%f ",  phi[25][i]);
        	fprintf(phi_25, "\n");
    	}
	fclose(phi_25);

	FILE *phi_50 = fopen("phi_50.dat", "w");
	
	for(int i = 0; i < N_PUNTOS; i++)
	{
            	fprintf(phi_50, "%f ",  phi[50][i]);
        	fprintf(phi_50, "\n");
    	}
	fclose(phi_50);

/******** ARCHIVO DE SONIDO ********/
	
	/* ARRAY DE POSICION DEL PUNTO MEDIO */

	int indice_pto_medio = 64;
	float phi_m[N_T];

	for( int j = 0; j < N_T; j++){ 
		phi_m[j] = phi[j][indice_pto_medio];
		//printf("phi_m[%d]=%f", j, phi_m[j]);
        	//printf("\n");
	}
		
	/* GENERA ARCHIVO */
	
	FILE *phi_med = fopen("phi_med.dat", "w");
	
	for(int j = 0; j < N_T; j++)
	{
            	fprintf(phi_med, "%f ",  phi_m[j]);
        	fprintf(phi_med, "\n");
    	}
	fclose(phi_med);

/******** SEGUNDA PARTE ********/

	/* DEFINO ARRAY 2D (slns de las iteraciones) */
	
	float (*psi)[N_PUNTOS] = malloc(sizeof(float[N_T][N_PUNTOS]));

	/* ARRAY DEL TIEMPO */

	float t[N_T];
	float dt = 1/N_T;
	
	for( int k = 0; k < N_T; k++){
		if(k == 0){
			t[k] = 0.0;
		}else{
			t[k] = t[k-1] + dt;
		}
	//printf("t[%d]=%f", j, t[j]);
	//printf("\n");
	}

	/* CONDICIONES DE FRONTERA */

	for( int j = 0; j < N_T; j++){ 
		for( int i = 0; i < N_PUNTOS; i++){
			if( i == 0 ){	
				phi[j][i] = 0.0; 
			}else if( i == N_PUNTOS-1 ){
				phi[j][i] = sinf((2*PI*c/L)*t[j]);
			}
		//printf("phi[%d][%d]=%f", j, i, phi[j][i]);
        	//printf("\n");
		}
	}
	
	/* ITERACIONES */

	for( int j = 0; j < N_T; j++){ 
		for( int i = 1; i < N_PUNTOS-1; i++){
			if( j == 0 ){	// EN t = 0 LA SLN SON LAS CONDICIONES INICIALES	
				psi[j][i] = phi_i[i]; 
			}else if( j == 1 ){  // LA PRIMERA ITERACION SON LAS CONDICIONES INICIALES
				psi[j][i] = psi[0][i] + (alpha*alpha/2.0) * (psi[0][i+1] - 2.0 * psi[0][i] + psi[0][i-1]);
			}else{
				psi[j][i] = (2.0 * (1.0-alpha*alpha)) * psi[j][i] - psi[j-1][i] + (alpha*alpha) * (psi[j][i+1] + psi[j][i-1]);
			}
	
		//printf("phi[%d][%d]=%f", j, i, phi[j][i]);
		//printf("%f", phi[1][i]-phi[0][i]);
        	//printf("\n");
		}
	}

	/* GENERA ARCHIVOS DE DATOS */ 

	FILE *psi_0 = fopen("psi_0.dat", "w");
	
	for(int i = 0; i < N_PUNTOS; i++)
	{
            	fprintf(psi_0, "%f ",  psi[0][i]);
        	fprintf(psi_0, "\n");
    	}
	fclose(psi_0);
	
	FILE *psi_12 = fopen("psi_12.dat", "w");
	
	for(int i = 0; i < N_PUNTOS; i++)
	{
            	fprintf(psi_12, "%f ",  psi[12][i]);
        	fprintf(psi_12, "\n");
    	}
	fclose(psi_12);

	FILE *psi_25 = fopen("psi_25.dat", "w");
	
	for(int i = 0; i < N_PUNTOS; i++)
	{
            	fprintf(psi_25, "%f ",  psi[25][i]);
        	fprintf(psi_25, "\n");
    	}
	fclose(psi_25);

	FILE *psi_50 = fopen("psi_50.dat", "w");
	
	for(int i = 0; i < N_PUNTOS; i++)
	{
            	fprintf(psi_50, "%f ",  psi[50][i]);
        	fprintf(psi_50, "\n");
    	}
	fclose(psi_50);


/******** TERCERA PARTE ********/

	/* LEE LOS DATOS */	

	FILE *datos_T; // datos es el archivo 
	datos_T = fopen("cond_ini_tambor.dat","r");
	
	/* GUARDA LOS DATOS EN UN ARRAY 29 x 29 */
		
	float Datos_T[N_PUNTOS_T][N_PUNTOS_T]; // Datos es la matriz de datos
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																											
	for(int j = 0; j < N_PUNTOS_T; j++){
        	for (int i = 0; i < N_PUNTOS_T; i++){
            		fscanf(datos_T," %f" , &Datos_T[j][i]);
			//printf("Datos_T[%d][%d]=%f", j, i, Datos_T[j][i]);
        		//printf("\n");
		}
	}
	
	/* ARRAYS DE CONDICIONES INICIALES */

	float xi_x[N_PUNTOS_T][N_PUNTOS_T];
	 
	for(int j = 0; j < N_PUNTOS_T; j++){
        	for (int i = 0; i < N_PUNTOS_T; i++){
            		xi_x[j][i] = Datos_T[j][i];
			//printf("xi_x[%d][%d]=%f", j, i, xi_x[j][i]);
        		//printf("\n");
		}
	}

	float xi_y[N_PUNTOS_T][N_PUNTOS_T];
	
	for(int i = 0; i < N_PUNTOS_T; i++){
        	for (int j = 0; j < N_PUNTOS_T; j++){
            		xi_y[j][i] = Datos_T[j][i];
			//printf("xi_y[%d][%d]=%f", j, i, xi_y[j][i]);
        		//printf("\n");
		}
	}

	/* DEFINO ARRAY 3D (slns de las iteraciones) */
	
	float *** xi;	
	xi = (float ***) malloc (sizeof (float ***)*N_T);

	for( int j = 0; j < N_PUNTOS_T; j++){ 
		xi[j] = (float **) malloc (sizeof(float ***)*N_PUNTOS_T);
		for( int i = 0; i < N_PUNTOS_T; i++){	
			xi[j][i] = (float *) malloc(sizeof(float)*N_PUNTOS_T);
		}
	}
	
	// IMPRIME EL ARRAY 3D

	/*for(int k = 0; k < N_T; k++){
		for( int j = 0; j < N_PUNTOS; j++){ 
			for( int i = 0; i < N_PUNTOS; i++){
				printf("xi[%d][%d][%d]=%f", k, j, i, xi[k][j][i]);
        			printf("\n");
			}
		}
	}*/


	/* CONDICIONES DE FRONTERA */

	for(int k = 0; k < N_T; k++){
		for( int j = 0; j < N_PUNTOS_T; j++){ 
			for( int i = 0; i < N_PUNTOS_T; i++){
				if( i == 0  || j == 0 ){	
					xi[k][j][i] = 0.0; 
				}else if(i == N_PUNTOS_T-1  || j == N_PUNTOS_T-1 ){
					xi[k][j][i] = 0.0;
				}
			}
		}
	}
		
	float betha = alpha; 
	//printf("%E", alpha);

	/* CONDICIONES INICIALES */

	for(int k = 0; k < N_T; k++){
		for(int j = 0; j < N_PUNTOS_T; j++){
        		for (int i = 0; i < N_PUNTOS_T; i++){
				if( k == 0 ){
					xi[k][j][i] = xi_x[j][i];
				}
            		 
			//printf("xi[%d][%d][%d]=%f", k, j, i, xi[k][j][i]);
        		//printf("\n");
			}
		}
	}

	for(int k = 0; k < N_T; k++){
		for(int i = 0; i < N_PUNTOS_T; i++){
        		for (int j = 0; j < N_PUNTOS_T; j++){
            			if( k == 0 ){
					xi[k][j][i] = xi_y[j][i];
				}
			//printf("xi[%d][%d][%d]=%f", k, j, i, xi[0][j][i]);
			//printf("%f", xi[k][j][i]-xi_y[j][i]);
        		//printf("\n");
			}
		}
	}

	/* ITERACIONES */ 

	for(int k = 1; k < N_T; k++){
		for( int j = 1; j < N_PUNTOS_T-1; j++){ 
			for( int i = 1; i < N_PUNTOS_T-1; i++){
				if( k == 1 ){ // PRIMER PASO 
					xi[k][j][i] = 2*xi[0][j][i] + (betha*betha)*(xi[0][j+1][i] + xi[0][j][i+1] - 4*xi[0][j][i] + xi[0][j-1][i] + xi[0][j][i-1]);
				}else{ // PASO PARA EL TIEMPO ARBITRARIO
					xi[k][j][i] = 2*(1-2*(betha*betha))*xi[k][j][i] - xi[k-1][j][i] + (betha*betha)*(xi[k][j+1][i] + xi[k][j][i+1] + xi[k][j-1][i] + xi[k][j][i-1]);
				}
			//printf("xi[%d][%d][%d]=%f", k, j, i, xi[k][j][i]);
        		//printf("\n");
			}
		}
	}

	/* GENERA ARCHIVOS DE DATOS */ 

	// TIEMPO 0

	FILE *xi_x0 = fopen("xi_x0.dat", "w"); 
	
	for(int j = 0; j < N_PUNTOS_T; j++){
        	for (int i = 0; i < N_PUNTOS_T; i++){
			fprintf(xi_x0, "%f ",  xi[0][j][i]);
        		
		}
	fprintf(xi_x0, "\n");
	}
	fclose(xi_x0);

	FILE *xi_y0 = fopen("xi_y0.dat", "w"); 
	
	for(int j = 0; j < N_PUNTOS_T; j++){
        	for (int i = 0; i < N_PUNTOS_T; i++){
			fprintf(xi_y0, "%f ",  xi[0][j][i]);
		}
	fprintf(xi_y0, "\n");
	}
	fclose(xi_y0);

	// TIEMPO 1/8T 

	FILE *xi_x12 = fopen("xi_x12.dat", "w"); 
	
	for(int j = 0; j < N_PUNTOS_T; j++){
        	for (int i = 0; i < N_PUNTOS_T; i++){
			fprintf(xi_x12, "%f ",  xi[12][j][i]);
		}
	fprintf(xi_x12, "\n");
	}
	fclose(xi_x12);

	FILE *xi_y12 = fopen("xi_y12.dat", "w"); 
	
	for(int j = 0; j < N_PUNTOS_T; j++){
        	for (int i = 0; i < N_PUNTOS_T; i++){
			fprintf(xi_y12, "%f ",  xi[12][j][i]);
		}
	fprintf(xi_y12, "\n");
	}
	fclose(xi_y12);

	// TIEMPO 1/4T 

	FILE *xi_x25 = fopen("xi_x25.dat", "w"); 
	
	for(int j = 0; j < N_PUNTOS_T; j++){
        	for (int i = 0; i < N_PUNTOS_T; i++){
			fprintf(xi_x25, "%f ",  xi[25][j][i]);
		}
	fprintf(xi_x25, "\n");
	}
	fclose(xi_x25);

	FILE *xi_y25 = fopen("xi_y25.dat", "w"); 
	
	for(int j = 0; j < N_PUNTOS_T; j++){
        	for (int i = 0; i < N_PUNTOS_T; i++){
			fprintf(xi_y25, "%f ",  xi[25][j][i]);
		}
	fprintf(xi_y25, "\n");
	}
	fclose(xi_y25);

	// TIEMPO 1/2T 

	FILE *xi_x50 = fopen("xi_x50.dat", "w"); 
	
	for(int j = 0; j < N_PUNTOS_T; j++){
        	for (int i = 0; i < N_PUNTOS_T; i++){
			fprintf(xi_x50, "%f ",  xi[50][j][i]);
		}
	fprintf(xi_x50, "\n");
	}
	fclose(xi_x50);

	FILE *xi_y50 = fopen("xi_y50.dat", "w"); 
	
	for(int j = 0; j < N_PUNTOS_T; j++){
        	for (int i = 0; i < N_PUNTOS_T; i++){
			fprintf(xi_y50, "%f ",  xi[50][j][i]);
		}
	fprintf(xi_y50, "\n");
	}
	fclose(xi_y50);
	
	return 0;

}


