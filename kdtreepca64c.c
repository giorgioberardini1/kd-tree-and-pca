/**************************************************************************************
*
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2018/19
*
* Progetto dell'algoritmo di Product Quantization for Nearest Neighbor Search
* in linguaggio assembly x86-32 + SSE
*
* Fabrizio Angiulli, aprile 2019
*
**************************************************************************************/

/*
*
* Software necessario per l'esecuzione:
*
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
*
* entrambi sono disponibili come pacchetti software
* installabili mediante il packaging tool del sistema
* operativo; per esempio, su Ubuntu, mediante i comandi:
*
*    sudo apt-get install nasm
*    sudo apt-get install gcc
*
* potrebbe essere necessario installare le seguenti librerie:
*
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
*
* Per generare il file eseguibile:
*
* nasm -f elf32 kdtreepca32.nasm && gcc -O0 -m32 -msse kdtreepca32.o kdtreepca32c.c -o kdtreepca32c && ./kdtreepca32c
*
* oppure
*
* ./runkdtreepca32
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	MATRIX		float*

#define min(a,b) \
({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _b : _a; })

  struct nodo{
    int row;
    struct nodo* sx;
    struct nodo* dx;

  };

  #define	KDTREE	struct	nodo* // modificare con il tipo di dato utilizzato





  typedef struct {
    char* filename; //nome del file, estensione .ds per il data set, estensione .qs per l'eventuale query set
    MATRIX ds; //data set
    MATRIX qs; //query set
    int n; //numero di punti del data set
    int k; //numero di dimensioni del data/query set
    int nq; //numero di punti del query set
    int h; //numero di componenti principali da calcolare 0 se PCA non richiesta
    int kdtree_enabled; //1 per abilitare la costruzione del K-d-Tree, 0 altrimenti
    KDTREE kdtree; //riferimento al K-d-Tree, NULL se costruzione non richiesta
    float r; //raggio di query, -1 se range query non richieste
    int silent; //1 per disabilitare le stampe, 0 altrimenti
    int display; //1 per stampare i risultati, 0 altrimenti
    MATRIX U; //matrice U restituita dall'algoritmo PCA
    MATRIX V; //matrice V restituita dall'algoritmo PCA

    int k_orig; //dopo la PCA k si aggiorna ma serve l'originale per trovare il nuovo Q'
    int indiceQA; //indice del Query Answer corrente
    int capacitaQA; // capacità Query Answer

    //STRUTTURE OUTPUT MODIFICABILI
    int* QA; //risposte alle query in forma di coppie di interi (id_query, id_vicino)
    int nQA; //numero di risposte alle query

    float * medieDS; //vettore delle medie di ogni dimensione del DS
  } params;


  struct nodo* buildTree(int* D, int l, params* input, int start, int end);
  void rangeQuery(struct nodo* albero, int dir, float* q, params* input, float** h, int qi, int l, float val, float* allocato,float allocatoFiglio, int n_padding,int k_padding );
  void stampaAlbero(struct nodo *a, params* input, int* count, int show);
  int treeTest(struct nodo *a, int l, params* input);
  int contains(struct nodo *nodo, int indice);
  void espandi(params* input);

  float seleziona(float *a, int i, int n);
  void insertion_sort(float *a, int n);
  int partition(float *a, int n, float x);

  void stampa(float *a,int n);


  // PROCEDURE ASSEMBLY
  extern int prova(params* input);
  extern void euclideanDistanceAssembly(float *p, float* q, int k,float* res);
  extern void prodottoScalareAssembly(float* v, int n,float* res);
  extern void centraDatasetAssembly(float* dataset,int n,int k,float* medieDS,float n_orig);
  extern void normalizzaVAssembly(float * v, int k, float norma);
  extern void aggiornaU(float* dataset, float*u,float v, int n);

  extern void azzeraVettoreAssembly(float* v,int n);
  extern void aggiornaMatrice(float * U, float *u, int n);
  extern void prodottoDsUV2(float* dataset , float* u,int n,float* v,float UtU);
  extern void aggiornaDataset(float* dataset,float* u,float* v,int n,int k,int nXk);

  extern void queryPointASM(float * queryPoint, float *querySet, float *medieDS, int k);
  extern void queryPointASMNoPCA(float * queryPoint, float *querySet,int k);
  extern void prodottoQPV(float * queryPoint, float* V, int k,float* res);



  int padding(int dim){

    return dim%8==0 ? dim : dim+(8-dim%8);

  }


  float seleziona(float *a, int i, int n)
  {
    if(n == 1) {
      return a[0];
    }

    int n_meds = 0;
    for(int i = 0; i < n; i += 5) {
      int l = min(5, n - i);
      insertion_sort(a + i, l);

      float tmp = a[i/5];
      a[i/5] = a[i + l/2];
      a[i + l/2] = tmp;



      n_meds++;
    }

    float median_of_medians;
    if(n_meds > 1) {
      median_of_medians = seleziona(a, n_meds/2, n_meds);
    }
    else {
      median_of_medians = a[0];
    }

    int k = partition(a, n, median_of_medians);

    if(k == i) {
      return median_of_medians;
    }
    else if (i < k) {
      //scarto gli elementi da k a n
      return seleziona(a, i, k);
    }
    else {
      return seleziona(a + k, i - k, n - k);
    }
  }


  void insertion_sort(float *a, int n){
    for(int j = 1; j < n; j++) {
      float key = a[j];

      int i = j - 1;
      while ((i >= 0) && (a[i] > key)) {
        a[i + 1] = a[i];
        i--;
      }
      a[i+1] = key;
    }
  }


  int partition(float *a, int n, float x){

    for(int i = 0; i < n; i++) {
      if(a[i] == x) {
        a[i] = a[n-1];
        a[n-1] = x; //puoi fare break
      }
    }

    int i = 0;
    for(int j = 0; j < (n-1); j++) {
      if(a[j] <= x) {
        float tmp = a[j];
        a[j] = a[i];
        a[i] = tmp;
        i++;
      }
    }
    a[n-1] = a[i];
    a[i] = x;

    return i;
  }

  /*
  *
  *	Le funzioni sono state scritte assumento che le matrici siano memorizzate
  * 	mediante un array (float*), in modo da occupare un unico blocco
  * 	di memoria, ma a scelta del candidato possono essere
  * 	memorizzate mediante array di array (float**).
  *
  * 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
  * 	matrici per righe (row-major order) o per colonne (column major-order).
  *
  * 	L'assunzione corrente è che le matrici siano in row-major order.
  *
  */


  void* get_block(int size, int elements) {
    return _mm_malloc(elements*size,32);
  }


  void free_block(void* p) {
    _mm_free(p);
  }


  MATRIX alloc_matrix(int rows, int cols) {
    return (MATRIX) get_block(sizeof(float),rows*cols);
  }


  void dealloc_matrix(MATRIX mat) {
    free_block(mat);
  }


  /*
  *
  * 	load_data
  * 	=========
  *
  *	Legge da file una matrice di N righe
  * 	e M colonne e la memorizza in un array lineare in row-major order
  *
  * 	Codifica del file:
  * 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
  * 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
  * 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
  *
  *****************************************************************************
  *	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
  * 	della matrice.
  *****************************************************************************
  *
  */
 
  MATRIX load_data_orizzontale(char* filename, int *n, int *k) {
    FILE* fp;
    int rows, cols, status, i;

    fp = fopen(filename, "rb");

    if (fp == NULL){
      printf("'%s': bad data file name!\n", filename);
      exit(0);
    }

    status = fread(&cols, sizeof(int), 1, fp);
    status = fread(&rows, sizeof(int), 1, fp);

    MATRIX data = alloc_matrix(padding(rows),padding(cols));
    

    i = 0; 
    int j=0; 


    for(;i<rows;i++){
      for(j=0;j<cols;j++){
      
          status = fread(data+(j+i*(padding(cols))), sizeof(float), 1, fp);    
      }
    }



    for(i=rows;i<padding(rows);i++){
      for(j=cols;j<padding(cols);j++){     
          data[j+i*(padding(cols))]=0;
        
      }
    }

    fclose(fp);


    *n = rows;
    *k = cols;



    return data;
  }

  MATRIX load_data(char* filename, int *n, int *k) {
    FILE* fp;
    int rows, cols, status, i;

    fp = fopen(filename, "rb");

    if (fp == NULL){
      printf("'%s': bad data file name!\n", filename);
      exit(0);
    }

    status = fread(&cols, sizeof(int), 1, fp);
    status = fread(&rows, sizeof(int), 1, fp);

    MATRIX data = alloc_matrix(padding(rows), padding(cols));

    //lettura verticale
    for(int j=0; j<padding(rows); j++){
      for(int i=0; i<padding(cols); i++ ){
         if(i>=cols || j>=rows){
          data[j+i*padding(rows)]=0;
        }else{
          status = fread(data+(i*padding(rows)+j), sizeof(float), 1, fp);

        }


      }

    }

    fclose(fp);


    *n = rows;
    *k = cols;
    return data;
  }

  /*
  *
  * 	save_data
  * 	=========
  *
  *	Salva su file un array lineare in row-major order
  *	come matrice di N righe e M colonne
  *
  * 	Codifica del file:
  * 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
  * 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
  * 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
  *
  */

  void save_data(char* filename, void* X, int n, int k) {
    FILE* fp;
    int i;
    fp = fopen(filename, "wb");
    if(X != NULL){
      fwrite(&k, 4, 1, fp);
      fwrite(&n, 4, 1, fp);
      for (i = 0; i < n; i++) {
        fwrite(X, 4, k, fp);
        X += 4*k;
      }
    }
    fclose(fp);
  }

  void trasposta(params* input){
    int n= input->n;
    int k= input->k;

    MATRIX dsT= alloc_matrix(n, k);

    for(int i=0; i<k; i++){
      for(int j=0; j<n; j++){
        dsT[j*k+i] = input->ds[j+i*n];
      }
    }
    _mm_free(input->ds);
    input->ds= dsT;
  
  }

  void save_data_Verticale(char* filename, void* X, int n, int k) {

    FILE* fp;
    int i;
    int j;

    fp = fopen(filename, "wb");
    if(X != NULL){
      fwrite(&k, 4, 1, fp);
      fwrite(&n, 4, 1, fp);
      for (i = 0; i <padding(n); i++) {
        for(j =0; j<padding(k); j++){
         
          if(i>=n || j>=k){
            break;
          }

          fwrite(X+4*(j*padding(n)+i), 4, 1, fp);

        }
      }
    }
    fclose(fp);


  }




  void stampa(float *a,int n){
    for(int i=0;i<n;i++){
      printf("%f ", *(a+i));
    }
    printf("%s\n", "");
  }
  void stampaI(int *a,int n){
    for(int i=0;i<n;i++){
      printf("%d ", *(a+i));
    }
    printf("%s\n", "");
  }


  /*
  *	PCA
  * 	=====================
  */


  void pca(params* input) {


	float * ds = input->ds;

    float soglia=0.00000001;

    int k_padding = padding(input->k);
    int n_padding = padding(input->n);
    int h_padding = padding(input->h);


    input->medieDS = _mm_malloc(sizeof(float)*k_padding, 32);

    //Centro D rispetto alla media di ogni colonna:                    
    for(int i=0; i<input->k; i++){

      centraDatasetAssembly(input->ds+i*n_padding, n_padding*4, k_padding*4, input->medieDS+i, input->n);

      int start = input->n+i*(input->n+(n_padding-input->n));
  		int end = n_padding*i+n_padding;
  		for(int p=start;p<end;p++)
			  input->ds[p]=0;

    }

     float *u = _mm_malloc(n_padding*sizeof(float),32);

    queryPointASMNoPCA(u, ds, n_padding*4);

    //|v|= 1xk
    float * v = _mm_malloc(k_padding*sizeof(float), 32);

    for(int i =input->k;i<k_padding;i++){
      v[i]=0;
    }

    //alloco le matrici di output U (score) e V (load)
    MATRIX U = alloc_matrix(n_padding,h_padding);
    MATRIX V = alloc_matrix(k_padding,h_padding);

    for(int j=0; j<input->h; j++){

      float t=0;
      float tPrimo=0;
      float normV= 0;
      float prodottoUtU=0;
      float prodottoVtV=0;


      do{

        prodottoUtU = 0;

        prodottoScalareAssembly(u,n_padding*4,&prodottoUtU);

         for(int i=0; i<k_padding; i+=4){
          prodottoDsUV2(ds+n_padding*i,u,n_padding*4,v+i,prodottoUtU);
    
        }


        normV= 0;
        prodottoScalareAssembly(v, k_padding*4,&normV);

        normalizzaVAssembly(v, k_padding*4, normV);
       
        t = prodottoUtU;

        prodottoVtV= 0;
        prodottoScalareAssembly(v, k_padding*4,&prodottoVtV);


        azzeraVettoreAssembly(u,n_padding*4);

        for(int p=0; p<k_padding; p++){
            
          aggiornaU(ds+p*n_padding,u,v[p],n_padding*4);
         
        }
        normalizzaVAssembly(u, n_padding*4, prodottoVtV);

        tPrimo=0;
       
        prodottoScalareAssembly(u,n_padding*4,&tPrimo);


      }while(fabsf(tPrimo-t)>=(soglia*tPrimo));

      aggiornaMatrice(U+j*n_padding, u, n_padding*4);
     

      aggiornaMatrice(V+j*k_padding, v, k_padding*4);




      //Aggiorno il DataSet D= D-u*vT
     
      aggiornaDataset(ds,u,v,n_padding*4,k_padding*4,n_padding*k_padding*4);
     
	  }
   

    _mm_free(u);
    _mm_free(v);

    input->U = U;
    input->V = V;


    int start = input->h*k_padding;
    int end = k_padding*h_padding;

    for(int i=start;i<end;i++)
      V[i]=0;

    //Uso U come nuovo DataSet
    _mm_free(input->ds);

    input->ds = input->U;

    //successivamente mi servirà ancora il k originale
    input->k_orig = input->k;

    input->k = input->h;

  }




  void kdtree(params* input) {
    float* data = input->ds;


    struct nodo* albero = NULL;

    int * vm =_mm_malloc(input->n*sizeof(int), 32);
    for(int i=0; i<input->n; i++){
      vm[i] = i*padding(input->k);
    }

    input->kdtree = buildTree(vm, 0, input, 0, input->n);

    _mm_free(vm);

  }


  float euclideanDistance(float* p, float* q,params* input){
    // k può essere input->k o input->k_orig
    float res = 0;
    for(int i=0; i<input->k; i++){
      res+=powf(p[i]-q[i], 2);
    }
    float risC = sqrtf(res);

    return risC;

  }

  struct nodo* buildTree(int* D, int l, params* input, int start, int end){

    if(end <= start)
    return NULL;


    struct nodo* albero = NULL;
    albero = malloc(sizeof(struct nodo));

    int size = end - start;
    int c = l%input->k;
    float* data = input->ds;

    float * v= _mm_malloc(size *sizeof(float), 32);

    for (int i = start; i < end; i++) {
      int j = D[i];
      v[i-start]=data[j+c];
    }

    // calcolo mediano
     float res = seleziona(v, size/2, size);
    
    _mm_free(v);
    
    int lt = start; int gt = end-1; int indMed;

    int limite = end;
    for(int i=start; i<limite; ++i){
      int j = D[i];
      float elem = data[j+c];
      if(elem<res){
        D[lt] = j;
        lt++;
      }
      else{
        int temp = D[gt];
        D[gt] = D[i];
        D[i] = temp;
        i--;
        limite--;
        gt--;

      }
      if(res == elem){
        indMed = gt+1;
      }
    }

    gt++; //per riportare gt oltre la soglia
    if(D[gt]!=res){ //per portare il mediano alla soglia
      int temp = D[gt];
      D[gt] = D[indMed];
      D[indMed] = temp;
    }
    gt++;

    albero->row = D[lt];

    albero->sx = NULL;
    albero->sx = buildTree(D, l+1, input, start, lt);

    albero->dx = NULL;
    albero->dx = buildTree(D, l+1, input, gt, end);

  }

  void stampaAlbero(struct nodo *radice, params* input, int* numeroNodi, int show){
    //show stampa visita inOrder
    if(radice!=NULL){
      float* data = input->ds;

      stampaAlbero(radice->sx, input, numeroNodi,show);
      *numeroNodi=*numeroNodi+1;
      if(show){
        printf("[ ");
        for(int i=0; i<input->k; i++){
          printf("%f ",input->ds[i + radice->row]);
        }
        printf(" ]\n");
      }
      stampaAlbero(radice->dx, input, numeroNodi,show);
    }//if
  }//stampaAlbero

  int treeTest(struct nodo *albero, int l, params* input){
    //verifica che nodo di dimensione c abbia componente c-esimo minore/maggiore
    //del figlio sx/dx

    float* data = input->ds;
    int dim = l%input->k;
    if(albero == NULL)
    return 1;

    if((albero->dx==NULL) && (albero->sx==NULL))//se sono su una foglia
    return 1;

    if((albero->sx == NULL)){ //hai solo il figlio destro
      if(data[albero->row + dim]>data[albero->dx->row + dim])
      return 0;
    }
    else if(albero->dx==NULL){ //se ho solo figlio destro
      if(data[albero->row + dim]<=data[albero->sx->row + dim])
      return 0;
    }
    else{ //se hai entrambi i figli
      if((data[albero->row + dim]>data[albero->dx->row + dim]) || (data[albero->row + dim]<=data[albero->sx->row + dim]))
      return 0;
    }
    return treeTest(albero->sx, l+1, input)*treeTest(albero->dx, l+1, input);
  }

  int contains(struct nodo *nodo, int indice){
    if(nodo==NULL)
    return 0;
    if(nodo->row == indice)
    return 1;

    return contains(nodo->sx,indice) + contains(nodo->dx, indice);

  }

  /*
  *	Range Query Search
  * 	======================
  */
  float distance(float* q,float* puntoRegione, float* hmax, float* hmin, params* input){
    
    float * p = puntoRegione;
    
    for(int i=0; i<input->k; i++){
      if(q[i]<=hmin[i]){
        p[i] = hmin[i];
      }
      else if(q[i]>=hmax[i]){
        p[i] = hmax[i];
      }
      else p[i] = q[i];
    }

    for(int i=input->k;i<padding(input->k);i++){
          p[i]=0;
      }



    float distance=0;
    euclideanDistanceAssembly(p, q, padding(input->k)*4,&distance);
<<<<<<< HEAD
     stampa(q,20); 
=======

>>>>>>> f487ddbc67fcc0d97727044519c026dda52c2334
    return distance;
  }




  void espandi(params* input){

    
    int * array= _mm_malloc(input->capacitaQA*2*sizeof(int), 32);
    
    //copia in un array le QA
    // for(int i=0;i<input->capacitaQA; i++){
    //   array[i] = input->QA[i];
    // }


    memcpy(array,input->QA,sizeof(int)*input->capacitaQA);

    _mm_free(input->QA);
   
    input->QA=array;
    input->capacitaQA=input->capacitaQA*2;


  }



  void range_query_pca(params* input){
    float* data = input->ds; //ottengo dataset
    float* querySet = input->qs; //ottengo queryset di dimensione nQ
    struct nodo* albero = input->kdtree;

    int n_padding = padding(input->n);
    int k_padding = padding(input->k);
    int k_padding_orig = padding(input->k_orig);



    float* puntoRegione = _mm_malloc(k_padding*sizeof(float),32);

    input->capacitaQA = n_padding*2;
    input->QA = _mm_malloc(n_padding*2*sizeof(int), 32);


   
    float **h = (float**) malloc(2*sizeof(float*));  
    int max = 0; 
    int min = 1; 
    h[max] = (float*) malloc(input->k*sizeof(float*));  
    h[min] = (float*) malloc(input->k*sizeof(float*));  

    for(int j=0; j<input->k; j++){ 
        h[max][j]=data[j];
        h[min][j]=data[j];
    }

    //trova hmax e hmin
    float elem=0; 
    for(int i=0; i<input->n; i++){ 
      for(int j=0; j<input->k; j++){ 
          elem = data[i*k_padding + j];
          if(elem<h[min][j])
            h[min][j] = elem;
          else if(elem>h[max][j])
            h[max][j] = elem;
        }
    }

    float * queryPoint = _mm_malloc(sizeof(float)*k_padding_orig, 32);
    float * queryPointTmp = _mm_malloc(k_padding*sizeof(float), 32);
    float * puntQueryPoint = queryPoint; 
    for(int qi = 0; qi<input->nq; qi++){

      queryPoint=puntQueryPoint;  

<<<<<<< HEAD
     
=======
      
>>>>>>> f487ddbc67fcc0d97727044519c026dda52c2334
      queryPointASM(queryPoint, (querySet+qi*k_padding_orig), input->medieDS, k_padding_orig*4);

      //se la PCA è abilitata
      //bisogna allocare un altro queryPoint altrimenti per i prodotti riga colonna successivi al primo si perde l'effettivo queryPoint

      for(int i=0; i<k_padding; i++){

       prodottoQPV(queryPoint, (input->V+i*k_padding_orig), k_padding_orig*4,queryPointTmp+i);
<<<<<<< HEAD

=======
     
>>>>>>> f487ddbc67fcc0d97727044519c026dda52c2334
      }
      queryPoint = queryPointTmp;


<<<<<<< HEAD
      //    if(qi==30){
      //  stampa(queryPoint,20); 
      // exit(0);
      // }


     
=======
>>>>>>> f487ddbc67fcc0d97727044519c026dda52c2334


      float raggio = input->r;


      //se la regione non interseca il raggio del punto di query
      if(distance(queryPoint,puntoRegione, h[max], h[min], input) > raggio){
        return;
      }


      //se arrivo qui la regione interseca il punto di query


      float * point = data+albero->row;
<<<<<<< HEAD
      
=======
     
>>>>>>> f487ddbc67fcc0d97727044519c026dda52c2334
      //se la distanza tra il punto di query e il punto è minore del raggio
      float distance=0;
      euclideanDistanceAssembly(queryPoint, point, k_padding*4, &distance);


      if(distance<=raggio){
        // if(euclideanDistance(queryPoint, point, input)<=(raggio)){

        //se l'array ha raggiunto la capacità raddoppia

        if(input->indiceQA==input->capacitaQA){

          espandi(input);

        }

        //nota: indiceQA rappresenta la posizione già incrementata
        input->QA[input->indiceQA]=qi;
        input->indiceQA++;

        //input->QA[input->indiceQA]=albero->row;
        input->QA[input->indiceQA]=(albero->row/(k_padding));
        input->indiceQA++;
        input->nQA++;
      }


      //se vai a dx il figlio avrà hmin in posizione 0 uguale a point[0]
      //se vai a sx il figlio avrà hmax in posizione 0 uguale a point[0]
      float hFiglio = point[0];

      float hpuntoRegione = puntoRegione[0]; 


      if(albero->sx)
      //con 0 la funzione capisce che sto arrivando da sx
        rangeQuery(albero->sx, 0, queryPoint, input, h, qi, 1, hFiglio,puntoRegione,hpuntoRegione,n_padding,k_padding); //indiceQA va messa in param


      if(albero->dx)
      //con 1 la funzione capisce che sto arrivando da dx
        rangeQuery(albero->dx, 1, queryPoint, input, h, qi, 1, hFiglio,puntoRegione,hpuntoRegione,n_padding,k_padding);

    }

  }





  void range_query(params* input) {
    if(input->h>0){
      range_query_pca(input);
      return;
    }


    float* data = input->ds;     //ottengo dataset
    float* querySet = input->qs; //ottengo queryset di dimensione nQ
    struct nodo* albero = input->kdtree;

    int n_padding = padding(input->n);
    int k_padding = padding(input->k);

    float* puntoRegione = _mm_malloc(k_padding*sizeof(float),32);

    input->QA = _mm_malloc((n_padding*2)*sizeof(int), 32);
    input->capacitaQA = n_padding*2;

    float **h = (float**) malloc(2*sizeof(float*));  
    int max = 0; 
    int min = 1; 
    h[max] = (float*) malloc(input->k*sizeof(float*));  
    h[min] = (float*) malloc(input->k*sizeof(float*));  

    for(int j=0; j<input->k; j++){ 
        h[max][j]=data[j];
        h[min][j]=data[j];
    }

    //trova hmax e hmin
    float elem=0; 
    for(int i=0; i<input->n; i++){ 
      for(int j=0; j<input->k; j++){ 
          elem = data[i*k_padding + j];
          if(elem<h[min][j])
            h[min][j] = elem;
          else if(elem>h[max][j])
            h[max][j] = elem;
        }
    }

    for(int qi = 0; qi<input->nq; qi++){

		float* queryPoint = querySet+qi*k_padding;
		
      float raggio = input->r;

      //se la regione non interseca il raggio del punto di query
      if(distance(queryPoint,puntoRegione, h[max], h[min], input) > raggio)
      return;


      //se arrivo qui la regione interseca il punto di query


      float * point = data+albero->row;


      //se la distanza tra il punto di query e il punto è minore del raggio
      float distance=0;

      euclideanDistanceAssembly(queryPoint, point, k_padding*4, &distance);

      if(distance<=raggio){
        // if(euclideanDistance(queryPoint, point, input)<=(raggio)){

        //se l'array ha raggiunto la capacità raddoppia

        if(input->indiceQA==input->capacitaQA){

          espandi(input);

        }
        //nota: indiceQA rappresenta la posizione già incrementata
        input->QA[input->indiceQA]=qi;
        input->indiceQA++;


        //input->QA[input->indiceQA]=albero->row;
        input->QA[input->indiceQA]=(albero->row/(k_padding));
        input->indiceQA++;
        input->nQA++;

      }


      //se vai a dx il figlio avrà hmin in posizione 0 uguale a point[0]
      //se vai a sx il figlio avrà hmax in posizione 0 uguale a point[0]
      float hFiglio = point[0];
      float puntoRegioneFiglio = puntoRegione[0]; 

      if(albero->sx)
      //con 0 la funzione capisce che sto arrivando da sx
        rangeQuery(albero->sx, 0, queryPoint, input, h, qi, 1, hFiglio,puntoRegione,puntoRegioneFiglio,n_padding,k_padding); //indiceQA va messa in param
      
      if(albero->dx)
      //con 1 la funzione capisce che sto arrivando da dx
        rangeQuery(albero->dx, 1, queryPoint, input, h, qi, 1, hFiglio,puntoRegione,puntoRegioneFiglio, n_padding,k_padding);
     
    }
   // _mm_free(point);

  }

void rangeQuery(struct nodo* albero, int dir, float* queryPoint, params* input, float** h, int qi, int l, float valPadre, float* puntoRegione,float puntoRegionePadre, int n_padding,int k_padding){
    //qi indice Query Set
    //l è il livello corrente
   
    float prevVal;
    float prevValpuntoRegione;
    float* data = input->ds;
    //calcola dimensione
    //prendi il valore
    int c = (l-1)%input->k;

   prevValpuntoRegione=puntoRegione[c]; 
   puntoRegione[c]=puntoRegionePadre; 


    prevVal = h[dir][c]; 
    h[dir][c] = valPadre; 


    float * point = data+albero->row;
  
    if(queryPoint[c]>=h[0][c]){
      puntoRegione[c]=h[0][c]; 


    }else if(queryPoint[c]<=h[1][c]){
      puntoRegione[c]=h[1][c]; 
      
    }else{
        puntoRegione[c]=queryPoint[c]; 
    }
    float distance=0;
    euclideanDistanceAssembly(puntoRegione, queryPoint, padding(input->k)*4,&distance);

    if(distance>input->r){
   
      h[dir][c] = prevVal; 
      puntoRegione[c] = prevValpuntoRegione; 
      return; 
    }

    distance=0;
    euclideanDistanceAssembly(queryPoint, point, k_padding*4, &distance);
    if(distance<=input->r){


      if(input->indiceQA==input->capacitaQA){

        espandi(input);

      }
      input->QA[input->indiceQA]=qi;
    

      input->indiceQA++;
      
      input->QA[input->indiceQA]=(albero->row/(k_padding));
      //input->QA[input->indiceQA]=albero->row;

      input->indiceQA++;
      input->nQA++;

    }

    float hFiglio = point[l%input->k];
    float hpuntoRegione = puntoRegione[l%input->k];



    if(albero->sx){
      rangeQuery(albero->sx, 0, queryPoint, input, h, qi, l+1, hFiglio,puntoRegione,hpuntoRegione,n_padding,k_padding);
    }

    h[dir][c] = prevVal; 
    puntoRegione[c] = prevValpuntoRegione; 


    if(albero->dx){
      rangeQuery(albero->dx, 1, queryPoint, input, h, qi, l+1, hFiglio,puntoRegione,hpuntoRegione,n_padding,k_padding);

    }

  }

  int main(int argc, char** argv) {

    char fname[256];
    int i, j, k;
    clock_t t;
    float time;
    char* dsname;

    //
    // Imposta i valori di default dei parametri
    //

    params* input = malloc(sizeof(params));

    input->filename = NULL;
    input->h = 0;
    input->kdtree = NULL;
    input->r = -1;
    input->silent = 0;
    input->display = 1;
    input->QA = NULL;
    input->nQA = 0;
    input->indiceQA = 0;
    input->capacitaQA=0;


    //
    // Visualizza la sintassi del passaggio dei parametri da riga comandi
    //

    if (argc <= 1 && !input->silent) {
      printf("Usage: %s <data_name> [-pca <h>] [-kdtree [-rq <r>]]\n", argv[0]);
      printf("\nParameters:\n");
      printf("\t-d: display query results\n");
      printf("\t-s: silent\n");
      printf("\t-pca <h>: h-component PCA enabled\n");
      printf("\t-kdtree: kdtree building enabled\n");
      printf("\t-rq <r>: range query search with radius r enabled\n");
      printf("\n");
      exit(0);
    }

    //
    // Legge i valori dei parametri da riga comandi
    //


    int par = 1;
    while (par < argc) {
      if (par == 1) {
        input->filename = argv[par];
        par++;
      } else if (strcmp(argv[par],"-s") == 0) {
        input->silent = 1;
        par++;
      } else if (strcmp(argv[par],"-d") == 0) {
        input->display = 1;
        par++;
      } else if (strcmp(argv[par],"-pca") == 0) {
        par++;
        if (par >= argc) {
          printf("Missing h value!\n");
          exit(1);
        }
        input->h = atoi(argv[par]);
        par++;
      } else if (strcmp(argv[par],"-kdtree") == 0) {
        input->kdtree_enabled = 1;
        par++;
        if (par < argc && strcmp(argv[par],"-rq") == 0) {
          par++;
          if (par >= argc) {
            printf("Missing radius value!\n");
            exit(1);
          }
          input->r = atof(argv[par]);
          if(input->r < 0){
            printf("Range query radius must be non-negative!\n");
            exit(1);
          }
          par++;
        }
      } else{
        printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
        par++;
      }
    }

    //
    // Legge i dati e verifica la correttezza dei parametri
    //

    if(input->filename == NULL || strlen(input->filename) == 0){
      printf("Missing input file name!\n");
      exit(1);
    }


    sprintf(fname, "%s.ds", input->filename);
    dsname = basename(strdup(input->filename));
    char *nomeFile = dsname;
    input->ds = load_data(fname, &input->n, &input->k);
    //inizializza k_orig
    input->k_orig = input->k;

    if(input->h < 0){
      printf("Invalid value of PCA parameter h!\n");
      exit(1);
    }
    if(input->h > input->k){
      printf("Value of PCA parameter h exceeds data set dimensions!\n");
      exit(1);
    }

    if(input->r >= 0){

      sprintf(fname, "%s.qs", input->filename);
      int k;
      input->qs = load_data_orizzontale(fname, &input->nq, &k);

      if(input->k != k){
        printf("Data set dimensions and query set dimensions are not compatible!\n");
        exit(1);
      }
    }
    // Visualizza il valore dei parametri


	
    if(!input->silent){
      printf("Input file name: '%s'\n", input->filename);
      printf("Data set size [n]: %d\n", input->n);
      printf("Number of dimensions [k]: %d\n", input->k);
      if(input->h > 0){
        printf("PCA search enabled\n");
        printf("Number of principal components [h]: %i\n",input->h);
      }else{
        printf("PCA search disabled\n");
      }
      if(input->kdtree_enabled){
        printf("Kdtree building enabled\n");
        if(input->r >= 0){
          printf("Range query search enabled\n");
          printf("Range query search radius [r]: %f\n",input->r);
        }else{
          printf("Range query search disabled\n");
        }
      }else{
        printf("Kdtree building disabled\n");

      }
    }

    //
    // Calcolo PCA
    //



    if(input->h > 0){

      t = clock();
      pca(input);
      t = clock() - t;
      time = ((float)t)/CLOCKS_PER_SEC;
     
<<<<<<< HEAD
      sprintf(fname, "%s.U",dsname);
      save_data_Verticale(fname, input->U, input->n, input->h);

      sprintf(fname, "%s.V", dsname);
=======
      sprintf(fname, "%s.U", input->filename);
      save_data_Verticale(fname, input->U, input->n, input->h);

      sprintf(fname, "%s.V", input->filename);
>>>>>>> f487ddbc67fcc0d97727044519c026dda52c2334
      save_data_Verticale(fname, input->V, input->k_orig, input->h);
     

    }else{
      time = -1;
      
    }

    if (!input->silent)
    printf("\nPCA time = %.3f secs\n", time);
    else
    printf("%.3f\n", time);
    //
    // Costruzione K-d-Tree
    //


    char * nomeFile2 = (char*) malloc(20*sizeof(char));
    strcpy(nomeFile2,nomeFile);

    if(input->h>0){

      char nomeFileCompleto[] = ".U";
      strcat(nomeFile2, nomeFileCompleto);
      //trasposta(input);
      input->ds = load_data_orizzontale(nomeFile2, &input->n, &input->k);
    }
    else{

      char nomeFileCompleto[] = ".ds";
      strcat(input->filename, nomeFileCompleto);
      _mm_free(input->ds);
      
      input->ds = load_data_orizzontale(input->filename, &input->n, &input->k);
    }
<<<<<<< HEAD
=======
    input->ds = load_data_orizzontale(nomeFile2, &input->n, &input->k);
>>>>>>> f487ddbc67fcc0d97727044519c026dda52c2334


    if(input->kdtree_enabled){

      t = clock();
      kdtree(input);
      t = clock() - t;

      time = ((float)t)/CLOCKS_PER_SEC;
     
    }else{
      time = -1;
      
    }
    if (!input->silent)
    printf("\nIndexing time = %.3f secs\n", time);
    else
    printf("%.3f\n", time);

    // Range query search
    //

    if(input->r >= 0){
      t = clock();
      range_query(input);
      t = clock() - t;
      time = ((float)t)/CLOCKS_PER_SEC;
     
    }else{
      time = -1;
      
    }
    if (!input->silent)
    printf("\nQuerying time = %.3f secs\n", time);
    else
    printf("%.3f\n", time);

    dealloc_matrix(input->V);
    
    // Salva il risultato delle query
    // da modificare se si modifica il formato delle matrici di output
    //
	
    if(input->r >= 0){
      if(!input->silent && input->display) {
       
         printf("\nQuery Answer:\n");
        for(i = 0; i < input->nq; i++){
         printf("query %d: [ ", i);

         for(j = 0; j < input->nQA; j++)
         if(input->QA[j*2] == i){
         printf("%d ", input->QA[j*2+1]);

           }
         printf("]\n");
         }
         printf("\n");
      }
      sprintf(fname, "%s.qa", dsname);
      save_data(fname, input->QA, input->nQA, 2);
      printf("\nQuery Answer NUMERO:%d\n ", input->nQA);
<<<<<<< HEAD
=======

    }
>>>>>>> f487ddbc67fcc0d97727044519c026dda52c2334

    }
   
    if (!input->silent)
    printf("\nDone.\n");
<<<<<<< HEAD


    // int x; 
    // int y; 
    // MATRIX QA = alloc_matrix(1000,1435); 
    // QA=load_data("test1_15.qa",&x,&y);
    // int cont = 0; 
    // for(int i=0;i<1000*1435;i++){
    //   // printf("%f ",QA[i]); 
    //   if(QA[i]!=-1){
    //     cont++; 
    //   }
    // }

    // printf("QUERY ANS FASSETTI %i",cont);




=======
    
    stampa(input->qs, 20);
>>>>>>> f487ddbc67fcc0d97727044519c026dda52c2334

    return 0;
  }
