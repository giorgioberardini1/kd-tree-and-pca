#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define MATRIX float *

#define min(a, b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _b : _a; })

struct node
{
  int rowIndex; //row index of the dataset
  struct node *sx;
  struct node *dx;
};

#define KDTREE struct node * 

typedef struct
{
  char *filename;     //.ds extension for dataset, .qs estension for the a possible query set
  MATRIX ds;          //data set
  MATRIX qs;          //query set
  int n;              //dataset's number of points
  int k;              //dataset's number of dimension
  int nq;             //query set's number of points
  int h;              //number of principal component (PCA option required)
  int kdtree_enabled; 
  KDTREE kdtree;      //reference to K-d-Tree, NULL if the kd-tree bulding isn't required
  float r;            //query radius, -1  if range query isn't required 
  int silent;        
  int display;        
  MATRIX U;           //matrix U computed by PCA algorithm
  MATRIX V;           //matrix V computed by PCA algorithm
  int k_orig;         //initial value of k (after PCA changes)
  int indexQA;        //actual index of the current Query Answer
  int capacityQA;     //actual capacity of Query Answer's data structure
  int *QA;            //reference to array of query answer in the following format (id_query, id_neighbour)
  int nQA;            //number of query answer 
  float *meanDS;      //array of float, each element is the mean of each dataset's column
} params;

struct node *buildTree(int *D, int l, params *input, int start, int end);
void rangeQuery(struct node *albero, int dir, float *q, params *input, float **h, int qi, int l, float val, float *regionPoint, float regionPointChild, int n_padding, int k_padding);
void printTree(struct node *root, params *input, int *count, int show);
int treeTest(struct node *a, int l, params *input);
int contains(struct node *node, int index);
void expand(params *input);
float mom(float *a, int i, int n); //median's of median algotithm
void insertion_sort(float *a, int n);
int partition(float *a, int n, float x);





//ASSEMBLY FUNCTIONS
extern void euclideanDistanceAssembly(float *p, float *q, int k, float *res);
extern void scalarProductAssembly(float *v, int n, float *res);
extern void centerDatasetAssembly(float *dataset, int n, int k, float *meanDS, float n_orig);
extern void normalizeVAssembly(float *v, int k, float norma);
extern void updateU(float *dataset, float *u, float v, int n);
extern void resetVectorAssembly(float *v, int n);
extern void updateMatrix(float *U, float *u, int n);
extern void prodottoDsUV2(float *dataset, float *u, int n, float *v, float UtU);
extern void updateDataset(float *dataset, float *u, float *v, int n, int k, int nXk);
extern void queryPointASM(float *queryPoint, float *querySet, float *meanDS, int k);
extern void queryPointASMNoPCA(float *queryPoint, float *querySet, int k);
extern void productQPV(float *queryPoint, float *V, int k, float *res);

int padding(int dim)
{
  return dim % 4 == 0 ? dim : dim + (4 - dim % 4);
}

float mom(float *a, int i, int n)
{
  if (n == 1)
  {
    return a[0];
  }

  int n_meds = 0;
  for (int i = 0; i < n; i += 5)
  {
    int l = min(5, n - i);
    insertion_sort(a + i, l);

    float tmp = a[i / 5];
    a[i / 5] = a[i + l / 2];
    a[i + l / 2] = tmp;

    n_meds++;
  }

  float median_of_medians;
  if (n_meds > 1)
  {
    median_of_medians = mom(a, n_meds / 2, n_meds);
  }
  else
  {
    median_of_medians = a[0];
  }

  int k = partition(a, n, median_of_medians);

  if (k == i)
  {
    return median_of_medians;
  }
  else if (i < k)
  {
    
    return mom(a, i, k);
  }
  else
  {
    return mom(a + k, i - k, n - k);
  }
}

void insertion_sort(float *a, int n)
{
  for (int j = 1; j < n; j++)
  {
    float key = a[j];

    int i = j - 1;
    while ((i >= 0) && (a[i] > key))
    {
      a[i + 1] = a[i];
      i--;
    }
    a[i + 1] = key;
  }
}

int partition(float *a, int n, float x)
{

  for (int i = 0; i < n; i++)
  {
    if (a[i] == x)
    {
      a[i] = a[n - 1];
      a[n - 1] = x; 
    }
  }

  int i = 0;
  for (int j = 0; j < (n - 1); j++)
  {
    if (a[j] <= x)
    {
      float tmp = a[j];
      a[j] = a[i];
      a[i] = tmp;
      i++;
    }
  }
  a[n - 1] = a[i];
  a[i] = x;

  return i;
}

void *get_block(int size, int elements)
{
  return _mm_malloc(elements * size, 16);
}

void free_block(void *p)
{
  _mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols)
{
  return (MATRIX)get_block(sizeof(float), rows * cols);
}

void dealloc_matrix(MATRIX mat)
{
  free_block(mat);
}

MATRIX load_data_horizontal(char *filename, int *n, int *k)
{

  //load data from file considering and stores it in row-major format 

  FILE *fp;
  int rows, cols, status, i;

  fp = fopen(filename, "rb");

  if (fp == NULL)
  {
    printf("'%s': bad data file name!\n", filename);
    exit(0);
  }

  status = fread(&cols, sizeof(int), 1, fp);
  status = fread(&rows, sizeof(int), 1, fp);

  MATRIX data = alloc_matrix(padding(rows), padding(cols));

  i = 0;
  int j = 0;

  for (; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
    {

      status = fread(data + (j + i * (padding(cols))), sizeof(float), 1, fp);
    }
  }

  for (i = rows; i < padding(rows); i++)
  {
    for (j = cols; j < padding(cols); j++)
    {
      data[j + i * (padding(cols))] = 0;
    }
  }

  fclose(fp);

  *n = rows;
  *k = cols;

  return data;
}

MATRIX load_data(char *filename, int *n, int *k)
{
  FILE *fp;
  int rows, cols, status, i;

  fp = fopen(filename, "rb");

  if (fp == NULL)
  {
    printf("'%s': bad data file name!\n", filename);
    exit(0);
  }

  status = fread(&cols, sizeof(int), 1, fp);
  status = fread(&rows, sizeof(int), 1, fp);

  MATRIX data = alloc_matrix(padding(rows), padding(cols));

  //vertical reading
  for (int j = 0; j < padding(rows); j++)
  {
    for (int i = 0; i < padding(cols); i++)
    {
      if (i >= cols || j >= rows)
      {
        data[j + i * padding(rows)] = 0;
      }
      else
      {
        status = fread(data + (i * padding(rows) + j), sizeof(float), 1, fp);
      }
    }
  }

  fclose(fp);

  *n = rows;
  *k = cols;

  return data;
}



void save_data(char *filename, void *X, int n, int k)
{

  //save data in row-major format
  FILE *fp;
  int i;
  fp = fopen(filename, "wb");
  if (X != NULL)
  {
    fwrite(&k, 4, 1, fp);
    fwrite(&n, 4, 1, fp);
    for (i = 0; i < n; i++)
    {
      fwrite(X, 4, k, fp);
      X += 4 * k;
    }
  }
  fclose(fp);
}



void save_data_vertical(char *filename, void *X, int n, int k)
{

  FILE *fp;
  int i;
  int j;

  fp = fopen(filename, "wb");
  if (X != NULL)
  {
    fwrite(&k, 4, 1, fp);
    fwrite(&n, 4, 1, fp);
    for (i = 0; i < padding(n); i++)
    {
      for (j = 0; j < padding(k); j++)
      {


        if (i >= n || j >= k)
        {
          break;
        }

        fwrite(X + 4 * (j * padding(n) + i), 4, 1, fp);
      }
    }
  }
  fclose(fp);
}


void pca(params *input)
{

  float *ds = input->ds;

  float treshold = 0.00000001;

  int k_padding = padding(input->k); //k after padding
  int n_padding = padding(input->n);
  int h_padding = padding(input->h);

  input->meanDS = _mm_malloc(sizeof(float) * k_padding, 16);

  //Center Dataset with the mean of each column
  for (int i = 0; i < k_padding; i++)
  {

    centerDatasetAssembly(input->ds + i * n_padding, n_padding * 4, k_padding * 4, input->meanDS + i, input->n);
    int start = input->n + i * (input->n + (n_padding - input->n));
    int end = n_padding * i + n_padding;
    for (int p = start; p < end; p++)
      input->ds[p] = 0;
  }

  //u is the first column
  float *u = _mm_malloc(n_padding * sizeof(float), 16);

  queryPointASMNoPCA(u, ds, n_padding * 4);

  float *v = _mm_malloc(k_padding * sizeof(float), 16);

  for (int i = input->k; i < k_padding; i++)
  {
    v[i] = 0;
  }

  //allocate output matrices U (score) and V (load)
  MATRIX U = alloc_matrix(n_padding, h_padding);
  MATRIX V = alloc_matrix(k_padding, h_padding);

  for (int j = 0; j < input->h; j++)
  {

    float t = 0;
    float t2 = 0;
    float normV = 0;
    float uTuProduct = 0;
    float vTvProduct = 0;

    do
    {

      uTuProduct = 0;

      scalarProductAssembly(u, n_padding * 4, &uTuProduct);

      //Compute load vector
      for (int i = 0; i < k_padding; i += 4)
      {
        prodottoDsUV2(ds + n_padding * i, u, n_padding * 4, v + i, uTuProduct);
      }

      normV = 0;
      scalarProductAssembly(v, k_padding * 4, &normV);

      //normalize load vector
      normalizeVAssembly(v, k_padding * 4, normV);

      t = uTuProduct;

      //calculate vT * v
      vTvProduct = 0;
      scalarProductAssembly(v, k_padding * 4, &vTvProduct);

      resetVectorAssembly(u, n_padding * 4);

      for (int p = 0; p < input->k; p++)
      {
        updateU(ds + p * n_padding, u, v[p], n_padding * 4);
      }

      normalizeVAssembly(u, n_padding * 4, vTvProduct);

      //calculate t2=uT*u 
      t2 = 0;
      scalarProductAssembly(u, n_padding * 4, &t2);


    } while (fabsf(t2 - t) >= (treshold * t2));

    //insert u as j-th column of U
    updateMatrix(U + j * n_padding, u, n_padding * 4);

    //insert v as j-th column of V
    updateMatrix(V + j * k_padding, v, k_padding * 4);

    //Update Dataset D = D - u*vT
    updateDataset(ds, u, v, n_padding * 4, k_padding * 4, n_padding * k_padding * 4);

  } //for h

  _mm_free(u);
  _mm_free(v);

  input->U = U;
  input->V = V;

  int start = input->h * input->k;
  int end = k_padding * h_padding;

  for (int i = start; i < end; i++)
    V[i] = 0;

  
  _mm_free(input->ds);

  //U is the new dataset
  input->ds = input->U;

  //store the original value of k 
  input->k_orig = input->k;

  input->k = input->h;
}

//build the tree
void kdtree(params *input)
{
  float *data = input->ds;

  // struct node *root = NULL;

  int *vm = _mm_malloc(input->n * sizeof(int), 16);


  for (int i = 0; i < input->n; i++)
  {
    vm[i] = i * padding(input->k);
  }

  
  input->kdtree = buildTree(vm, 0, input, 0, input->n);

  _mm_free(vm);
}

//not used
float euclideanDistance(float *p, float *q, params *input)
{
  
  float res = 0;
  for (int i = 0; i < input->k; i++)
  {
    res += powf(p[i] - q[i], 2);
  }
  float risC = sqrtf(res);

  return risC;
}



//buld the tree recursively
struct node *buildTree(int *D, int l, params *input, int start, int end)
{


  //D -> represent the index of the column's elements
  //l -> is the current level
  //start -> is the start index of D 
  //end -> is the end index of D 

  //for memory managment purpose we prefferred to allocate D once (outside buildTree function) 
  //and reuses it

  if (end <= start)
    return NULL;

  struct node *root = malloc(sizeof(struct node));

  int size = end - start;
  int c = l % input->k;
  float *dataset = input->ds;

  float *v = _mm_malloc(size * sizeof(float), 16);

  for (int i = start; i < end; i++)
  {
    int j = D[i];
    v[i - start] = dataset[j + c];
  }

  //find median using medians of median algorithm 
  float res = mom(v, size / 2, size);

  _mm_free(v);

  int lt = start; 
  int gt = end - 1;
  int indMed;

  int limit = end;
  for (int i = start; i < limit; ++i)
  {
    int j = D[i];
    float elem = dataset[j + c];
    if (elem < res)
    {
      D[lt] = j;
      lt++;
    }
    else
    {
      int temp = D[gt];
      D[gt] = D[i];
      D[i] = temp;
      i--;
      limit--;
      gt--;
    }
    if (res == elem)
    {
      indMed = gt + 1;
    }
  }

  gt++; 
  if (D[gt] != res)
  { 
    int temp = D[gt];
    D[gt] = D[indMed];
    D[indMed] = temp;
  }
  gt++;

  root->rowIndex = D[lt];

  root->sx = buildTree(D, l + 1, input, start, lt);

  root->dx = buildTree(D, l + 1, input, gt, end);
}

void printTree(struct node *root, params *input, int *numberOfNodes, int show)
{
  //show stampa visita inOrder
  if (root != NULL)
  {
    float *data = input->ds;

    printTree(root->sx, input, numberOfNodes, show);
    *numberOfNodes = *numberOfNodes + 1;
    if (show)
    {
      printf("[ ");
      for (int i = 0; i < input->k; i++)
      {
        printf("%f ", input->ds[i + root->rowIndex]);
      }
      printf(" ]\n");
    }
    printTree(root->dx, input, numberOfNodes, show);
  } //if
} //printTree



int treeTest(struct node *root, int l, params *input)
{
  
  //check if the node of dimension c has c-th component less/greater of the left/right child
  
  float *data = input->ds;
  int dim = l % input->k;
  if (root == NULL)
    return 1;

  if ((root->dx == NULL) && (root->sx == NULL)) //leaf node
    return 1;

  if ((root->sx == NULL))//only rigth child
  { 
    if (data[root->rowIndex + dim] > data[root->dx->rowIndex + dim])
      return 0;
  }
  else if (root->dx == NULL) //only left child
  { 
    if (data[root->rowIndex + dim] <= data[root->sx->rowIndex + dim])
      return 0;
  }
  else
  { //se hai entrambi i figli
    if ((data[root->rowIndex + dim] > data[root->dx->rowIndex + dim]) || (data[root->rowIndex + dim] <= data[root->sx->rowIndex + dim]))
      return 0;
  }
  return treeTest(root->sx, l + 1, input) * treeTest(root->dx, l + 1, input);
}

int contains(struct node *node, int index)
{
  if (node == NULL)
    return 0;
  if (node->rowIndex == index)
    return 1;

  return contains(node->sx, index) + contains(node->dx, index);
}




float distance(float *q, float *pointRegion, float *hmax, float *hmin, params *input)
{

  float *p = pointRegion;

  for (int i = 0; i < input->k; i++)
  {
    if (q[i] <= hmin[i])
    {
      p[i] = hmin[i];
    }
    else if (q[i] >= hmax[i])
    {
      p[i] = hmax[i];
    }
    else
      p[i] = q[i];
  }

  for (int i = input->k; i < padding(input->k); i++)
  {
    p[i] = 0;
  }

  float distance = 0;
  euclideanDistanceAssembly(p, q, padding(input->k) * 4, &distance);

  return distance;
}




void expand(params *input)
{

  int *array = _mm_malloc(input->capacityQA * 2 * sizeof(int), 16);

  memcpy(array, input->QA, sizeof(int) * input->capacityQA);

  _mm_free(input->QA);

  input->QA = array;
  input->capacityQA = input->capacityQA * 2;
}




void range_query_pca(params *input)
{

  float *data = input->ds;    
  float *querySet = input->qs; 
  struct node *root = input->kdtree;

  int n_padding = padding(input->n);
  int k_padding = padding(input->k);
  int k_padding_orig = padding(input->k_orig);

  float *regionPoint = _mm_malloc(k_padding * sizeof(float), 16);

  input->QA = _mm_malloc(n_padding * 2 * sizeof(int), 16);
  input->capacityQA = n_padding * 2;

  float **h = (float **)malloc(2 * sizeof(float *));
  int max = 0;
  int min = 1;
  h[max] = (float *)malloc(input->k * sizeof(float *));
  h[min] = (float *)malloc(input->k * sizeof(float *));

  for (int j = 0; j < input->k; j++)
  {
    h[max][j] = data[j];
    h[min][j] = data[j];
  }

  //find hmax e hmin
  float elem = 0;
  for (int i = 0; i < input->n; i++)
  { 
    for (int j = 0; j < input->k; j++)
    { 
      elem = data[i * k_padding + j];
      if (elem < h[min][j])
        h[min][j] = elem;
      else if (elem > h[max][j])
        h[max][j] = elem;
    }
  }

  

  float *queryPoint = _mm_malloc(sizeof(float) * k_padding_orig, 16);
  float *queryPointTmp = _mm_malloc(k_padding * sizeof(float), 16);
  float *puntQueryPoint = queryPoint;
  for (int qi = 0; qi < input->nq; qi++)
  {
    queryPoint = puntQueryPoint;

    queryPointASM(queryPoint, (querySet + qi * k_padding_orig), input->meanDS, k_padding_orig * 4);

    for (int i = 0; i < k_padding; i++)
    {

      productQPV(queryPoint, (input->V + i * k_padding_orig), k_padding_orig * 4, queryPointTmp + i);
    }

    queryPoint = queryPointTmp;

    float raggio = input->r;

   
    if (distance(queryPoint, regionPoint, h[max], h[min], input) > raggio)
    {
      return;
    }

    

    float *point = data + root->rowIndex;

    float distance = 0;

    euclideanDistanceAssembly(queryPoint, point, k_padding * 4, &distance);

    if (distance <= raggio)
    {
      // if(euclideanDistance(queryPoint, point, input)<=(raggio)){

      

      if (input->indexQA == input->capacityQA)
      {

        expand(input);
      }

      //add the point to the query answer's array
      input->QA[input->indexQA] = qi;
      input->indexQA++;

      
      input->QA[input->indexQA] = (root->rowIndex / (k_padding));
      input->indexQA++;
      input->nQA++;
    }


    float hChild = point[0];
    float hAllocate = regionPoint[0];

    if (root->sx)
      
      rangeQuery(root->sx, 0, queryPoint, input, h, qi, 1, hChild, regionPoint, hAllocate, n_padding, k_padding); 

    if (root->dx)
     
      rangeQuery(root->dx, 1, queryPoint, input, h, qi, 1, hChild, regionPoint, hAllocate, n_padding, k_padding);
  }
}

void range_query(params *input)
{

  if (input->h > 0)
  {
    range_query_pca(input);
    return;
  }

  float *dataset = input->ds;     
  float *querySet = input->qs;
  struct node *root = input->kdtree;

  int n_padding = padding(input->n); 
  int k_padding = padding(input->k);

  float radius = input->r;
  float *regionPoint = _mm_malloc(k_padding * sizeof(float), 16);

  input->capacityQA = n_padding * 2;
  input->QA = _mm_malloc(n_padding * 2 * sizeof(int), 16);

  float **h = (float **)malloc(2 * sizeof(float *));
  int max = 0;
  int min = 1;
  h[max] = (float *)malloc(input->k * sizeof(float *));
  h[min] = (float *)malloc(input->k * sizeof(float *));

  for (int j = 0; j < input->k; j++)
  {
    h[max][j] = dataset[j];
    h[min][j] = dataset[j];
  }

  //find hmax and hmin
  float elem = 0;
  for (int i = 0; i < input->n; i++)
  { 
    for (int j = 0; j < input->k; j++)
    { 
      elem = dataset[i * k_padding + j];
      if (elem < h[min][j])
        h[min][j] = elem;
      else if (elem > h[max][j])
        h[max][j] = elem;
    }
  }

  
  for (int qi = 0; qi < input->nq; qi++)
  {

    float *queryPoint = querySet + qi * k_padding;

 
    //if the region doesn't intersect the radius around the query point
    if (distance(queryPoint, regionPoint, h[max], h[min], input) > radius)
      return;

    float *point = dataset + root->rowIndex;

    float distance = 0;

    euclideanDistanceAssembly(queryPoint, point, k_padding * 4, &distance);
  
    //if the distance beetween query point and point is less than radius
    if (distance <= radius)
    {
      // if(euclideanDistance(queryPoint, point, input)<=(raggio)){

      //if the array has reached the capacity -> double it 

      if (input->indexQA == input->capacityQA)
      {

        expand(input);
      }

      
      input->QA[input->indexQA] = qi;
      input->indexQA++;

      //input->QA[input->indexQA]=albero->rowIndex;
      input->QA[input->indexQA] = (root->rowIndex / (k_padding));
      input->indexQA++;
      input->nQA++;
    }

    
    float hChild = point[0];
    float regionPointChild = regionPoint[0];

    if (root->sx)
      
      rangeQuery(root->sx, 0, queryPoint, input, h, qi, 1, hChild, regionPoint, regionPointChild, n_padding, k_padding); 

    if (root->dx)
      rangeQuery(root->dx, 1, queryPoint, input, h, qi, 1, hChild, regionPoint, regionPointChild, n_padding, k_padding);
  }
  // _mm_free(point);
}





//dir = 0 if left child 
//dir = 1 if right child 
void rangeQuery(struct node *currentNode, int dir, float *queryPoint, params *input, float **h, int qi, int l, float parentValue, float *regionPoint, float regionPointFather, int n_padding, int k_padding)
{

  //qi Query Set index
  //l is the current level
  //c is the 

  float prevVal;
  float prevValAllocate;
  float *data = input->ds;

  //calculate cutting size
  int c = (l - 1) % input->k;

  prevValAllocate = regionPoint[c];
  regionPoint[c] = regionPointFather;

  prevVal = h[dir][c];
  h[dir][c] = parentValue;

  float *point = data + currentNode->rowIndex;

  if (queryPoint[c] >= h[0][c])
  {
    regionPoint[c] = h[0][c];
  }
  else if (queryPoint[c] <= h[1][c])
  {
    regionPoint[c] = h[1][c];
  }
  else
  {
    regionPoint[c] = queryPoint[c];
  }

  float distance = 0;
  euclideanDistanceAssembly(regionPoint, queryPoint, padding(input->k) * 4, &distance);

  if (distance > input->r)
  {

    h[dir][c] = prevVal;
    regionPoint[c] = prevValAllocate;
    return;
  }

  distance = 0;
  euclideanDistanceAssembly(queryPoint, point, k_padding * 4, &distance);
  if (distance <= input->r)
  {

    if (input->indexQA == input->capacityQA)
    {
      expand(input);
    }

    input->QA[input->indexQA] = qi;

    input->indexQA++;

    
    input->QA[input->indexQA] = (currentNode->rowIndex / (k_padding));

    input->indexQA++;
    input->nQA++;
  }

  float hChild= point[l % input->k];
  float hAllocate = regionPoint[l % input->k];

  if (currentNode->sx)
  {
    rangeQuery(currentNode->sx, 0, queryPoint, input, h, qi, l + 1, hChild, regionPoint, hAllocate, n_padding, k_padding);
  }

  h[dir][c] = prevVal;
  regionPoint[c] = prevValAllocate;

  if (currentNode->dx)
  {
    rangeQuery(currentNode->dx, 1, queryPoint, input, h, qi, l + 1, hChild, regionPoint, hAllocate, n_padding, k_padding);
  }
}

int main(int argc, char **argv)
{

  char fname[256];
  int i, j, k;
  clock_t t;
  float time;
  char *dsname;


  params *input = malloc(sizeof(params));

  input->filename = NULL;
  input->h = 0;
  input->kdtree = NULL;
  input->r = -1;
  input->silent = 0;
  input->display = 1;
  input->QA = NULL;
  input->nQA = 0;
  input->indexQA = 0;
  input->capacityQA = 0;


  if (argc <= 1 && !input->silent)
  {
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


  int par = 1;
  while (par < argc)
  {
    if (par == 1)
    {
      input->filename = argv[par];
      par++;
    }
    else if (strcmp(argv[par], "-s") == 0)
    {
      input->silent = 1;
      par++;
    }
    else if (strcmp(argv[par], "-d") == 0)
    {
      input->display = 1;
      par++;
    }
    else if (strcmp(argv[par], "-pca") == 0)
    {
      par++;
      if (par >= argc)
      {
        printf("Missing h value!\n");
        exit(1);
      }
      input->h = atoi(argv[par]);
      par++;
    }
    else if (strcmp(argv[par], "-kdtree") == 0)
    {
      input->kdtree_enabled = 1;
      par++;
      if (par < argc && strcmp(argv[par], "-rq") == 0)
      {
        par++;
        if (par >= argc)
        {
          printf("Missing radius value!\n");
          exit(1);
        }
        input->r = atof(argv[par]);
        if (input->r < 0)
        {
          printf("Range query radius must be non-negative!\n");
          exit(1);
        }
        par++;
      }
    }
    else
    {
      printf("WARNING: unrecognized parameter '%s'!\n", argv[par]);
      par++;
    }
  }

  if (input->filename == NULL || strlen(input->filename) == 0)
  {
    printf("Missing input file name!\n");
    exit(1);
  }

  sprintf(fname, "%s.ds", input->filename);
  dsname = basename(strdup(input->filename));
  char *nomeFile = dsname;
  input->ds = load_data(fname, &input->n, &input->k);

  //initialize k_orig 
  input->k_orig = input->k;

  if (input->h < 0)
  {
    printf("Invalid value of PCA parameter h!\n");
    exit(1);
  }
  if (input->h > input->k)
  {
    printf("Value of PCA parameter h exceeds data set dimensions!\n");
    exit(1);
  }

  if (input->r >= 0)
  {

    sprintf(fname, "%s.qs", input->filename);
    int k;
    input->qs = load_data_horizontal(fname, &input->nq, &k);

    if (input->k != k)
    {
      printf("Data set dimensions and query set dimensions are not compatible!\n");
      exit(1);
    }
  }

  

  if (!input->silent)
  {
    printf("Input file name: '%s'\n", input->filename);
    printf("Data set size [n]: %d\n", input->n);
    printf("Number of dimensions [k]: %d\n", input->k);
    if (input->h > 0)
    {
      printf("PCA search enabled\n");
      printf("Number of principal components [h]: %i\n", input->h);
    }
    else
    {
      printf("PCA search disabled\n");
    }
    if (input->kdtree_enabled)
    {
      printf("Kdtree building enabled\n");
      if (input->r >= 0)
      {
        printf("Range query search enabled\n");
        printf("Range query search radius [r]: %f\n", input->r);
      }
      else
      {
        printf("Range query search disabled\n");
      }
    }
    else
    {
      printf("Kdtree building disabled\n");
    }
  }


  if (input->h > 0)
  {

    t = clock();
    pca(input);
    t = clock() - t;
    time = ((float)t) / CLOCKS_PER_SEC;

   
    sprintf(fname, "%s.U", dsname);
    save_data_vertical(fname, input->U, input->n, input->h);
    sprintf(fname, "%s.V", dsname);
    save_data_vertical(fname, input->V, input->k_orig, input->h);
  }
  else
  {
    time = -1;
  }

  if (!input->silent)
    printf("\nPCA time = %.3f secs\n", time);
  else
    printf("%.3f\n", time);


  char *nomeFile2 = (char *)malloc(20 * sizeof(char));
  strcpy(nomeFile2, nomeFile);

  if (input->h > 0)
  {

    char completeFileName[] = ".U";
    strcat(nomeFile2, completeFileName);
    //trasposta(input);
    input->ds = load_data_horizontal(nomeFile2, &input->n, &input->k);
  }
  else
  {

    char completeFileName[] = ".ds";
    strcat(input->filename, completeFileName);
    _mm_free(input->ds);

    input->ds = load_data_horizontal(input->filename, &input->n, &input->k);
  }

  if (input->kdtree_enabled)
  {

    t = clock();
    kdtree(input);
    t = clock() - t;

    time = ((float)t) / CLOCKS_PER_SEC;
  }
  else
  {
    time = -1;
  }
  if (!input->silent)
    printf("\nIndexing time = %.3f secs\n", time);
  else
    printf("%.3f\n", time);

  // Range query search

  if (input->r >= 0)
  {
    t = clock();
    range_query(input);
    t = clock() - t;
    time = ((float)t) / CLOCKS_PER_SEC;
  }
  else
  {
    time = -1;
  }
  if (!input->silent)
    printf("\nQuerying time = %.3f secs\n", time);
  else
    printf("%.3f\n", time);

  dealloc_matrix(input->V);

  if (input->r >= 0)
  {
    if (!input->silent && input->display)
    {
    

      printf("\nQuery Answer:\n");
      for (i = 0; i < input->nq; i++)
      {
        printf("query %d: [ ", i);

        for (j = 0; j < input->nQA; j++)
          if (input->QA[j * 2] == i)
          {

            printf("%d ", input->QA[j * 2 + 1]);
          }
        printf("]\n");
      }
      printf("\n");
    }
    sprintf(fname, "%s.qa", dsname);
    save_data(fname, input->QA, input->nQA, 2);
   
  }

  if (!input->silent)
    printf("\nDone.\n");
  return 0;
}
