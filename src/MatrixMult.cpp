#include "GraphMatRuntime.cpp"



template <class E>
class MatrixMult : public GraphProgram<double, double, double, E> {
//T::message_type, U::message_reduction_type, V::vertex_property_type, E::edge_type{
  public:

  MatrixMult() {
   this->order = OUT_EDGES;
   this->activity = ALL_VERTICES;
  }

  void reduce_function(double& a, const double& b) const {
    a += b;     
  }
  void process_message(const double& message, const E edge_val, const double& vertexprop, double& res) const {
    res = edge_val * message;
    //printf ("message %lf edge_val %d vertexprop %f, res %f \n", message, edge_val,vertexprop,res);

  }
  bool send_message(const double& vertexprop, double& message) const {
      message = vertexprop;
      return true;
  }
  void apply(const double& message_out, double& vertexprop) {
    //double v=vertexprop;
    vertexprop = message_out; 
    //printf("apply %lf -> %lf \n", v, vertexprop);
  }

};


template <class edge>
void matrix_col_mult(Graph<double, edge> &A , SparseVector<edge>& b_row) {
  for (int i = 0; i <b_row.length; i++) { 
    if(b_row.exists(i)){
      A.setVertexproperty(i+1,(double)b_row.getValue(i));
   }else{
    A.setVertexproperty(i+1,(double)0);   
   }
  }
  MatrixMult<edge> mm;
  auto dg_tmp = graph_program_init(mm, A);
 
//  struct timeval start, end;
//  gettimeofday(&start, 0);

  A.setAllActive();
 
  run_graph_program(&mm, A, 1, &dg_tmp);

//  gettimeofday(&end, 0);
//  double time = (end.tv_sec-start.tv_sec)*1e3+(end.tv_usec-start.tv_usec)*1e-3;
//  printf("Time  = %.3f ms \n", time);

  graph_program_clear(dg_tmp);
#if 0
  for (int i = 1; i <= (unsigned long long int)A.getNumberOfVertices(); i++) { 
    printf("%d : %lf \n", i, A.getVertexproperty(i));
  }
#endif

}

template <class edge>
void run_matrix_mult(const char* filename_A, const char* filename_B,int nthreads) {

  Graph<double, edge> A;
  Graph<double, edge> B;

  struct timeval start, end;
  gettimeofday(&start, 0);
  
  A.ReadMTX(filename_A, nthreads*4); //nthread pieces of matrix
  B.ReadMTX(filename_B, nthreads*4); //nthread pieces of matrix
  

 for(int j=0;j<B.nparts;j++){ 
 MatrixDC<edge>* Bm = B.matT[j];
  for(int i=0; i<Bm->nzx ;i++){
    SparseVector<edge>  r =Bm->getRow(i);
    matrix_col_mult(A,r );
  }
}
 
  gettimeofday(&end, 0);
  double time = (end.tv_sec-start.tv_sec)*1e3+(end.tv_usec-start.tv_usec)*1e-3;
  printf("Time  = %.3f ms \n", time);

 return; 
}

int main(int argc, char* argv[]) {

  const char* input_filename_A = argv[1];
  const char* input_filename_B = argv[2];
  if (argc < 3) {
    printf("Correct format: %s A.mtx B.mtx\n", argv[0]);
    return 0;
  }

#pragma omp parallel
  {
#pragma omp single
    {
      nthreads = omp_get_num_threads();
      printf("num threads got: %d\n", nthreads);
    }
  }
  
  run_matrix_mult<int>(input_filename_A, input_filename_B, nthreads);

  
}

