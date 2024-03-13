#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Tpetra_Version.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

typedef void *laplacian_t;
extern "C" laplacian_t laplacian_weighted(long long *vl, unsigned nel,
                                          unsigned nv, MPI_Comm comm,
                                          int verbose);
extern "C" uint laplacian_size(laplacian_t L);
extern "C" ulong *laplacian_rows(laplacian_t L);
extern "C" ulong *laplacian_columns(laplacian_t L);
extern "C" double *laplacian_values(laplacian_t L);
extern "C" void laplacian_free(laplacian_t *L);

using Teuchos::ArrayView;
using Teuchos::Comm;
using Teuchos::MpiComm;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;

using local_t = int;
using global_t = long long;
using scalar_t = double;

using Map_t = Tpetra::Map<local_t, global_t>;
using Vector_t = Tpetra::Vector<scalar_t, local_t, global_t>;
using Matrix_t = Tpetra::CrsMatrix<scalar_t, local_t, global_t>;

using MatrixAdapter_t = Zoltan2::XpetraCrsMatrixAdapter<Matrix_t>;
using MultiVectorAdapter_t = Zoltan2::XpetraMultiVectorAdapter<Vector_t>;
using PartitioningProblem = Zoltan2::PartitioningProblem<MatrixAdapter_t>;
using PartitioningSolution = Zoltan2::PartitioningSolution<MatrixAdapter_t>;
using part_t = MatrixAdapter_t::part_t;

static void check_solution(const RCP<Matrix_t> &matrix,
                           const MatrixAdapter_t &adapter,
                           const PartitioningSolution &solution, const int rank,
                           const int verbose) {
  RCP<Vector_t> product, initial;
  {
    product = Tpetra::createVector<scalar_t, local_t, global_t>(
        matrix->getRangeMap());
    initial = Tpetra::createVector<scalar_t, local_t, global_t>(
        matrix->getDomainMap());
    initial->randomize();
    matrix->apply(*initial, *product);
  }

  if (rank == 0 && verbose) {
    fprintf(stderr, "Redistributing matrix and the vector ...\n");
    fflush(stdout);
  }

  RCP<Matrix_t> matrix2;
  adapter.applyPartitioningSolution(*matrix, matrix2, solution);

  RCP<Vector_t> product2, initial2;
  {
    MultiVectorAdapter_t vector_adapter(initial);
    vector_adapter.applyPartitioningSolution(*initial, initial2, solution);
    product2 = Tpetra::createVector<scalar_t, local_t, global_t>(
        matrix2->getRangeMap());
    matrix2->apply(*initial2, *product2);
  }

  scalar_t norm = product->norm2();
  scalar_t norm2 = product2->norm2();

  if (rank == 0 && verbose) {
    fprintf(stderr, "norm = %e, norm2 = %e\n", norm, norm2);
    fflush(stdout);
  }
}

extern "C" int Zoltan2_partMesh(int *part, long long *vl, unsigned nel, int nv,
                                MPI_Comm comm_, int verbose) {
  int rank, size;
  {
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &size);
  }

  const long long nel_ = nel;
  long long num_global_elements = 0, element_offset;
  double imbalance_tol = 0;
  {
    MPI_Allreduce(&nel_, &num_global_elements, 1, MPI_LONG_LONG, MPI_SUM,
                  comm_);
    MPI_Scan(&nel_, &element_offset, 1, MPI_LONG_LONG, MPI_SUM, comm_);

    long long num_elements_per_rank = num_global_elements / size;
    imbalance_tol = (num_elements_per_rank + 1.0) / num_elements_per_rank;

    if (rank == 0 && verbose) {
      fprintf(stderr,
              "num_global_elements = %lld element_offset = %lld "
              "imbalance_tol = %lf\n",
              num_global_elements, element_offset, imbalance_tol);
      fflush(stdout);
    }
  }

  RCP<Map_t> map;
  {
    RCP<const Comm<int>> comm(new MpiComm<int>(comm_));
    map = rcp(new Map_t(num_global_elements, nel, 0, comm));
    assert(map->isContiguous());
  }

  RCP<Matrix_t> matrix;
  {
    laplacian_t L = laplacian_weighted(vl, nel, nv, comm_, verbose);
    local_t size = laplacian_size(L);
    ulong *rows = laplacian_rows(L);
    ulong *columns = laplacian_columns(L);
    double *values = laplacian_values(L);

    matrix = rcp(new Matrix_t(map, 27));
    for (local_t id = 0; id < size; id++) {
      global_t row = rows[id];
      global_t column = columns[id];
      scalar_t value = values[id];
      matrix->insertGlobalValues(row, ArrayView<global_t>(&column, 1),
                                 ArrayView<scalar_t>(&value, 1));
    }

    laplacian_free(&L);
  }

  matrix->fillComplete();

  if (rank == 0 && verbose) {
    fprintf(stderr,
            "global: num_rows = %lld, num_non_zeros = %lld, num_procs = %lld\n",
            matrix->getGlobalNumRows(), matrix->getGlobalNumEntries(), size);
    fprintf(stderr, "local: num_rows = %lld, num_non_zeros = %lld, nel = %d\n",
            matrix->getLocalNumRows(), matrix->getLocalNumEntries(), nel);
    fflush(stderr);
  }

  ParameterList params;
  {
    params.set("partitioning_approach", "partition");
    params.set("algorithm", "parmetis");
    params.set("imbalance_tolerance", imbalance_tol);
    params.set("num_global_parts", size);
  }

  MatrixAdapter_t adapter(matrix);
  PartitioningProblem problem(&adapter, &params);

  try {
    problem.solve();
  } catch (std::exception &e) {
    fprintf(stderr, "Exception returned from solve() : %s\n", e.what());
    fflush(stderr);
    return 1;
  }

  PartitioningSolution solution = problem.getSolution();

  const part_t *assignments = solution.getPartListView();
  for (local_t id = 0; id < nel; id++) part[id] = assignments[id];

  return 0;
}
