#ifndef INC_MINIMAL_TEST_H_
#define INC_MINIMAL_TEST_H_

#include <memory>

#include <utilities.h>
#include <ts_triangulation.h>
#include <tspline_function.h>
#include <isoparametric_function.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/solver_relaxation.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/lac/sparse_direct.h>

#include <fstream>
#include <iostream>

#include <utility>
#include <chrono>

#include <minimal_surface.h>


namespace Minimal_Surface {
  using namespace dealt;
  using namespace dealii;


  class Minimal_RHS1 : public Function<2> {
  public:
    Minimal_RHS1() : Function<2>(1) {}
    ~Minimal_RHS1() = default;

    virtual double value(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;
  };
  class Minimal_NC1 : public Function<2> {
  public:
    Minimal_NC1() : Function<2>(1) {}
    ~Minimal_NC1() = default;

    virtual double value(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;
  };

class Minimal_DC : public Function<2> {
  public:
    Minimal_DC() : Function<2>(1) {}
    ~Minimal_DC() = default;

    virtual double value(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;
  };
  class Minimal_Test 
  {
    enum Boundary {
      None,
      Dirichlet_0,
      Dirichlet_1,
      Dirichlet_s,
      Neumann_Case_1

    };
    
  private:
    int ref;
    int order;

    IPF_Data<2, 2>            data;
    TS_Triangulation<2, 2>    tria;


    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;

    OutputSetup               problem_out;

    Vector<double>            system_rhs;
    Vector<double>            current_solution;
    Vector<double>            newton_update;

    unsigned int              cycle;

    Minimal_DC                DC_fcn;
    Minimal_NC1               nc_fcn;
    Minimal_RHS1               rhs_fcn;


    ProblemShape              problem_shape;
    RefinementStrategy        refinement_strategy;
  public:
    Minimal_Test(
      int ref,
      int order,
      const ProblemShape shape,
      const RefinementStrategy strategy = RefinementStrategy::Uniform
    );
    void run();

    
  private:
    IPF_Data<2, 2> get_IPF_data();


    void   setup_system();
    void   assemble_system();
    void   impose_boundary_condition();
    double compute_residual_for_steplength(double alpha);
    void   solve_system();

    void   output_system();
    double determine_step_length_const() const {return 1.;};
    void   estimate_and_mark();   
    void   print_numerical_solution(std::string addition = "");
  };
}




#endif /* INC_MINIMAL_TEST_H_ */