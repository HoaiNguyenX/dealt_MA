



#ifndef INC_NONLINEAR_3D_H_
#define INC_NONLINEAR_3D_H_


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

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>

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

#include <nonlinear_2d.h>



namespace Nonlinear {

  using namespace dealii;
  using namespace dealt;

  class Nonlinear3D_RHS1 : public Function<3> {
  public:
    Nonlinear3D_RHS1() : Function<3>(1) {}
    ~Nonlinear3D_RHS1() = default;

    virtual double value(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override;
        
    virtual Tensor<1, 3> gradient(
      const Point<3>&       p,
      const unsigned int    component = 0
    ) const override; 
  };

  class Nonlinear3D_NC1 : public Function<3> {
  public:
    Nonlinear3D_NC1() : Function<3>(1) {}
    ~Nonlinear3D_NC1() = default;
    virtual double value(
      const Point<3>& p,
      const unsigned int     component = 0
    ) const override;
  };

  class Nonlinear3D_NC2 : public Function<3> {
  public:
    Nonlinear3D_NC2() : Function<3>(1) {}
    ~Nonlinear3D_NC2() = default;
    virtual double value(
      const Point<3>& p,
      const unsigned int     component = 0
    ) const override;
  };



  class Nonlinear3D_Benchmark
  {
    enum Boundary {
      None, 
      Dirichlet_0,
      Neumann_Case_1,
      Neumann_Case_2,
      Dirichlet
    };

  private:
    int ref;
    int order;
    int exponent; // exponent k of equation. k=3

    IPF_Data<3, 3>            data;
    TS_Triangulation<3, 3>    tria;


    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;

    OutputSetup               problem_out;

    Vector<double>            system_rhs;
    Vector<double>            current_solution;
    Vector<double>            newton_update;

    unsigned int              cycle;

    double                    L2 = 0;
    double                    H1 = 0;


    Function<3>*              rhs_fcn;
    Nonlinear3D_NC1           nc1_fcn;
    Nonlinear3D_NC2           nc2_fcn;

    ProblemCase               problem_case;
    RefinementStrategy        refinement_strategy;
    StepLengthStrategy        step_length_strategy;


  public:
    Nonlinear3D_Benchmark(
      int ref,
      int order,
      int exponent,
      const ProblemCase problem_case,
      const RefinementStrategy strategy = RefinementStrategy::Uniform,
      const StepLengthStrategy steplength_strategy = StepLengthStrategy::Constant
    );
    
    void run();

  private:
    IPF_Data<3> get_IPF_data() const;


    void   setup_system();
    void   assemble_system();
    void   impose_boundary_condition();
    double compute_residual_for_steplength(double alpha) const;
    void   solve_system();
    void   output_system();
    double determine_step_length_const() const {return 1;};
    double determine_step_length_LS() const;
    void   estimate_and_mark();   
    void   print_numerical_solution(std::string addition = "");
    void   get_neumann_data(
      std::map< types::boundary_id,
                      const Function<3>* >& neumann_data
      );
  };

} // namespace Nonlinear

#endif /* INC_NONLINEAR_3D_H_ */