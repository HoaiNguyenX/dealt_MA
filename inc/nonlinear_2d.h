/*
 * nonlinear.h
 *
 *  Created on: 07.04.2025
 *      Author: nguyen
 */


#ifndef INC_NONLINEAR_2D_H_
#define INC_NONLINEAR_2D_H_

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

#include <fstream>
#include <iostream>

#include <utility>
#include <chrono>


namespace Nonlinear {

  using namespace dealii;
  using namespace dealt;

  class Nonlinear2D_RHS : public Function<2> {
  public:
    Nonlinear2D_RHS() : Function<2>(1) {}
    ~Nonlinear2D_RHS() = default;

    virtual double value(
      const Point<2>&       p,
      const unsigned int    component = 0
    ) const override;
        
    virtual Tensor<1, 2> gradient(
      const Point<2>&       p,
      const unsigned int    component = 0
    ) const override; 
  };

  class Nonlinear2D_NC : public Function<2> {
  public:
    Nonlinear2D_NC() : Function<2>(1) {}
    ~Nonlinear2D_NC() = default;
    virtual double value(
      const Point<2>& p,
      const unsigned int     component = 0
    ) const override;
  };

  enum ProblemCase{
    Case_1,
    Case_2
  };
  
  enum RefinementStrategy{
    Adaptive, 
    Uniform
  };

  enum StepLengthStrategy{
    Constant,
    LineSearch
  };


  class Nonlinear2D_Benchmark
  {
    enum Boundary {
      None, 
      Neumann_Case_1,
      Dirichlet_0,
      Dirichlet
    };

  private:
    int ref;
    int order;
    int exponent; // exponent k of equation. k=3

    IPF_Data<2, 2>            data;
    TS_Triangulation<2, 2>    tria;


    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;

    OutputSetup               problem_out;

    Vector<double>            system_rhs;
    Vector<double>            current_solution;
    Vector<double>            newton_update;

    unsigned int              cycle;

    Nonlinear2D_RHS             rhs_fcn;
    Nonlinear2D_NC              nc_fcn;

    RefinementStrategy        refinement_strategy;
    StepLengthStrategy        step_length_strategy;


  public:
    Nonlinear2D_Benchmark(
      int ref,
      int order,
      int exponent,
      const ProblemCase problem_case,
      const RefinementStrategy strategy = RefinementStrategy::Uniform,
      const StepLengthStrategy steplength_strategy = StepLengthStrategy::Constant
    );
    
    void run();

  private:
    IPF_Data<2> get_IPF_data() const;


    void   setup_system();
    void   assemble_system();
    void   impose_boundary_condition();
    double compute_residual_for_steplength(double alpha) const;
    void   solve_system();
    void   output_system();
    double determine_step_length_const() const {return 0.1;};
    double determine_step_length_LS() const;
    void   estimate_and_mark();   
    void   print_numerical_solution(std::string addition = "");

  };

} // namespace Nonlinear

#endif /* INC_NONLINEAR_2D_H_ */