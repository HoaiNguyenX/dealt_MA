/*
 * minimal_surface.h
 *
 *  Created on: 01.03.2025
 *      Author: nguyen
 */


#ifndef INC_MINIMAL_SURFACE_H_
#define INC_MINIMAL_SURFACE_H_

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


namespace Minimal_Surface {

  using namespace dealt;
  using namespace dealii;

  class Minimal_DC_circle : public Function<2> {
  public:
    Minimal_DC_circle() : Function<2>(1) {}
    ~Minimal_DC_circle() = default;

    virtual double value(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;
  };

  class Minimal_SOL_catenoid : public Function<2> {
  public:
    Minimal_SOL_catenoid() : Function<2>(1) {}
    ~Minimal_SOL_catenoid() = default;

    virtual double value(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;

    virtual Tensor<1, 2> gradient(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;
  };


  class Minimal_DC_catenoid : public Function<2> {
  public:
    Minimal_DC_catenoid() : Function<2>(1) {}
    ~Minimal_DC_catenoid() = default;

    virtual double value(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;
  };

  class Minimal_DC_square : public Function<2> {
  public:
    Minimal_DC_square() : Function<2>(1) {}
    ~Minimal_DC_square() = default;

    virtual double value(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;
  };

  enum RefinementStrategy{
    Adaptive, 
    Uniform
  };

  enum ProblemShape{
    Halfcircle,
    Tilted_Halfcircle,
    Annulus,
    Square
  };

  enum StepLengthStrategy{
    Constant,
    LineSearch
  };

  class Minimal_Benchmark 
  {
    enum Boundary {
      None,
      Dirichlet_0,
      Dirichlet_c,
      Dirichlet_a,
      Dirichlet_s
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

    Minimal_SOL_catenoid      cat_sol_fcn;

    Minimal_DC_circle         circ_DC_fcn;
    Minimal_DC_square         square_DC_fcn;   
    Minimal_DC_catenoid       cat_DC_fcn;
    
    ProblemShape              problem_shape;
    RefinementStrategy        refinement_strategy;
    StepLengthStrategy        step_length_strategy;
  public:
    Minimal_Benchmark(
      int ref,
      int order,
      const ProblemShape shape,
      const RefinementStrategy strategy = RefinementStrategy::Uniform,
      const StepLengthStrategy steplength_strategy = StepLengthStrategy::Constant
    );
    void run();

    
  private:
    IPF_Data<2, 2> get_IPF_data(ProblemShape shape);
    IPF_Data<2, 2> get_IPF_data_halfcircle();
    IPF_Data<2, 2> get_IPF_data_tiltedhalfcircle();
    IPF_Data<2, 2> get_IPF_data_annulus();
    IPF_Data<2, 2> get_IPF_data_square();

    void   setup_system();
    void   assemble_system();
    void   impose_boundary_condition(ProblemShape shape);
    double compute_residual_for_steplength(double alpha);
    void   solve_system();
    void   output_system();
    double determine_step_length_const() const {return 1.;};
    double determine_step_length_LS() const;
    void   estimate_and_mark();   
    void   print_numerical_solution(std::string addition = "");

  };

} // namespace Minimal_Surface


#endif /* INC_MINIMAL_SURFACE_H_ */
  