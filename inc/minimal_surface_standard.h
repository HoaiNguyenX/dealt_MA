/*
 * minimal_surface_standard.h
 * A copy of step-15 from deal.ii for comparison
 *
 *  Created on: 28.02.2025
 *      Author: nguyen
 */

#ifndef INC_MINIMAL_SURFACE_STANDARD_H_
#define INC_MINIMAL_SURFACE_STANDARD_H_

#include <memory>
#include <utility>
#include <chrono>

#include <utilities.h>
#include <ts_triangulation.h>
#include <tspline_function.h>
#include <isoparametric_function.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/convergence_table.h>
 
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
 
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
 
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
 
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
 
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
 
 
#include <fstream>
#include <iostream>
 
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include <cmath>
#include <array>
#include <fstream>
#include <iostream>
#include <filesystem>


namespace Minimal_Surface {
  using namespace dealii;
  using namespace dealt;

  class BoundaryValues : public Function<2>
  {
  public:
    virtual double value(
      const Point<2>  &p,
      const unsigned int component = 0
    ) const override;
  };
  
  class Minimal_RHS : public Function<2> {
  public:
    Minimal_RHS() : Function<2>(1) {}
    ~Minimal_RHS() = default;

    virtual double value(
      const Point<2>&    p,
      const unsigned int component = 0
    ) const override;
  };

  class Minimal_SOL : public Function<2> {
  public:
    Minimal_SOL() : Function<2>(1) {}
    ~Minimal_SOL() = default;

    virtual double value(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;

    virtual Tensor<1, 2> gradient(
      const Point<2>&     p,
      const unsigned int  component = 0
    ) const override;
  };

  class Minimal_Benchmark_Standard {

  private:
    Triangulation<2> triangulation;
 
    DoFHandler<2> dof_handler;
    const FE_Q<2> fe;
 
    AffineConstraints<double> zero_constraints;
    AffineConstraints<double> nonzero_constraints;
 
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
 
    Vector<double> current_solution;
    Vector<double> newton_update;
    Vector<double> system_rhs;

    Minimal_RHS    rhs_fcn;

    double H1 = 1;
    double L2 = 1;


    int cycle = 0;
    int newton_iteration = 0;
    int old_level = 0;

    OutputSetup problem_out;


  public:
    Minimal_Benchmark_Standard(int order);
    void run(
      const unsigned int ref
    );

  private:
    void   setup_system();
    void   assemble_system();
    void   solve();
    void   refine_mesh();
    double compute_residual(const double alpha) const;
    double determine_step_length() const;
    void   output_results(const unsigned int refinement_cycle) const;
    void   print_mesh_file() const;
    void compute_h1_error();
    void process_results(const unsigned int cycle);

  };

} // namespace Minimal_Surface_Standard

#endif 