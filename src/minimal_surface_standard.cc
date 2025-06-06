#include <minimal_surface_standard.h>


namespace Minimal_Surface {

  double BoundaryValues::value(
    const Point<2> &p,
    const unsigned int /*component*/
    ) const {
    double a = 0.5;
    double out = 0.;
    if (p[0] >= 0. && p[0] < 0.25)
      out = a*p[0];
    else if (p[0] >= 0.25 && p[0] < 0.5)
      out = a*(0.5 - p[0]);
    else if(p[0] >= 0.5 && p[0] <= 0.75)
      out = a*(p[0] - 0.5);
    else if (p[0] > 0.75 && p[0] <= 1)
      out = a*(1 - p[0]);
    //out = 1.;       // inner radius = 1.54308
    //out = 1.55643;    // inner radius = 1.1
    return out;
  }

  double Minimal_RHS::value(
    const Point<2> &p,
    const unsigned int /*component*/
    ) const {
    return 0;
  }

  double Minimal_SOL::value(
    const Point<2> &p,
    const unsigned int /*component*/
    ) const {
    const double scaling = 1.0;

    double out = 0.;
    double radius = std::sqrt(p[0]*p[0] + p[1]*p[1]);
    out = -scaling * std::acosh(radius/scaling) + 2.;
    return out;
  }

  Tensor<1, 2> Minimal_SOL::gradient(
    const Point<2> &p,
    const unsigned int /*component*/
    ) const {
    Tensor<1, 2> out;
    double scaling = 1.0;
    double r = std::sqrt(p[0]*p[0] + p[1]*p[1]);
    out[0] = -1. * scaling*p[0] / (r* std::sqrt(r*r - scaling*scaling));
    out[1] = -1. * scaling*p[1] / (r* std::sqrt(r*r - scaling*scaling));
     
    return out;
  }


  Minimal_Benchmark_Standard::Minimal_Benchmark_Standard(
    int order
  )
    : dof_handler(triangulation)
    , fe(order)
    , rhs_fcn()
  {

    std::string name = "minimal_standard_square/";



    // Set the output file name for annulus
      const std::vector<std::string>& columns =
        {"Level", "Cycles", "Cells", "DoFs", "k_{Newton}", "update_norm", "last_norm"};
      const std::vector<std::string>& tex_captions = 
        {"Level", "Cycles", "# cells", "\\# dofs", "\\k_{Newton}", "update_norm", "last_norm"};
      const std::vector<bool>& scientific =
        {false, false, false, false, false, true, true};
      const std::vector<unsigned int> precision = 
        {0, 0, 0, 0, 0, 2, 2};
      const std::vector<std::string>& super_column_names = 
        {"Grid Info", "Newton"};
      const std::vector<std::vector<std::string>> super_columns =
        {{"Level", "Cycles", "Cells", "DoFs"}, {"k_{Newton}", "update_norm", "last_norm"}};
      
      // const std::vector<std::string>& columns =
      //   {"Level", "Cycles", "Cells", "DoFs", "k_{Newton}", "update_norm", "last_norm", "L2", "H1"};
      // const std::vector<std::string>& tex_captions = 
      //   {"Level", "Cycles", "# cells", "\\# dofs", "\\k_{Newton}", "update_norm", "last_norm", "$L_2$-error", "$H^1$-error"};
      // const std::vector<bool>& scientific =
      //   {false, false, false, false, false, true, true, true, true};
      // const std::vector<unsigned int> precision = 
      //   {0, 0, 0, 0, 0, 2, 2, 2, 2};
      // const std::vector<std::string>& super_column_names = 
      //   {"Grid Info", "Newton", "Errors"};
      // const std::vector<std::vector<std::string>> super_columns =
      //   {{"Level", "Cycles", "Cells", "DoFs"}, {"k_{Newton}", "update_norm", "last_norm"}, {"L2", "H1"}};



    problem_out = OutputSetup(  name
                                  , fe.degree
                                  , columns
                                  , tex_captions
                                  , scientific
                                  , precision
                                  , super_column_names
                                  , super_columns
                                  ); 
    problem_out.table.set_auto_fill_mode(true);
  }



  void Minimal_Benchmark_Standard::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    current_solution.reinit(dof_handler.n_dofs());
 
    zero_constraints.clear();
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<2>(),
                                             zero_constraints);
    // VectorTools::interpolate_boundary_values(dof_handler,
                                            //  1,
                                            //  Functions::ZeroFunction<2>(),
                                            //  zero_constraints);                                         
    DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);
    zero_constraints.close();
 
    nonzero_constraints.clear();

    // VectorTools::interpolate_boundary_values(dof_handler,
                                            //  1,
                                            //  Functions::ZeroFunction<2>(),
                                            //  nonzero_constraints);                                         
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             BoundaryValues(),
                                             nonzero_constraints);
    DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
    nonzero_constraints.close();
 
    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
 
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, zero_constraints);
 
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
  }
 
 
  
  void Minimal_Benchmark_Standard::assemble_system()
  {
    const QGauss<2> quadrature_formula(fe.degree + 1);
 
    system_matrix = 0;
    system_rhs    = 0;
 
    FEValues<2> fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points |
                              update_JxW_values);
 
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();
 
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
 
    std::vector<Tensor<1, 2>> old_solution_gradients(n_q_points);
 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        cell_rhs    = 0;
 
        fe_values.reinit(cell);
 
        fe_values.get_function_gradients(current_solution,
                                         old_solution_gradients);
 
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double coeff =
              1.0 / std::sqrt(1 + old_solution_gradients[q] *
                                    old_solution_gradients[q]);
 
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) +=
                    (((fe_values.shape_grad(i, q)      // ((\nabla \phi_i
                       * coeff                         //   * a_n
                       * fe_values.shape_grad(j, q))   //   * \nabla \phi_j)
                      -                                //  -
                      (fe_values.shape_grad(i, q)      //  (\nabla \phi_i
                       * coeff * coeff * coeff         //   * a_n^3
                       * (fe_values.shape_grad(j, q)   //   * (\nabla \phi_j
                          * old_solution_gradients[q]) //      * \nabla u_n)
                       * old_solution_gradients[q]))   //   * \nabla u_n)))
                     * fe_values.JxW(q));              // * dx
 
                cell_rhs(i) -= (fe_values.shape_grad(i, q)  // \nabla \phi_i
                                * coeff                     // * a_n
                                * old_solution_gradients[q] // * \nabla u_n
                                * fe_values.JxW(q));        // * dx
              }
          }
 
        cell->get_dof_indices(local_dof_indices);
        zero_constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
        
        
      }
      //system_matrix.print_pattern(std::cout);
  }
 
 
 
 
  
  void Minimal_Benchmark_Standard::solve()
  {
    SolverControl            solver_control(system_rhs.size(),
                                 system_rhs.l2_norm() * 1e-6);
    SolverCG<Vector<double>> solver(solver_control);
 
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);
 
    solver.solve(system_matrix, newton_update, system_rhs, preconditioner);
    std::cout << "        norm of newton_update = "
              << newton_update.l2_norm() << std::endl;
    zero_constraints.distribute(newton_update);
 
    const double alpha = determine_step_length();
    current_solution.add(alpha, newton_update);
  }
 
 
 
  
  void Minimal_Benchmark_Standard::refine_mesh()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 
    KellyErrorEstimator<2>::estimate(
      dof_handler,
      QGauss<2 - 1>(fe.degree + 1),
      std::map<types::boundary_id, const Function<2> *>(),
      current_solution,
      estimated_error_per_cell);
      std::cout << " cycle " << cycle << std::endl; 
    if (cycle < 4)
    {
      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.30,
                                                    0.0);
    }
    else
    {
      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.20,
                                                    0.0);
    }
    
 
    triangulation.prepare_coarsening_and_refinement();
 
    SolutionTransfer<2> solution_transfer(dof_handler);
    const Vector<double>  coarse_solution = current_solution;
    solution_transfer.prepare_for_coarsening_and_refinement(coarse_solution);
 
    triangulation.execute_coarsening_and_refinement();
 
    setup_system();
 
    solution_transfer.interpolate(coarse_solution, current_solution);
 
    nonzero_constraints.distribute(current_solution);
  }
 
 
 
 
  
  double Minimal_Benchmark_Standard::compute_residual(const double alpha) const
  {
    Vector<double> residual(dof_handler.n_dofs());
 
    Vector<double> evaluation_point(dof_handler.n_dofs());
    evaluation_point = current_solution;
    evaluation_point.add(alpha, newton_update);
 
    const QGauss<2> quadrature_formula(fe.degree + 1);
    FEValues<2>     fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points |
                              update_JxW_values);
 
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();
 
    Vector<double>              cell_residual(dofs_per_cell);
    std::vector<Tensor<1, 2>> gradients(n_q_points);
 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_residual = 0;
        fe_values.reinit(cell);
 
        fe_values.get_function_gradients(evaluation_point, gradients);
 
 
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double coeff =
              1. / std::sqrt(1 + gradients[q] * gradients[q]);
 
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_residual(i) -= (fe_values.shape_grad(i, q) // \nabla \phi_i
                                   * coeff                    // * a_n
                                   * gradients[q]             // * \nabla u_n
                                   * fe_values.JxW(q));       // * dx
          }
 
        cell->get_dof_indices(local_dof_indices);
        zero_constraints.distribute_local_to_global(cell_residual,
                                                    local_dof_indices,
                                                    residual);
      }
 
    return residual.l2_norm();
  }
 
 
  
  double Minimal_Benchmark_Standard::determine_step_length() const
  {
    return 0.1;
  }
 
 
 
 
  
  void Minimal_Benchmark_Standard::output_results(
    const unsigned int refinement_cycle) const
  {
    DataOut<2> data_out;
 
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(current_solution, "solution");
    data_out.add_data_vector(newton_update, "update");
    data_out.build_patches();
 
    const std::string filename =
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);
  }
 



  void Minimal_Benchmark_Standard::print_mesh_file(
  ) const {
    FullMatrix<double> cell_list(triangulation.n_active_cells(), 4);
    FullMatrix<double> vert_list(triangulation.n_used_vertices(), 2);
    unsigned int ind = 0;
    for (const auto& cell : triangulation.active_cell_iterators()){
      for (unsigned int v = 0; v < 4; v++){
        cell_list(ind, v) = cell -> vertex_index(v);
        const Point<2>& vertex = cell -> vertex(v);
        vert_list(cell->vertex_index(v), 0) = vertex(0); 
        vert_list(cell->vertex_index(v), 1) = vertex(1); 
      }
      ind++;
    }

    const std::string mesh = problem_out.degree.string() 
                              + "l"
                              + std::to_string(triangulation.n_levels() - 1)
                              + "_mesh.txt";
    std::ofstream mesh_out(mesh, std::ios::out | std::ios::trunc);
    mesh_out << 4 
             << " " 
             << triangulation.n_active_cells() 
             << " " 
             << triangulation.n_used_vertices() 
             << " " 
             << 0 
             << " " 
             << 2 
             << std::endl;

    for (unsigned int i = 0; i < triangulation.n_active_cells(); i++){
      for (unsigned int d = 0; d < 4; d++)
        mesh_out << cell_list(i, d) << " ";

      mesh_out << std::endl;
    }

    for (unsigned int i = 0; i < triangulation.n_used_vertices(); i++){
      for (unsigned int d = 0; d < 2; d++)
        mesh_out << vert_list(i, d) << " ";
      
      mesh_out << std::endl;
    } 

    mesh_out.close();

  } // print_mesh_file

  void Minimal_Benchmark_Standard::process_results(
    const unsigned int cycle
  ) {
  	//.vtu for visualization
    std::cout << "        Printing to files" << std::endl;


    // Minimal_SOL sol_fcn;
    // Vector<float> difference_per_cell(triangulation.n_active_cells());
    // VectorTools::integrate_difference(dof_handler,
    //                                   current_solution,
    //                                   sol_fcn,
    //                                   difference_per_cell,
    //                                   QGauss<2>(fe.degree + 1),
    //                                   VectorTools::L2_norm);
    // const double L2_error =
    //   VectorTools::compute_global_error(triangulation,
    //                                     difference_per_cell,
    //                                     VectorTools::L2_norm);
 
    // VectorTools::integrate_difference(dof_handler,
    //                                   current_solution,
    //                                   sol_fcn,
    //                                   difference_per_cell,
    //                                   QGauss<2>(fe.degree + 1),
    //                                   VectorTools::H1_seminorm);
    // const double H1_semi_error =
    //   VectorTools::compute_global_error(triangulation,
    //                                     difference_per_cell,
    //                                     VectorTools::H1_seminorm);

    // const double H1_error = 
    //   std::sqrt(L2_error * L2_error + H1_semi_error * H1_semi_error);


  	Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

  	KellyErrorEstimator<2>::estimate(
                     // map,
                     dof_handler,
  	 								 QGauss<2 - 1>(2 * fe.degree + 1),
  									 std::map<types::boundary_id, const Function<2> *>(),
  									 current_solution,
  									 estimated_error_per_cell);

  	DataOut<2> data_out;
  	data_out.attach_dof_handler(dof_handler);
  	data_out.add_data_vector(current_solution, "solution");
    data_out.add_data_vector(newton_update, "update");
    data_out.add_data_vector(estimated_error_per_cell, "estimated_error");
  	data_out.build_patches();

    

    const std::string folder_vtg = problem_out.vtg.string();
    const std::string folder_svg = problem_out.svg.string();
    const std::string level      = std::to_string(triangulation.n_levels() - 1);

    std::string vtu_name = folder_vtg
                          + "/step0_standard_grid_l"
                          + level
                          + ".vtu";
    std::string svg_name = folder_svg
                          + "/step0_standard_grid_l"
                          + level 
                          + ".svg";


    GridOutFlags::Svg svg_flags;
    svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
    // svg_flags.label_level_number  = true;
    // svg_flags.label_cell_index    = true;
    // svg_flags.label_boundary_id   = true;
  
    std::ofstream vtu_out(vtu_name);

  	data_out.write_vtu(vtu_out);
    GridOut       grid_out;
    grid_out.set_flags(svg_flags);
  
    if (triangulation.n_levels()-1 < 15) {
  	  std::ofstream svg_out(svg_name);
      grid_out.write_svg(triangulation, svg_out);
  	  //print_numerical_solution();
      print_mesh_file();
    }


  	const unsigned int n_active_cells = triangulation.n_active_cells();
  	const unsigned int n_dofs         = dof_handler.n_dofs();
  
  	// output in console
  	std::cout << "Refinement " << cycle << ':' << std::endl
  	<< "   Number of active cells:       " << n_active_cells << std::endl
  	<< "   Number of degrees of freedom: " << n_dofs << std::endl;
  
  	// cast unsigned int to int
    if ( triangulation.n_levels() != old_level ) {

  	problem_out.add_single_value_to_table("Level", triangulation.n_levels() - 1);
    problem_out.add_single_value_to_table("Cycles", cycle);
    problem_out.add_single_value_to_table("Cells", n_active_cells);
    problem_out.add_single_value_to_table("DoFs", n_dofs);
    problem_out.add_single_value_to_table("k_{Newton}", newton_iteration);
    problem_out.add_single_value_to_table("update_norm", newton_update.l2_norm());
    problem_out.add_single_value_to_table("last_norm", estimated_error_per_cell.l2_norm());
    // problem_out.add_single_value_to_table("L2", L2_error);
    // problem_out.add_single_value_to_table("H1", H1_error);


  	// output in .text 
    problem_out.write_table_text();
    problem_out.write_table_tex();
    std::cout << "        Outputting table level " << triangulation.n_levels() - 1 << std::endl;
    old_level = triangulation.n_levels();
    }
  } // process_results

 


  void Minimal_Benchmark_Standard::run(
    const unsigned int ref
  )
  {
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(1);

    // const Point<2> center(0., 0.);
    // const double inner_radius = 1.1, outer_radius = 3.762196;
    // //const double inner_radius = 1.54308, outer_radius = 3.762196;
    // GridGenerator::hyper_shell(
    //   triangulation,
    //   center,
    //   inner_radius,
    //   outer_radius,
    //   4,
    //   true
    // );
    // triangulation.refine_global(1);

    std::cout << "    Setting system... " << std::endl;
  	setup_system();
    nonzero_constraints.distribute(current_solution);

    //refinement cycle
    while (triangulation.n_levels() < ref + 1 )
    {
  		std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle != 0)
      {
        std::cout << "    Refining grid..." << std::endl;
        refine_mesh();
      }

      std::cout << "  Initial residual: " << compute_residual(0) << std::endl;
      newton_iteration = 0;
      while (newton_iteration < 50)
      {
        assemble_system();
        solve();
        double current_residual = compute_residual(0);
        std::cout << "  Residual in step " << newton_iteration <<": " << current_residual << std::endl;
        if (current_residual < 1e-8)
          break;
        newton_iteration++;
      }
      std::cout << "    Processing results..." << std::endl;
  	  process_results(cycle++);
  	}
    

    //problem_out.write_table_text(std::cout);
	}

}