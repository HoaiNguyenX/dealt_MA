#include <nonlinear_2d.h>


namespace Nonlinear {

  Nonlinear2D_Benchmark::Nonlinear2D_Benchmark(
    int ref,
    int order,
    int exponent,
    const ProblemCase problem_case,
    const RefinementStrategy strategy,
    const StepLengthStrategy steplength
  ) : ref(ref)
    , order(order)
    , exponent(exponent)
    , data(this->get_IPF_data())
    , tria(data)
    , cycle(0)
    , refinement_strategy(strategy)
    , problem_case(problem_case)
    , step_length_strategy(steplength)
  {
    tria.degree_elevate_global(order);
    //tria.refine_global(2);
    tria.refine_bezier_elements();
    const auto& cell = tria.begin_active();
    const auto& CPs = tria.get_control_points(cell);
    tria.coarsen_bezier_elements();

    tria.prepare_assembly();
    

    std::string name = "nonlinear2d";
    if (strategy == RefinementStrategy::Adaptive)
      name += "_adaptive/";
    else
      name += "_uniform/";
    
    name += "k" + std::to_string(exponent) + "_";

    if (problem_case == ProblemCase::Case_1)
    {
      name += "case_1/";
      rhs_fcn = new Nonlinear2D_RHS1();
    }
    else if (problem_case == ProblemCase::Case_2)
    {
      name += "case_2/"; 
      rhs_fcn = new Nonlinear2D_RHS2();
    }

    else if (problem_case == ProblemCase::Case_3)
    {
      name += "case_3/";
      rhs_fcn = new Nonlinear2D_RHS3();
    }

    
    else
      AssertThrow(false, ExcNotImplemented());


    const std::vector<std::string>& columns =
      {"Level", "Cycles", "Cells", "DoFs", "k_{Newton}", "update_norm", "initial_norm", "last_norm"};
    const std::vector<std::string>& tex_captions = 
      {"Level", "Cycles", "\\# cells", "\\# dofs", "k_{Newton}", "update_norm", "initial_norm", "last_norm"};
    const std::vector<bool>& scientific =
      {false, false, true, true, true, true, true, true};
    const std::vector<unsigned int> precision = 
      {0, 0, 0, 0, 0, 2, 2, 2};
    const std::vector<std::string>& super_column_names = 
      {"Grid Info", "Newton"};
    const std::vector<std::vector<std::string>> super_columns =
      {{"Level", "Cycles", "Cells", "DoFs"}, {"k_{Newton}", "update_norm", "initial_norm", "last_norm"}};

    problem_out = OutputSetup(  name
                              , data.max_degree() + order
                              , columns
                              , tex_captions
                              , scientific
                              , precision
                              , super_column_names
                              , super_columns
                              ); 
    problem_out.table.set_auto_fill_mode(true);




    
    // Set boundary indicators
    for (auto& face : tria.active_face_iterators())
    { 
      std::cout << "face: " << face -> index() << " with coords ";
      std::cout << face->center()[0] << ", "<< face->center()[1]<<std::endl;
      if (!face -> at_boundary())
        continue;
      else if (problem_case == ProblemCase::Case_1)
      {
        const Point<2>& c = face -> center();
        if (std::fabs(c(0)) < 1e-15)
          face -> set_boundary_id(Boundary::Neumann_Case_1);
        else
          face -> set_boundary_id(Boundary::Dirichlet_0);
      }
      else if (problem_case == ProblemCase::Case_2)
      {
        const Point<2>& c = face -> center();
        if (std::fabs(c(0)) < 1e-15)
          face -> set_boundary_id(Boundary::Neumann_Case_2);
        else
          face -> set_boundary_id(Boundary::Dirichlet_0);
      }
      else if (problem_case == ProblemCase::Case_3)
      {
        const Point<2>& c = face -> center();
        if (std::fabs(c(0)) < 1e-15)
          face -> set_boundary_id(Boundary::Neumann_Case_3);
        else
          face -> set_boundary_id(Boundary::Dirichlet_0);
      }
      else
        face -> set_boundary_id(Boundary::None);
    } // for ( face )
  } // constructor

    void Nonlinear2D_Benchmark::setup_system()
  {
    std::cout << "Setting up system ... " << std::endl;

    unsigned int n_global_dofs = tria.n_active_splines();
    system_rhs.reinit(n_global_dofs);
    current_solution.reinit(n_global_dofs);
    newton_update.reinit(n_global_dofs);   
  
    const auto& IEN_array = tria.get_IEN_array();
    Assert(IEN_array.size() > 0, ExcInternalError());

    sparsity_pattern.reinit(
        n_global_dofs,
        n_global_dofs,
        n_global_dofs );
  
    for (const auto& [_, arr] : IEN_array)
      for (unsigned int i : arr)
        for (unsigned int j : arr)
          sparsity_pattern.add(i, j);
  
    // Free superfluous space
    sparsity_pattern.compress();
    system_matrix.reinit(sparsity_pattern);
  } // setup_system

  void Nonlinear2D_Benchmark::assemble_system()
  {

    system_matrix = 0;
    system_rhs    = 0;

    // Setup initial tables that store the bernstein values / grads / hessians.
  
    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = degrees[0]  + 2;
    degrees[1] = degrees[1]  + 2;

    TSValues<2> ts_values(
        &tria,
        degrees,
        update_values |
        update_gradients |
        update_quadrature_points |
        update_JxW_values);

    TSFaceValues<2> face_values(
        &tria,
        degrees,
        update_values |
        update_quadrature_points |
        update_JxW_values);

    const unsigned int dofs_per_cell = ts_values.n_dofs_per_cell();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);


    const int k = exponent;

    for (const auto& cell : tria.active_cell_iterators())
    {
      // Get the Bernstein values on the cell
      ts_values.reinit(cell);

      // Get the values of the old solution for the current cell
      std::vector< double > old_value(ts_values.n_quadrature_points_per_cell());
      ts_values.get_function_values(current_solution, old_value);

      // Get the gradients of the old solution for the current cell
      std::vector< Tensor<1, 2> > old_gradient(ts_values.n_quadrature_points_per_cell());
      ts_values.get_function_gradients(current_solution, old_gradient);

      // Reset the cell matrix
      cell_matrix       = 0;
      cell_rhs          = 0;

      // Quadrature sum:
      for (const unsigned int q : ts_values.quadrature_point_indices())
      {
        const Point<2>& mapped_q = ts_values.quadrature_point(q);
        const double rhs = rhs_fcn -> value(mapped_q);


        // Build the cell matrix and rhs
        for (const unsigned int i : ts_values.dof_indices())
        {
          for(const unsigned int j : ts_values.dof_indices())
            cell_matrix(i,j) +=
                      ((ts_values.shape_grad(i, q)         // ((\nabla \phi_i
                       * ts_values.shape_grad(j, q))       //   * \nabla \phi_j)                                        
                       +                                   //  +
                      (ts_values.shape_value(i, q)         //  (\phi_i
                       * k * std::pow(old_value[q], k-1)   //  * k * u_n^k-1
                       * ts_values.shape_value(j, q)) )    //     * \phi_j)
                     * ts_values.JxW(q);                   // * dx


          cell_rhs(i) +=  (ts_values.shape_value(i, q)              //   (\phi_i
                              * (rhs - std::pow(old_value[q], k))   //      * f - u_n^k
                            - ts_values.shape_grad(i, q)            //    - \nabla \phi_i   
                              * old_gradient[q])                    //      * \nabla u_n)
                          * ts_values.JxW(q);                       // * dx
        } // for ( i )      
      } // for ( q )

      // Check for neumann conditions
      if (cell ->at_boundary())
      {
        for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++)
        {
          if (cell -> face(f) -> at_boundary()
              && (cell -> face(f) -> boundary_id() == Boundary::Neumann_Case_1
              || cell -> face(f) -> boundary_id() == Boundary::Neumann_Case_2
              || cell -> face(f) -> boundary_id() == Boundary::Neumann_Case_3))
          {
            face_values.reinit(cell, f);

            for (const unsigned int q : face_values.quadrature_point_indices())
            {
              double g_value;
              if ( problem_case == ProblemCase::Case_1)
                g_value = nc1_fcn.value(face_values.quadrature_point(q));
              else if ( problem_case == ProblemCase::Case_2)
                g_value = nc2_fcn.value(face_values.quadrature_point(q));
              else if ( problem_case == ProblemCase::Case_3)
                g_value = nc3_fcn.value(face_values.quadrature_point(q));
              else
                AssertThrow(false, ExcNotImplemented());

              for (const unsigned int i : face_values.dof_indices())
                cell_rhs(i) +=  (face_values.shape_value(i, q)    //  \phi_i
                                  *  g_value)                     //  * g
                                * face_values.JxW(q);             // * dx
              
            } // for ( q )
          } // if ( Neumann)

        } // for ( f )
      } // if

      // Add the values to the system
      std::vector< unsigned int > local_dof_indices =
              tria.get_IEN_array(cell);

      system_matrix.add(local_dof_indices, cell_matrix);
      system_rhs.add(local_dof_indices, cell_rhs);

    } // for ( cell )


    // Impose zero Dirichlet boundary condition to the system for
    // all Dirichlet boundary dofs. 
    unsigned int n_global_dofs = tria.n_active_splines();
    const auto& boundary_dofs = tria.get_boundary_dofs();

    if (boundary_dofs.find(Boundary::Dirichlet_0) 
          != boundary_dofs.end())
    {    
      for (const auto& dof : boundary_dofs.at(Boundary::Dirichlet_0)){
        for (unsigned int i = 0; i < n_global_dofs; i++){
          system_matrix.set(i, dof, 0.);
          system_matrix.set(dof, i, 0.);
        }
        system_matrix.set(dof, dof, 1.);
        system_rhs(dof) =  0.;
      }
    }
     
  } // assemble system


  void Nonlinear2D_Benchmark::get_neumann_data(
      std::map< types::boundary_id,
                const Function<2>* >& neumann_data
  ) {
      if (this -> problem_case == ProblemCase::Case_1)
        neumann_data = {{Boundary::Neumann_Case_1, &nc1_fcn}};
      else if (problem_case == ProblemCase::Case_2)
        neumann_data = {{Boundary::Neumann_Case_2, &nc2_fcn}};
      else if (problem_case == ProblemCase::Case_3)
        neumann_data = {{Boundary::Neumann_Case_3, &nc3_fcn}};
      else
        AssertThrow(false, ExcNotImplemented());
  } // get_neumann_data


  void Nonlinear2D_Benchmark::impose_boundary_condition()
  {
    std::cout << "    Imposing boundary condition to system ..." << std::endl;
    const auto& boundary_dofs = tria.get_boundary_dofs();

    // For the case the initial mesh has no boundary dofs
    if (boundary_dofs.size() == 0)
      return;

    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = degrees[0]  + 1;
    degrees[1] = degrees[1]  + 1;

    std::map<
        types::global_dof_index, 
        double
      > boundary_values;

    // Secondly set the zero Dirichlet boundary values
    if (boundary_dofs.find(Boundary::Dirichlet_0) 
          != boundary_dofs.end())
      for (const auto& dof : boundary_dofs.at(Boundary::Dirichlet_0))
        boundary_values[dof] = 0.;
    
    MatrixTools::apply_boundary_values(
      boundary_values, 
      system_matrix,
      current_solution,
      system_rhs
    );
  } // impose_boundary_condition



  double Nonlinear2D_Benchmark::determine_step_length_LS() const
  { 
    // Room for improvement. LS via deal.II
    return 0.1;
  } // determine_step_length



  void Nonlinear2D_Benchmark::estimate_and_mark()
  {
    std::cout << "    Refining grid..." << std::endl;
    if (refinement_strategy == RefinementStrategy::Uniform) 
    {
      tria.coarsen_bezier_elements();
      tria.refine_global();   
    } 
    else // RefinementStrategy::Adaptive
    {
      std::vector< TriaIterator<CellAccessor<2, 2>> > mark;
      const std::vector<unsigned int>& degrees = tria.get_degree();
      Vector<double>  local_residuals(tria.n_active_cells());

      std::map< types::boundary_id,
              const Function<2>* >   neumann_data = {};
      get_neumann_data(neumann_data);
      // for (const auto& [boundary, fcn] : neumann_data)
      //   std::cout << " boundary: " << boundary << std::endl;
      // std::cout << "neumann_data size: "<< neumann_data.size() << std::endl;
      
      tria.nonlinear_residual_error_estimate(
                    {degrees[0]*degrees[0] + 1,
                     degrees[1]*degrees[1] + 1},
                     exponent,
                     rhs_fcn,
                     neumann_data,
                     current_solution,
                     local_residuals
                     );

      // tria.poisson_residual_error_estimate(
      //           {degrees[0]*degrees[0] + 1,
      //            degrees[1]*degrees[1] + 1},
      //            rhs_fcn,
      //            neumann_data,
      //            current_solution,
      //            local_residuals
      //            );
      
      tria.refine_fixed_number(local_residuals, 0.20);
    }

    tria.prepare_assembly();
  } // estimate_and_mark


  
  void Nonlinear2D_Benchmark::output_system()
  {
    std::cout << "Printing system matrix and rhs ... " << std::endl;
  
    const unsigned int level = tria.n_levels() - 1;
    const std::string name = problem_out.degree.string();
    const std::string level_name = name + "l" + std::to_string(level);


    
    std::string matrix = level_name + "_mat.dat" ;
    std::string vector = level_name + "_vec.dat" ;
    std::string soluti = level_name + "_sol.dat" ;

    // // Find some interesting spline index
    // if (tria.n_levels() - 1 == 9) {
    //   const double target_x = 0.625;
    //   const double target_y_min = 0.875;
    //   const double target_y_max = 1.;

    //   const auto& splines = tria.get_splines(); 

    //   int index = -1;
    //   for (const auto& t : splines) {
    //     const auto& A0 = t -> get_anchor(0);
    //     const auto& A1 = t -> get_anchor(1);
        
    //     if (A0.first == target_x && 
    //           A1.first == target_y_min && 
    //           A1.second == target_y_max)
    //       index = t -> get_level_index();
    //   }

    //   std::cout << "index = " << index << std::endl;
    // }
  
    if (tria.n_levels() - 1 < 17) {
      std::filebuf mat, vec, sol;
      mat.open(matrix.c_str(), std::ios::out);
      vec.open(vector.c_str(), std::ios::out);
      sol.open(soluti.c_str(), std::ios::out);
  
      std::ostream mat_out(&mat);
      std::ostream vec_out(&vec);
      std::ostream sol_out(&sol);
  
      system_matrix.print_formatted(mat_out, 16, true, 1, "0");
      system_rhs.print(vec_out, 16);
      current_solution.print(sol_out, 16);
  
      mat.close();
      vec.close();
      sol.close();
  
      const auto& kv = data.kv;
      const int nx = kv[0].size();
      const int ny = kv[1].size();
  
      const unsigned int N1 = 100; 
      const unsigned int N2 = 100;
      const double xmin = kv[0][0];
      const double ymin = kv[1][0]; 
      const double xmax = kv[0][nx-1];
      const double ymax = kv[1][ny-1];
      std::vector<Point<2>> evals; 
      FullMatrix< double > E(N1 * N2, 2);
      unsigned int ind = 0;
      for (double j = 0.; j < N2; j++) {
        for (double i = 0.; i < N1; i++) {
          const double x = xmin + (i/(N1-1.)) * (xmax - xmin); 
          const double y = ymin + (j/(N2-1.)) * (ymax - ymin); 
          evals.push_back(Point<2>(x, y));
          E(ind, 0) = x;
          E(ind, 1) = y;
          ind++;
        }
      }

      std::filebuf e_f;
      e_f.open(level_name + "_evals.dat", std::ios::out);
      std::ostream out_e(&e_f);
      E.print_formatted(out_e, 16, true, 1, "0");
      e_f.close();

      
      // Print the IPF wireframe:
      // tria.print_grid(name);
      tria.generate_mesh_file<0>(level_name, false, 16);
      tria.generate_mesh_file<0>(level_name, true, 16);
      tria.print_IPF_wireframe(level_name);
      tria.printIPF(evals, level_name, 16, true, true);
      tria.coarsen_bezier_elements();
      tria.generate_mesh_file<0>(level_name, false, 16);
      tria.generate_mesh_file<0>(level_name, true, 16);
      tria.print_IPF_wireframe(level_name);
      tria.refine_bezier_elements();
    }
    

    // Write the grid to a seperate file: 
    const std::string& name_vtg = problem_out.vtg.string() + "physical_grid_l" + std::to_string(level) + ".vtu";

    // First: Make a copy of the triangulation: 
    Triangulation<2> physical_grid; 
    physical_grid.copy_triangulation(tria);

    // And transform it with the IPF from tria.
    // Note: This will make a linear representation
    // of the boundary and interior nodes.
    const IsoparametricManifold<2> geometry(tria.get_IPF()); 
    GridTools::transform(
      [&geometry](const Point<2>& p){
        return geometry.push_forward(p);
      },
      physical_grid 
    );

    // Generate the output object
    DataOut<2> data_out;
    data_out.attach_triangulation(physical_grid); 

    
    const std::vector<unsigned int>& degrees = tria.get_degree();
    Vector<double>  cell_errors(tria.n_active_cells());

    std::map< types::boundary_id,
              const Function<2>* >   neumann_data = {};
      get_neumann_data(neumann_data);
      
      tria.nonlinear_residual_error_estimate(
                          {degrees[0]*degrees[0] + 1, 
                           degrees[1]*degrees[1] + 1},
                           exponent,
                           rhs_fcn,
                           neumann_data,
                           current_solution,
                           cell_errors
                           );

    data_out.add_data_vector(cell_errors, "cell_errors");
    

    // Build patches
    data_out.build_patches(); 

    // Open a file
    std::ofstream vtu_out(name_vtg); 
    data_out.write_vtu(vtu_out);



    // print the difference pointwise for each cell
    std::cout << " ... output to files done!" << std::endl;
    





    
    GridOutFlags::Svg svg_flags;
    svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
    // svg_flags.label_level_number  = true;
    // svg_flags.label_cell_index    = true;
    svg_flags.label_boundary_id   = true;
    
    std::string svg_name = problem_out.svg.string() + "parametric_grid_l"
                          + std::to_string(tria.n_levels() - 1)
                          + ".svg";  
    //std::cout << "... printing to: " << svg_name << std::endl;
    std::ofstream out(svg_name);
    GridOut       grid_out;
    grid_out.set_flags(svg_flags);
    grid_out.write_svg(tria, out);
    out.close();
  } // output_system



  void Nonlinear2D_Benchmark::solve_system()
  {
    
    newton_update = system_rhs;
    SparseDirectUMFPACK A_direct;
    A_direct.solve(system_matrix, newton_update);


    // Add the update to the current solution
    double alpha = determine_step_length_const();
    current_solution.add(alpha, newton_update);

    // Add current solution values to the Tsplines data structure
    const auto& splines = tria.get_splines();
    for (unsigned int i = 0; i < tria.n_active_splines(); i++)
        splines[i] -> set_solution(current_solution[i]);

  } // solve_system


  void Nonlinear2D_Benchmark::print_numerical_solution(
    Vector<double> vector,
    std::string addition
  ){
    
    // Number of points to evaluate per direction
    const unsigned int N = 31;
    const unsigned int actual_N = std::pow(N, 2);

    // Declare a container to store the values 
    FullMatrix<double> B_num(3,actual_N);

    const auto& splines = tria.get_splines();
    const auto& mapping = tria.get_IPF();
    const auto& kv      = data.kv;
    const unsigned int n0 = kv[0].size() - 1;
    const unsigned int n1 = kv[1].size() - 1;

    {
      for (unsigned int i = 0; i < N; i++)
      {
        for (unsigned int j = 0; j < N; j++)
        {
          const Point<2> P(kv[0][0] + (kv[0][n0] -  kv[0][0]) * i / (N-1.), 
                            kv[1][0] + (kv[1][n1] -  kv[1][0]) * j / (N-1.));
          B_num(0, i*N + j) = P[0];
          B_num(1, i*N + j) = P[1];
          for (unsigned int k = 0; k < tria.n_active_splines(); k++)
          {
            B_num(2, i*N + j) += vector[k] * splines[k] -> value(P);
            //B_num(2, i*N + j) += vector[k];
          }
        }
      }
    }
    
    
    std::string name_num =  problem_out.dat.string() 
                            + "numsol_" + addition + "l" 
                            + std::to_string(tria.n_levels() - 1) 
                            + ".dat";
    //std::cout << name_num << std::endl;
    std::filebuf f0;
    f0.open(name_num.c_str(), std::ios::out);

    std::ostream out0(&f0);
    B_num.print_formatted(out0, 16, true, 1, "0");

    f0.close();
    
  }







  void Nonlinear2D_Benchmark::run()
  {
    std::cout << "Running benchmark ... " << std::endl;
    unsigned int old_level = 0;
    const std::vector<unsigned int>& degrees = tria.get_degree();
    setup_system();

    

    // Set initial solution u_0 to zero. 
    current_solution = 0.;

    // get Neumann data for estimating error
    Vector<double> current_residuals;
    std::map< types::boundary_id,
              const Function<2>* >   neumann_data = {};
    get_neumann_data(neumann_data);    

    // boundary condition on the first solution u_0 and system
    impose_boundary_condition();


    // Refinement cycle
    unsigned int cycle = 0;
    while (tria.n_levels() < this -> ref + 1)
    {
      std::cout << "  Refinement cycle " << cycle << ':' << std::endl;
      if (cycle != 0)
      {
        estimate_and_mark();
        setup_system();
        tria.transfer_solution(current_solution);
        impose_boundary_condition();
      }

      // Estimate residual before Newton
      current_residuals.reinit(tria.n_active_cells());
      tria.nonlinear_residual_error_estimate(
              {degrees[0]*degrees[0] + 1,
               degrees[1]*degrees[1] + 1},
               exponent,
               rhs_fcn,
               neumann_data,
               current_solution,
               current_residuals
               );

      std::cout << " On Refinement level: " << tria.n_levels() << std::endl;
      
      double initial_estimator = current_residuals.l2_norm();
      std::cout << "    Initial estimated norm: " 
                << initial_estimator << std::endl;

      std::string name_norm =  problem_out.dat.string() 
                            + "norm_l" 
                            + std::to_string(tria.n_levels() - 1) 
                            + ".dat";
      // outputs iteration and norms
      std::ofstream norm_file(name_norm.c_str(), std::ios::app); // Open file in append mode



      unsigned int newton_iteration = 0;
      double  current_norm = 1., current_residual = 1.;
      do {
        // solve system with zero boundary condition
        assemble_system();
        solve_system();
        current_residuals = system_rhs;


        current_residual = system_rhs.l2_norm();
        std::cout << "     ||sysem_rhs|| " << newton_iteration+1 
                  << ":   " << std::fixed << std::setprecision(8) 
                  << current_residual;

        current_norm = newton_update.l2_norm();
        std::cout << "      ||update_n+1|| " 
                  << ":   " << std::fixed << std::setprecision(8) 
                  << current_norm << std::endl;

        if (norm_file.is_open()) {
          norm_file << newton_iteration + 1 << " " // Add iteration number for reference
              << std::fixed << std::setprecision(8) << current_residual << " "
              << std::fixed << std::setprecision(8) << current_norm << "\n";
        } else 
            std::cerr << "Error: Unable to open file for writing.\n";
        
        newton_iteration++;
      } // do ( ... )
      while (newton_iteration < 300
             && current_norm > 1e-8
            );

      norm_file.close(); // Close the file

      current_residuals.reinit(tria.n_active_cells());
      tria.nonlinear_residual_error_estimate(
                    {degrees[0]*degrees[0] + 1,
                     degrees[1]*degrees[1] + 1},
                     exponent,
                     rhs_fcn,
                     neumann_data,
                     current_solution,
                     current_residuals
                     );

      double last_estimator = current_residuals.l2_norm();
      std::cout << "   Estimated Norm: " 
                << last_estimator << "\n\n" << std::endl;

      if (tria.n_levels() != old_level)
      {
        std::cout << "Outputting system after last newton iteration ..." << std::endl;

        print_numerical_solution(current_solution);
        std::cout << " Print table in cycle: " << cycle << std::endl; 
        problem_out.add_values_to_table(
          tria.n_levels() - 1,
          cycle,
          tria.n_active_cells(),
          tria.n_active_splines(),
          newton_iteration,
          current_norm,
          last_estimator,
          initial_estimator
        );
        // Output the table to a file preemptively
        problem_out.write_table_text();
        problem_out.write_table_tex();
        if (tria.n_levels() < 17)
          output_system();
      }
      old_level = tria.n_levels();
      cycle++;
    }
    std::cout << " ... done!" << std::endl;
  } // run


//===================================================================

//===================================================================
  double Nonlinear2D_RHS1::value(
    const Point<2>&     /*p*/, 
    const unsigned int /* component */
  ) const {
    double out = 0.;
    return out;
  } // value 

  double Nonlinear2D_RHS2::value(
    const Point<2>&     /*p*/, 
    const unsigned int /* component */
  ) const {
    double out = 1.;
    return out;
  } // value 

  double Nonlinear2D_RHS3::value(
    const Point<2>&     /*p*/, 
    const unsigned int /* component */
  ) const {
    double out = 0.;
    return out;
  } // value 


  double Nonlinear2D_NC1::value(
    const Point<2>&     /*p*/, 
    const unsigned int /* component */
  ) const {
    double out = 1.;

    return out;
  } // value 

  double Nonlinear2D_NC2::value(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    double out = 0.;
    out = std::sin(4. * numbers::PI * p[1]);
    return out;
  } // value 

  double Nonlinear2D_NC3::value(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    double out = 0.;

    return out;
  } // value 






  // =================================================================

  // =================================================================


  
  IPF_Data<2> Nonlinear2D_Benchmark::get_IPF_data() const
  {
    std::vector< std::vector< double > > kv;
    std::vector< Point<2 + 1> > cps;
    std::vector< unsigned int > deg;


     // define the knot vectors
    kv = std::vector< std::vector< double > >(2);
    kv[0] = {0, 0, 1, 1};
    kv[1] = {0, 0, 1, 1};


    // define the control points vector
    cps = std::vector< Point<3> >(4);
    cps[0] = Point<2 + 1>( 0.,  0.,  1.);
    cps[1] = Point<2 + 1>( 1.,  0.,  1.);
    cps[2] = Point<2 + 1>( 0.,  1.,  1.);
    cps[3] = Point<2 + 1>( 1.,  1.,  1.);

    // define the degrees
    deg = {1, 1};

    return IPF_Data<2>(cps, kv, deg);
  } // get_IPF_data 2-dim

} // namespace Nonlinear