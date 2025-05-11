#include <minimal_test.h>



namespace Minimal_Surface{


  Minimal_Test::Minimal_Test(
    int ref,
    int order,
    const ProblemShape shape,
    const RefinementStrategy strategy

    ) : ref(ref)
      , order(order)
      , data(this->get_IPF_data())
      , tria(data)
      , cycle(0)
      , problem_shape(shape)
      , refinement_strategy(strategy)
    {    
      // Assembling the file name and elevating the degree
      std::string name = "minimal_surface";
      if (shape == ProblemShape::Halfcircle)
      {
        tria.degree_elevate(1, 1);
        name += "_halfcircle";
      }
      else if (shape == ProblemShape::Tilted_Halfcircle)  
      {
        tria.degree_elevate(1, 1);
        name += "_tilted_halfcircle";
      }
      else if (shape == ProblemShape::Annulus)
      {
        tria.degree_elevate(1, 1);
        
        name += "_annulus";

      }
      else if (shape == ProblemShape::Square)
        name += "_square";
      else
        AssertThrow(false, ExcNotImplemented());



      if (strategy == RefinementStrategy::Adaptive)
        name += "_adaptive";
      else
        name += "_uniform";
      
      tria.degree_elevate_global(order);
      tria.refine_global();
      tria.refine_bezier_elements();
      const auto& cell = tria.begin_active();
      const auto& CPs = tria.get_control_points(cell);
      tria.coarsen_bezier_elements();

      tria.prepare_assembly();

      
      // Set the output file name for annulus
      const std::vector<std::string>& columns2 =
        {"Level", "Cycles", "Cells", "DoFs", "k_{Newton}", "update_norm", "initial_residual", "last_residual", "L2", "H1"};
      const std::vector<std::string>& tex_captions2 = 
        {"Level", "Cycles", "# cells", "# dofs", "k_{Newton}", "update_norm", "initial_residual", "last_residual", "$L_2$-error", "$H^1$-error"};
      const std::vector<bool>& scientific2 =
        {false, false, true, true, true, true, true, true, true, true};
      const std::vector<unsigned int> precision2 = 
        {0, 0, 0, 0, 0, 2, 2, 2, 2, 2};
      const std::vector<std::string>& super_column_names2 = 
        {"Grid Info", "Newton", "Errors"};
      const std::vector<std::vector<std::string>> super_columns2 =
        {{"Level", "Cycles", "Cells", "DoFs"}, {"k_{Newton}", "update_norm", "initial_residual", "last_residual"}, {"L2", "H1"}};
      
      const std::vector<std::string>& columns =
        {"Level", "Cycles", "Cells", "DoFs", "k_{Newton}", "update_norm", "initial_residual", "last_residual"};
      const std::vector<std::string>& tex_captions = 
        {"Level", "Cycles", "\\# cells", "\\# dofs", "k_{Newton}", "update_norm", "initial_residual", "last_residual"};
      const std::vector<bool>& scientific =
        {false, false, true, true, true, true, true, true};
      const std::vector<unsigned int> precision = 
        {0, 0, 0, 0, 0, 2, 2, 2};
      const std::vector<std::string>& super_column_names = 
        {"Grid Info", "Newton"};
      const std::vector<std::vector<std::string>> super_columns =
        {{"Level", "Cycles", "Cells", "DoFs"}, {"k_{Newton}", "update_norm", "initial_residual", "last_residual"}};

      if (shape == ProblemShape::Annulus)
      {
        problem_out = OutputSetup(  name
                                  , data.max_degree() + order
                                  , columns2
                                  , tex_captions2
                                  , scientific2
                                  , precision2
                                  , super_column_names2
                                  , super_columns2
                                  ); 
        //problem_out.table.set_auto_fill_mode(true);
      }
      else
        problem_out = OutputSetup(  name
                                  , data.max_degree() + order
                                  , columns
                                  , tex_captions
                                  , scientific
                                  , precision
                                  , super_column_names
                                  , super_columns
                                  ); 



      // Set boundary indicators
      for (auto& face : tria.active_face_iterators())
      { 
        std::cout << "face: " << face -> index() << " with coords ";
        std::cout << face->center()[0] << ", "<< face->center()[1]<<std::endl;
      if (!face -> at_boundary())
        continue;
      else if (problem_shape == ProblemShape::Square)
      {
        const Point<2>& c = face -> center();
        if (std::fabs(c(0)) < 1e-15)
          face -> set_boundary_id(Boundary::Dirichlet_0);
        else
          face -> set_boundary_id(Boundary::Dirichlet_0);
      }
        else
          face -> set_boundary_id(Boundary::None);
      } // for ( face )
    
    } // Minimal_Test



  void Minimal_Test::setup_system()
  {

    std::cout << "Setting up system ... " << std::endl;
    unsigned int n_global_dofs = tria.n_active_splines();
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
    system_rhs.reinit(n_global_dofs);


  } // setup_system


  void Minimal_Test::assemble_system()
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


    const int k = 1;

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
        // const double rhs = rhs_fcn -> value(mapped_q);
        const double rhs = mapped_q[0];
        double u_nk= 1;
        double u_nk_1= k;
        if(k!=0)
        {
          u_nk=std::pow(old_value[q], k);
        }
        if(k!=1)
        {
          u_nk_1=k*std::pow(old_value[q], k-1);
        }

        // Build the cell matrix and rhs
        for (const unsigned int i : ts_values.dof_indices())
        {
          for(const unsigned int j : ts_values.dof_indices())
          {
            cell_matrix(i,j) +=
                      (ts_values.shape_grad(i, q)          // ((\nabla \phi_i
                       * ts_values.shape_grad(j, q))       //   * \nabla \phi_j)                                        
                       +                                   //  +
                      (ts_values.shape_value(i, q)         //  (\phi_i
                       * u_nk_1                            //    k * u_n^k-1
                       * ts_values.shape_value(j, q))      //     * \phi_j)
                     * ts_values.JxW(q);                   // * dx
          } // for ( j )

          
          cell_rhs(i) +=  (ts_values.shape_value(i, q)     //   (\phi_i
                            * (rhs - u_nk)                 //      * f - u_n^k
                          - ts_values.shape_grad(i, q)     //    - \nabla \phi_i   
                            * old_gradient[q]              //      * \nabla u_n)
                          )
                          * ts_values.JxW(q);              // * dx
        } // for ( i )      
      } // for ( q )

      // Check for neumann conditions
      if (cell ->at_boundary())
      {
        for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; f++)
        {
          if (cell -> face(f) -> at_boundary()
              && cell -> face(f) -> boundary_id() == Boundary::Neumann_Case_1)
          {
            face_values.reinit(cell, f);

            for (const unsigned int q : face_values.quadrature_point_indices())
            {
              double g_value;
              if ( problem_shape == ProblemShape::Square)
                g_value = nc_fcn.value(face_values.quadrature_point(q));

              for (const unsigned int i : face_values.dof_indices())
                cell_rhs(i) +=  0*(face_values.shape_value(i, q)    //  \phi_i
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
        //std::cout << "setting boundary at dof: " << dof << std::endl;
        }
    }

  } // assemble_system



  
  void Minimal_Test::impose_boundary_condition()
  {
    
    std::cout << "Imposing boundary condition ..." << std::endl;
    const auto& boundary_dofs = tria.get_boundary_dofs();
    unsigned int n_global_dofs = tria.n_active_splines();
    // For the case the initial mesh has no boundary dofs
    if (boundary_dofs.size() == 0)
      return;
    // print number of dofs at boundary
    //std::cout << "boundary_dofs: " << std::endl;
    //for (const auto& [boundary, dofs] : boundary_dofs)
    //  std::cout << "Boundary " << boundary << ": " << dofs.size() << std::endl;
    //std::cout << "number of global dofs: " << tria.n_active_splines() << std::endl;

    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = degrees[0]  + 1;
    degrees[1] = degrees[1]  + 1;

    std::map<
        types::global_dof_index, 
        double
      > boundary_values;


    // Firstly set the nonzero Dirichlet boundary values
    if (this -> problem_shape == ProblemShape::Square)
    {
      if (boundary_dofs.find(Boundary::Dirichlet_s) 
          == boundary_dofs.end())
        return;
      std::map<
        types::boundary_id,
        const Function< 2 >*
      > boundary_fcns =
          {{Boundary::Dirichlet_s, &DC_fcn}};
      tria.project_boundary_values(
            boundary_fcns,
            degrees,
            boundary_values);

    }
    else
      AssertThrow(false, ExcNotImplemented());


    // Secondly set the zero Dirichlet boundary values
    if (boundary_dofs.find(Boundary::Dirichlet_0) 
          != boundary_dofs.end())
      for (const auto& dof : boundary_dofs.at(Boundary::Dirichlet_0))
        boundary_values[dof] = 0.;
    

    //TODO: Für circle sind die dofs noch Null obwohl Dirichlet_c gesetzt ist
    // Die Punkte, die in DC_c inputted werden, sind irgendwie alle Null???
    // Obwohl .vtg und matlab das richtige zeigt
    // print boundary values
    //Boundary bndry = Boundary::Dirichlet_a1;
    //if (boundary_dofs.find(Boundary::Dirichlet_a1) 
    //      != boundary_dofs.end())
    //{
    //  std::cout << "Boundary dofs of Dirichlet_a1: " << std::endl;
    //  const auto& splines = tria.get_splines(); 
    //  for (const auto& dof : boundary_dofs.at(bndry)){
    //    const auto& ts = splines.at(dof); 
    //    const auto& anchor = ts -> get_anchor();
    //    std::cout << dof << ": " << 0.5 * anchor.first + 0.5*anchor.second << ", value = " << boundary_values.at(dof) << std::endl;
    //  } // for ( dof )
    //}
    
    MatrixTools::apply_boundary_values(
      boundary_values, 
      system_matrix,
      current_solution,
      system_rhs
    );

    
  } // impose_boundary_condition



  void Minimal_Test::solve_system()
  {
    // std::cout << "Solving system ... " << std::endl;

    //SolverControl            solver_control(1e3,
    //                             1e-5);
    //SolverCG<Vector<double>> solver(solver_control);
    //
    //PreconditionJacobi<SparseMatrix<double>> preconditioner;
    //preconditioner.initialize(system_matrix);
    //
    //solver.solve(system_matrix, newton_update, system_rhs, preconditioner);

    // Solve with direct UMFPACK
    newton_update = system_rhs;
    SparseDirectUMFPACK A_direct;
    A_direct.solve(system_matrix, newton_update);


    // Add the update to the current solution
    double alpha = determine_step_length_const();
    //double alpha = determine_step_length_LS();
    //std::cout << "alpha = " << alpha << std::endl;
    current_solution.add(alpha, newton_update);

    // Add current solution values to the Tsplines data structure
    const auto& splines = tria.get_splines();
    for (unsigned int i = 0; i < tria.n_active_splines(); i++)
        splines[i] -> set_solution(current_solution[i]);
  } // solve_system









 

  void Minimal_Test::estimate_and_mark()
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

      tria.minimal_residual_error_estimate(
                          {degrees[0]*degrees[0] + 1,
                           degrees[1]*degrees[1] + 1},
                           current_solution,
                           local_residuals
                           );


      tria.refine_fixed_number(local_residuals, 0.10);
    }

    tria.prepare_assembly();
  } // estimate_and_mark








  void Minimal_Test::print_numerical_solution(std::string addition)
  {
    // Number of points to evaluate per direction
    const unsigned int N = 31;

    // Declare a container to store the values 
    FullMatrix<double> B_num(3,N*N);

    const auto& splines = tria.get_splines();
    const auto& mapping = tria.get_IPF();
    const auto& kv      = data.kv;
    const unsigned int n0 = kv[0].size() - 1;
    const unsigned int n1 = kv[1].size() - 1;   

    if (problem_shape == ProblemShape::Square)
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
            B_num(2, i*N + j) += current_solution[k] * splines[k] -> value(P);
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
  } // print_numerical_solution





  void Minimal_Test::run()
  {
    std::cout << "Running benchmark ... " << std::endl;
    unsigned int old_level = 0;

    setup_system();
    // Set initial solution u_0 to zero. 

    Vector<double> current_residuals;
    current_residuals = 0.;
    current_solution = 0.;
    const std::vector<unsigned int>& degrees = tria.get_degree();

    // boundary condition on the first solution u_0
    impose_boundary_condition();
    

    // Refinement cycle
    unsigned int cycle = 0;
    while (tria.n_levels() < this -> ref + 1)
    {
      std::cout << "  Refinement cycle " << cycle << ':' << std::endl;
      
      current_residuals.reinit(tria.n_active_cells());
      tria.minimal_residual_error_estimate(
                    {degrees[0]*degrees[0] + 1,
                     degrees[1]*degrees[1] + 1},
                     current_solution,
                     current_residuals
                     );


      std::cout << " n_levels: " << tria.n_levels() << std::endl;
      



      double initial_residual = current_residuals.l2_norm();
      std::cout << "    Initial norm of estimator: " << initial_residual << std::endl;


      std::string name_norm =  problem_out.dat.string() 
                            + "norm_cycle" 
                            + std::to_string(cycle) 
                            + ".dat";
      // outputs iteration and norms
      std::ofstream norm_file(name_norm.c_str(), std::ios::app); // Open file in append mode

      if (cycle != 0)
      {
        estimate_and_mark();
        setup_system();
        tria.transfer_solution(current_solution);
        impose_boundary_condition();

      }
      unsigned int newton_iteration = 0;
      double  current_norm = 1., current_residual = 1.;
      do {
        // solve system with zero boundary condition
        assemble_system();
        solve_system();
        current_residuals = system_rhs;
        


        current_residual = current_residuals.l2_norm();
        std::cout << "Residual of N_It " << newton_iteration+1 
                  << ":   " << std::fixed << std::setprecision(8) << current_residual;

        current_norm = newton_update.l2_norm();
        std::cout << "      L2 ||update_n+1|| " 
                  << ":   " << std::fixed << std::setprecision(8) << current_norm << std::endl;

        if (norm_file.is_open()) {
          norm_file << newton_iteration + 1 << " " // Add iteration number for reference
              << std::fixed << std::setprecision(8) << current_residual << " "
              << std::fixed << std::setprecision(8) << current_norm << "\n";
        } else 
            std::cerr << "Error: Unable to open file for writing.\n";
        
        
        
        newton_iteration++;

      } // do ( ... )
      while (newton_iteration < 200
             && current_norm > 1e-8
            );

      norm_file.close(); // Close the file
      current_residuals.reinit(tria.n_active_cells());
      tria.minimal_residual_error_estimate(
                    {degrees[0]*degrees[0] + 1,
                     degrees[1]*degrees[1] + 1},
                     current_solution,
                     current_residuals
                     );
      std::cout << "    Last norm of estimator: " 
                << current_residuals.l2_norm() << std::endl;
      std::cout << "    Last residual: " 
                << std::fixed << std::setprecision(5)
                << current_residual << "\n\n"<< std::endl;

      std::cout << "Outputting system after last newton iteration ..." << std::endl;
      if (cycle != 0 && tria.n_levels() != old_level)
      {
        print_numerical_solution();

        


      }
      old_level = tria.n_levels();
      cycle++;
    }
    std::cout << " ... done!" << std::endl;
  } // run




  // =================================================================

  // =================================================================
  double Minimal_RHS1::value(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    double out = 0.;//std::sin(2 * numbers::PI * (p[0] + p[1]));
    

    return out;
  } // value_DC_circle

  double Minimal_DC::value(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    double out = 0.;//std::sin(2 * numbers::PI * (p[0] + p[1]));
    

    return out;
  } // value_DC_circle


  double Minimal_NC1::value(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    double out = 1.;//std::sin(2 * numbers::PI * (p[0] + p[1]));
    

    return out;
  } // value_DC_circle










  // =================================================================

  // =================================================================
  IPF_Data<2, 2> Minimal_Test::get_IPF_data(
  ) {
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
    
  } // get_IPF_data

  
}