#include <minimal_surface.h>


namespace Minimal_Surface{


  Minimal_Benchmark::Minimal_Benchmark(
    int ref,
    int order,
    const ProblemShape shape,
    const RefinementStrategy strategy,
    const StepLengthStrategy steplength
    ) : ref(ref)
      , order(order)
      , data(this->get_IPF_data(shape))
      , tria(data)
      , cycle(0)
      , problem_shape(shape)
      , refinement_strategy(strategy)
      , step_length_strategy(steplength)
    {
      tria.degree_elevate(1, 1);
      tria.degree_elevate_global(order);
      tria.refine_bezier_elements();
      const auto& cell = tria.begin_active();
      const auto& CPs = tria.get_control_points(cell);
      tria.coarsen_bezier_elements();

      tria.prepare_assembly();

      // Assembling the file name
      std::string name = "minimal_surface";
      if (shape == ProblemShape::Halfcircle)
        name += "_halfcircle";
      else if (shape == ProblemShape::Tilted_Halfcircle)
        name += "_tilted_halfcircle";
      else if (shape == ProblemShape::Annulus)
        name += "_annulus";
      else if (shape == ProblemShape::Square)
        name += "_square";
      else
        AssertThrow(false, ExcNotImplemented());
      if (strategy == RefinementStrategy::Adaptive)
        name += "_adaptive";
      else
        name += "_uniform";
      
      problem_out = OutputSetup(name, data.max_degree() + order); 




      // Set boundary indicators
      for (auto& face : tria.active_face_iterators())
      { 
        std::cout << "face: " << face -> index() << " with coords ";
        std::cout << face->center()[0] << ", "<< face->center()[1]<<std::endl;
        if (!face -> at_boundary())
          continue;
        else if (shape == ProblemShape::Halfcircle || shape == ProblemShape::Tilted_Halfcircle)
          {
            const Point<2>& c = face -> center();
            if (std::fabs(c(1)- 1) < 1e-15)
              face -> set_boundary_id(Boundary::Dirichlet_c);
          }

        else if (shape == ProblemShape::Annulus) {
          const Point<2>& c = face -> center();
          if (std::fabs(c(1) + 1) < 1e-15)
          {
            face -> set_boundary_id(Boundary::Dirichlet_0);
          }
          if (std::fabs(c(1) - 1) < 1e-15)
          {
            face -> set_boundary_id(Boundary::Dirichlet_a);
          }
        }

        else if (shape == ProblemShape::Square) {
          const Point<2>& c = face -> center();
          if (std::fabs(c(1)) < 1e-15 || std::fabs(c(1)-1) < 1e-15)
            face -> set_boundary_id(Boundary::Dirichlet_s);
          if (std::fabs(c(0)) < 1e-15 || std::fabs(c(0)-1) < 1e-15)
            face -> set_boundary_id(Boundary::Dirichlet_0);
        }
        else
          face -> set_boundary_id(Boundary::None);
      } // for ( face )
    
    } // Minimal_Benchmark



  void Minimal_Benchmark::setup_system()
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


  void Minimal_Benchmark::assemble_system()
  {
    //std::cout << "Assembling system matrix ... " << std::endl;
  
    // Setup initial tables that store the bernstein values / grads / hessians.
  
    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = degrees[0]  + 1;
    degrees[1] = degrees[1]  + 1;

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
        update_gradients |
        update_quadrature_points |
        update_JxW_values);

    const unsigned int dofs_per_cell = ts_values.n_dofs_per_cell();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);



    for (const auto& cell : tria.active_cell_iterators())
    {
      // Get the Bernstein values on the cell
      ts_values.reinit(cell);

      // Get the old_solution_gradient for the current cell
      std::vector< Tensor<1, 2> > old_gradient(ts_values.n_quadrature_points_per_cell());
      ts_values.get_function_gradients(current_solution, old_gradient);

      // Reset the cell matrix
      cell_matrix       = 0;
      cell_rhs          = 0;

      // Quadrature sum:
      for (const unsigned int q : ts_values.quadrature_point_indices())
      {
        // coefficient a_n depending on old solution gradient for readability
        const double coeff = 1. / std::sqrt(1 + old_gradient[q] * old_gradient[q]);
        //std::cout << "old_gradient: " <<old_gradient[q][0] << " " << old_gradient[q][1] << std::endl;
        //std::cout << "coeff: " << coeff << std::endl;

        // Build the cell matrix and rhs
        for (const unsigned int i : ts_values.dof_indices())
        {
          
          for(const unsigned int j : ts_values.dof_indices())
          {
            cell_matrix(i,j) +=
                    (((ts_values.shape_grad(i, q)      // ((\nabla \phi_i
                       * coeff                         //   * a_n
                       * ts_values.shape_grad(j, q))   //   * \nabla \phi_j)
                      -                                //  -
                      (ts_values.shape_grad(i, q)      //  (\nabla \phi_i
                       * coeff * coeff * coeff         //   * a_n^3
                       * (ts_values.shape_grad(j, q)   //   * (\nabla \phi_j
                            * old_gradient[q])         //       * \nabla u_n)
                         * old_gradient[q]))           //     * \nabla u_n))
                     * ts_values.JxW(q));              // * dx
            //std::cout << "cell_matrix: " << cell_matrix(i,j) << std::endl;
            
          } // for ( j )

          cell_rhs(i) -= ((ts_values.shape_grad(i, q)             // \nabla \phi_i
                              * coeff                             //    * a_n
                              * old_gradient[q])                  //    * \nabla u_n
                            * ts_values.JxW(q));                  // * dx
        } // for ( i )
      } // for ( q )

      // Add the values to the system
      std::vector< unsigned int > local_dof_indices =
              tria.get_IEN_array(cell);
      

      system_matrix.add(local_dof_indices, cell_matrix);
      system_rhs.add(local_dof_indices, cell_rhs);
      
    } // for ( cell )


    // Impose zero Dirichlet boundary condition to the system for
    // all boundary dofs. Newton_update is zero at full boundary
    unsigned int n_global_dofs = tria.n_active_splines();
    const auto& boundary_dofs = tria.get_boundary_dofs();


    for (const auto& [boundary, dofs] : boundary_dofs){
      for (const auto& dof : dofs){
        for (unsigned int i = 0; i < n_global_dofs; i++){
          system_matrix.set(i, dof, 0.);
          system_matrix.set(dof, i, 0.);
        }
        system_matrix.set(dof, dof, 1.);
        system_rhs(dof) =  0.;
        //std::cout << "setting boundary at dof: " << dof << std::endl;
      }
    }

    // print system_matrix and system_rhs


    //std::cout << "System rhs: " << std::endl;
    //for (unsigned int i = 0; i < n_global_dofs; i++)
    //  std::cout << system_rhs(i) << " ";
    //std::cout << std::endl;

        

  } // assemble_system

  void Minimal_Benchmark::impose_boundary_condition(ProblemShape shape)
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
    //  std::cout << boundary << ": " << dofs.size() << std::endl;
    //std::cout << "number of global dofs: " << tria.n_active_splines() << std::endl;

    std::vector< unsigned int > degrees = tria.get_degree();
    degrees[0] = degrees[0]  + 1;
    degrees[1] = degrees[1]  + 1;

    std::map<
        types::global_dof_index, 
        double
      > boundary_values;


    // Firstly set the nonzero Dirichlet boundary values
    if (shape == ProblemShape::Square)
    {
      if (boundary_dofs.find(Boundary::Dirichlet_s) 
          == boundary_dofs.end())
        return;
      std::map<
        types::boundary_id,
        const Function< 2 >*
      > boundary_fcns =
          {{Boundary::Dirichlet_s, &square_DC_fcn}};
      tria.project_boundary_values(
            boundary_fcns,
            degrees,
            boundary_values
      );
    }
    else if (shape == ProblemShape::Annulus)
    {
      if (boundary_dofs.find(Boundary::Dirichlet_a) 
          == boundary_dofs.end())
        return;
      std::map<
        types::boundary_id,
        const Function< 2 >*
      > boundary_fcns =
          {{Boundary::Dirichlet_a, &cat_DC_fcn}};
      tria.project_boundary_values(
            boundary_fcns,
            degrees,
            boundary_values
      );
    }
    else if (shape == ProblemShape::Tilted_Halfcircle 
          || shape == ProblemShape::Halfcircle)
    {
      if (boundary_dofs.find(Boundary::Dirichlet_c) 
          == boundary_dofs.end())
        return;
      std::cout << "Bin in Dirichlet_c in impose" << std::endl;
      std::map<
        types::boundary_id,
        const Function< 2 >*
      > boundary_fcns =
          {{Boundary::Dirichlet_c, &circ_DC_fcn}};
      tria.project_boundary_values(
            boundary_fcns,
            degrees,
            boundary_values
      );
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
    //std::cout << "Boundary dofs of Dirichlet_c: " << std::endl;
    //const auto& splines = tria.get_splines(); 
    //for (const auto& dof : boundary_dofs.at(Boundary::Dirichlet_c)){
    //  const auto& ts = splines.at(dof); 
    //  const auto& anchor = ts -> get_anchor();
    //  std::cout << dof << ": " << 0.5 * anchor.first + 0.5*anchor.second << ", value = " << boundary_values.at(dof) << std::endl;
    //} // for ( dof )
    
    MatrixTools::apply_boundary_values(
      boundary_values, 
      system_matrix,
      current_solution,
      system_rhs
    );
    
  } // impose_boundary_condition




  double Minimal_Benchmark::determine_step_length_LS() const
  { 
    // Room for improvement. LS via deal.II
    return 0.1;
  }



  void Minimal_Benchmark::solve_system()
  {
    // std::cout << "Solving system ... " << std::endl;

    SolverControl            solver_control(system_rhs.size(),
                                 system_rhs.l2_norm() * 1e-6);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionJacobi<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix);

    solver.solve(system_matrix, newton_update, system_rhs, preconditioner);

    // print newton_update
    //std::cout << "newton_update: " << std::endl;
    //for (unsigned int i = 0; i < tria.n_active_splines(); i++)
    //  std::cout << newton_update[i] << " ";
    //std::cout << std::endl;

    // Add the update to the current solution
    double alpha = determine_step_length_const();
    current_solution.add(alpha, newton_update);

    // Add current solution values to the Tsplines data structure
    const auto& splines = tria.get_splines();
    for (unsigned int i = 0; i < tria.n_active_splines(); i++)
        splines[i] -> set_solution(current_solution[i]);


  } // solve_system


  double Minimal_Benchmark::compute_residual_for_steplength(double alpha)
  {
    // Compute the discrete residual for the current solution
    double n_global_dofs = tria.n_active_splines();
    Vector<double> evaluation_points(n_global_dofs);
    evaluation_points = current_solution;
    evaluation_points.add(alpha, newton_update);

    Vector<double> cell_residuals(n_global_dofs);
    const std::vector<unsigned int>& degrees = tria.get_degree();

    TSValues<2> ts_values(
        &tria,
        degrees,
        update_values |
        update_gradients |
        update_quadrature_points |
        update_JxW_values);

    const unsigned int nvc = GeometryInfo<2>::vertices_per_cell;
    //const unsigned int nvf = GeometryInfo<dim>::vertices_per_face;
    
    
    // Loop over all cells
    for (const auto& cell : tria.active_cell_iterators()) 
    {

      std::vector< unsigned int > local_dof_indices = tria.get_IEN_array(cell);
      ts_values.reinit(cell);

      // get old solution gradient
      std::vector< Tensor<1, 2> > sol_grad(ts_values.n_quadrature_points_per_cell());
      ts_values.get_function_gradients(evaluation_points, sol_grad);

      // Quadrature sum:
      for (const unsigned int q : ts_values.quadrature_point_indices())
      {
        // coefficient a_n depending on old solution gradient for visibilities
        const double coeff = 1. / std::sqrt(1 + sol_grad[q] * sol_grad[q]);

        // Fehler aber keine Ahnung warum bei ts_values.dof_indices()
        for (const unsigned int i : ts_values.dof_indices())
        {

          // root and squares cancel out -> res_cell = cell_width^2 * norm^2
          cell_residuals(i) += (ts_values.shape_grad(i, q)        // \nabla \phi_i
                               * coeff                        // * a_n
                               * sol_grad[q])                 // * \nabla u_n
                               * ts_values.JxW(q);            // * dx
        
        } // for ( i )

      } // for ( q )
      
      //const double cell_width = (cell->vertex(0)).distance(
      //                          cell->vertex(nvc - 1));
      //cell_residual *= cell_width;

    } // for (cell)

    // Compute the L2 norm of the residual
    double residual_norm = 0.;
    for (unsigned int i = 0; i < n_global_dofs; ++i)
      residual_norm += cell_residuals(i) * cell_residuals(i);

    return std::sqrt(residual_norm);
  } // compute_residual_for_steplength


  void Minimal_Benchmark::estimate_and_mark()
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



  void Minimal_Benchmark::output_system()
  {
    std::cout << "Printing system matrix and rhs ... " << std::endl;
  
    const unsigned int level = tria.n_levels() - 1;
    const std::string name = problem_out.degree.string();
    const std::string level_name = name + "l" + std::to_string(level);


    ///*
    std::string matrix = level_name + "_mat.dat" ;
    std::string vector = level_name + "_vec.dat" ;
    std::string soluti = level_name + "_sol.dat" ;

    // Find some interesting spline index
    if (tria.n_levels() - 1 == 9) {
      const double target_x = 0.625;
      const double target_y_min = 0.875;
      const double target_y_max = 1.;

      const auto& splines = tria.get_splines(); 

      int index = -1;
      for (const auto& t : splines) {
        const auto& A0 = t -> get_anchor(0);
        const auto& A1 = t -> get_anchor(1);
        
        if (A0.first == target_x && 
              A1.first == target_y_min && 
              A1.second == target_y_max)
          index = t -> get_level_index();
      }

      std::cout << "index = " << index << std::endl;
    }
  
    if (tria.n_levels() - 1 < 11) {
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

      const auto& bounding_box = tria.get_bounding_box();
      const auto& kv = data.kv;
      const int nx = kv[0].size();
      const int ny = kv[1].size();
  
      const unsigned int N1 = 100; 
      const unsigned int N2 = 100;
      const double xmin = kv[0][0];
      const double ymin = kv[1][0]; 
      const double xmax = kv[0][nx-1];
      const double ymax = kv[1][ny-1];
      std::cout << "xmin: " << xmin << ", xmax: " << xmax << std::endl;
      std::cout << "ymin: " << ymin << ", ymax: " << ymax << std::endl;
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
    //*/

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

    tria.minimal_residual_error_estimate(
                  {degrees[0]*degrees[0] + 1,
                   degrees[1]*degrees[1] + 1},
                   current_solution,
                   cell_errors
                   );
    data_out.add_data_vector(cell_errors, "cell_errors");


    // Build patches
    data_out.build_patches(); 

    // Open a file
    std::ofstream vtu_out(name_vtg); 
    data_out.write_vtu(vtu_out);

    // Print the table to a file preemptively
    // output to .txt
    problem_out.write_table_text();
    problem_out.write_table_tex();

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


  void Minimal_Benchmark::print_numerical_solution(std::string addition)
  {
    // Number of points to evaluate per direction
    const unsigned int N = 25;

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


  void Minimal_Benchmark::run()
  {
    std::cout << "Running benchmark ... " << std::endl;

    setup_system();
    // Set initial solution u_0 to zero. //TODO: Brauche ich das oder ist schon default?
    for (unsigned int i = 0; i < current_solution.size(); ++i)
      current_solution[i] = 0.;

    // boundary condition on the first solution u_0
    impose_boundary_condition(this -> problem_shape);
    

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
        impose_boundary_condition(this -> problem_shape);
      }


      if (tria.n_levels() < 15) {
      //output_system();
      std::cout << "    Initial residual: " << compute_residual_for_steplength(0) << std::endl;

      unsigned int newton_iteration = 0;
      double last_norm = 0., current_norm = 1., norm_diff = 1.;
      do {
        // solve system with zero boundary condition
        assemble_system();
        //system_matrix.print_pattern(std::cout);
        solve_system();
        print_numerical_solution();
        //output_system();


        //std::cout << "    Residual of N_It " << newton_iteration+1 
        //          << ": " << compute_residual_for_steplength(0) << std::endl;
        current_norm = newton_update.l2_norm();
        std::cout << "L2 ||update_n+1|| in N_It " << newton_iteration+1 
                  << ": " << current_norm << std::endl;

        //norm_diff = std::fabs(current_norm - last_norm);
        //std::cout << "    discr. ||update_{n+1} - update_n|| in N_It " << newton_iteration+1 
        //          << ": " << norm_diff << std::endl;

        last_norm = current_norm;
        newton_iteration++;

      } // do ( ... )
      while (newton_iteration < 30
             //&& norm_diff > 1e-4 
             && current_norm > 1e-2
            );

      std::cout << "Outputting system after last newton iteration ..." << std::endl;
      output_system();

      }
      
      cycle++;
    }
    std::cout << " ... done!" << std::endl;
  } // run




  // =================================================================

  // =================================================================
  double Minimal_DC_circle::value(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    double out = 1.;//std::sin(2 * numbers::PI * (p[0] + p[1]));
    std::cout << "Point: " << p[0] << ", " << p[1] << std::endl;

    return out;
  } // value_DC_circle



  double Minimal_DC_catenoid::value(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    double out = 1.;    // Height of the inner boundary
    return out;
  } // value_DC_Catenoid


  double Minimal_DC_square::value(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    double out = 0.;
    double a = 0.5;
    //if (p[0] >= 0. && p[0] < 0.5)
    //  out = a*p[0];
    //else if(p[0] >= 0.5 && p[0] <= 1.)
    //  out = a*(1 - p[0]);
    std::cout << "Point: " << p[0] << ", " << p[1] << std::endl;
    if (p[0] >= 0. && p[0] < 0.25)
      out = a*p[0];
    else if (p[0] >= 0.25 && p[0] < 0.5)
      out = a*(0.5 - p[0]);
    else if(p[0] >= 0.5 && p[0] <= 0.75)
      out = a*(p[0] - 0.5);
    else if (p[0] > 0.75 && p[0] <= 1)
      out = a*(1 - p[0]);


    return out;
  } // value_DC_square

  double Minimal_SOL_catenoid::value(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    double out = 0.;

    return out;
  } // value_SOL_catenoid

  Tensor<1, 2> Minimal_SOL_catenoid::gradient(
    const Point<2>&     p, 
    const unsigned int /* component */
  ) const {
    Tensor<1, 2> out;

    return out;
  } // gradient_SOL_catenoid










  // =================================================================

  // =================================================================
  IPF_Data<2, 2> Minimal_Benchmark::get_IPF_data(
    ProblemShape shape
  ) {
    switch(shape){
      case ProblemShape::Halfcircle:
        return get_IPF_data_halfcircle();
      case ProblemShape::Tilted_Halfcircle:
        return get_IPF_data_tiltedhalfcircle();
      case ProblemShape::Annulus:
        return get_IPF_data_annulus();
      case ProblemShape::Square:
        return get_IPF_data_square();
      default:
        AssertThrow(false, ExcNotImplemented());
    }
  } // get_IPF_data

  IPF_Data<2, 2> Minimal_Benchmark::get_IPF_data_halfcircle(
  ) {
    std::vector< std::vector< double > > kv;
    std::vector< Point<2 + 1> > cps;
    std::vector< unsigned int > deg;

    // define the knot vectors
    kv = std::vector< std::vector< double > >(2);
    kv[0] = {-1, -1, -1, 0, 1, 1, 1};
    kv[1] = {0, 0, 1, 1};

    // define the control points vector
    // Nach NURBSbook
    cps = std::vector< Point<3> >(8);
    cps[0] = Point<2 + 1>(-1. , 0., 1. );
    cps[1] = Point<2 + 1>(-0.5 , 0.5, 0.5);
    cps[2] = Point<2 + 1>( 0.5 , 0.5, 0.5);
    cps[3] = Point<2 + 1>( 1. , 0., 1. );
    cps[4] = Point<2 + 1>( 0. , 0., 1. );
    cps[5] = Point<2 + 1>(0. , 0., 0.5);
    cps[6] = Point<2 + 1>( 0. , 0., 0.5);
    cps[7] = Point<2 + 1>( 0. , 0., 1. );

    // define the degrees
    deg = {2, 1};

    return IPF_Data<2, 2>(cps, kv, deg);

  } // get_IPF_data_halfcircle

  IPF_Data<2, 2> Minimal_Benchmark::get_IPF_data_tiltedhalfcircle(
  ) {
    std::vector< std::vector< double > > kv;
    std::vector< Point<2 + 1> > cps;
    std::vector< unsigned int > deg;

    // define the knot vectors
    kv = std::vector< std::vector< double > >(2);
    kv[0] = {-1, -1, -1, 0, 1, 1, 1};
    kv[1] = {0, 0, 1, 1};

    // define the control points vector
    // Nach NURBSbook
    cps = std::vector< Point<3> >(8);
    double a = 1./std::sqrt(2); 
    cps[0] = Point<2 + 1>(-a ,-a , 1. );
    cps[1] = Point<2 + 1>(-a , 0., 0.5);
    cps[2] = Point<2 + 1>( 0., a , 0.5);
    cps[3] = Point<2 + 1>( a , a , 1. );
    cps[4] = Point<2 + 1>( 0., 0., 1. );
    cps[5] = Point<2 + 1>( 0., 0., 0.5);
    cps[6] = Point<2 + 1>( 0., 0., 0.5);
    cps[7] = Point<2 + 1>( 0., 0., 1. );



    // define the degrees
    deg = {2, 1};

    return IPF_Data<2, 2>(cps, kv, deg);

  } // get_IPF_data_tiltedhalfcircle

  IPF_Data<2, 2> Minimal_Benchmark::get_IPF_data_annulus(
  ) { 

    std::vector< std::vector< double > > kv;
    std::vector< Point<2 + 1> > cps;
    std::vector< unsigned int > deg;

    // define the knot vectors
    kv = std::vector< std::vector< double > >(2);
    //kv[0] = {-1, -1, -1 , -0.5, 0, 0, 0.5, 1, 1, 1};
    kv[0] = {-1, -1, -1 , 0, 1, 1, 1};
    kv[1] = {-1, -1, 1, 1};

    // define the control points vector
    cps = std::vector< Point<3> >(14);
    cps = std::vector< Point<3> >(8);
    unsigned int ind = 0; 
    const double r = 0.5;   // An annulus with outer radius 1 and inner radius 0 < r =< 1
    cps[ind++]  = Point<2 + 1>(-1. ,  0. , 1. );
    cps[ind++]  = Point<2 + 1>(-0.5,  0.5, 0.5);
    cps[ind++]  = Point<2 + 1>( 0.5,  0.5, 0.5);
    cps[ind++]  = Point<2 + 1>( 1. ,  0. , 1. );
    //cps[4]  = Point<2 + 1>( 0.5, -0.5, 0.5);
    //cps[5]  = Point<2 + 1>(-0.5, -0.5, 0.5);
    //cps[6]  = Point<2 + 1>(-1. ,  0. , 1. );

    cps[ind++]  = Point<2 + 1>(-1.*r ,  0.   , 1. );
    cps[ind++]  = Point<2 + 1>(-0.5*r,  0.5*r, 0.5);
    cps[ind++]  = Point<2 + 1>( 0.5*r,  0.5*r, 0.5);
    cps[ind++]  = Point<2 + 1>( 1.*r ,  0.   , 1. );
    //cps[11] = Point<2 + 1>( 0.5*r, -0.5*r, 0.5);
    //cps[12] = Point<2 + 1>(-0.5*r, -0.5*r, 0.5);
    //cps[13] = Point<2 + 1>(-1.*r ,  0.   , 1. );

    
    
    // define the degrees
    deg = {2, 1};

    return IPF_Data<2, 2>(cps, kv, deg);
  } // get_IPF_data_annulus

  IPF_Data<2, 2> Minimal_Benchmark::get_IPF_data_square(
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

    return IPF_Data<2, 2>(cps, kv, deg);
  } // get_IPF_data_square



}