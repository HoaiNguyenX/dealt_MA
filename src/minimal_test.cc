// #include <minimal_test.h>


// namespace Minimal_Surface {

//   using namespace dealt;
//   using namespace dealii;


//   Minimal_Test::Minimal_Test(
//     int ref,
//     int order,
//     const RefinementStrategy strategy
//     ) : ref(ref)
//       , order(order)
//       , data(this->get_IPF_data())
//       , tria(data)
//       , cycle(0)
//       , refinement_strategy(strategy)
//     {    
//       // Assembling the file name and elevating the degree
//       std::string name = "minimal_test";

//       if (strategy == RefinementStrategy::Adaptive)
//         name += "_adaptive";
//       else
//         name += "_uniform";
      
//       tria.degree_elevate_global(order);
//       tria.refine_bezier_elements();
//       const auto& cell = tria.begin_active();
//       const auto& CPs = tria.get_control_points(cell);
//       tria.coarsen_bezier_elements();

//       tria.prepare_assembly();

//       const std::vector<std::string>& columns =
//         {"Level", "Cycles", "Cells", "DoFs", "k_{Newton}", "update_norm", "residual"};
//       const std::vector<std::string>& tex_captions = 
//         {"Level", "Cycles", "\\# cells", "\\# dofs", "k_{Newton}", "update_norm", "residual"};
//       const std::vector<bool>& scientific =
//         {false, false, true, true, true, true, true,};
//       const std::vector<unsigned int> precision = 
//         {0, 0, 0, 0, 8, 8, 8};
//       const std::vector<std::string>& super_column_names = 
//         {"Grid Info", "Newton"};
//       const std::vector<std::vector<std::string>> super_columns =
//         {{"Level", "Cycles", "Cells", "DoFs"}, {"k_{Newton}", "update_norm", "residual"}};
//       problem_out = OutputSetup(  name
//                                 , data.max_degree() + order
//                                 , columns
//                                 , tex_captions
//                                 , scientific
//                                 , precision
//                                 , super_column_names
//                                 , super_columns
//                                 ); 
      



//       // Set boundary indicators
//       for (auto& face : tria.active_face_iterators())
//       { 
//         std::cout << "face: " << face -> index() << " with coords ";
//         std::cout << face->center()[0] << ", "<< face->center()[1]<<std::endl;
//         if (!face -> at_boundary())
//           continue;
//         const Point<2>& c = face -> center();
//         if (std::fabs(c(1)) < 1e-15 || std::fabs(c(1)-1) < 1e-15)
//           face -> set_boundary_id(Boundary::Dirichlet_1);
//         if (std::fabs(c(0)) < 1e-15 || std::fabs(c(0)-1) < 1e-15)
//           face -> set_boundary_id(Boundary::Dirichlet_0);
//       } // for ( face )  
//     } // Minimal_Benchmark


// void Minimal_Test::setup_system()
//   {

//     std::cout << "Setting up system ... " << std::endl;
//     unsigned int n_global_dofs = tria.n_active_splines();
//     system_rhs.reinit(n_global_dofs);
//     current_solution.reinit(n_global_dofs);
//     newton_update.reinit(n_global_dofs);
  
//     const auto& IEN_array = tria.get_IEN_array();
//     Assert(IEN_array.size() > 0, ExcInternalError());

//     sparsity_pattern.reinit(
//         n_global_dofs,
//         n_global_dofs,
//         n_global_dofs );
  
//     for (const auto& [_, arr] : IEN_array)
//       for (unsigned int i : arr)
//         for (unsigned int j : arr)
//           sparsity_pattern.add(i, j);
  
//     // Free superfluous space
//     sparsity_pattern.compress();
//     system_matrix.reinit(sparsity_pattern);
//   } // setup_system


// void Minimal_Test::assemble_system()
//   {
//     //std::cout << "Assembling system matrix ... " << std::endl;
  
//     // Setup initial tables that store the bernstein values / grads / hessians.
  
//     //reinit the system matrix and rhs since setpup is not called in newton loop
//     system_matrix = 0;
//     system_rhs    = 0;

//     std::vector< unsigned int > degrees = tria.get_degree();
//     degrees[0] = degrees[0]  + 1;
//     degrees[1] = degrees[1]  + 1;

//     TSValues<2> ts_values(
//         &tria,
//         degrees,
//         update_values |
//         update_gradients |
//         update_quadrature_points |
//         update_JxW_values);

//     const unsigned int dofs_per_cell = ts_values.n_dofs_per_cell();
//     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
//     Vector<double>     cell_rhs(dofs_per_cell);

//     const double power = 1.;

//     for (const auto& cell : tria.active_cell_iterators())
//     {
//       // Get the Bernstein values on the cell
//       ts_values.reinit(cell);

//       // Get the old_solution_gradient for the current cell
//       std::vector< Tensor<1, 2> > old_gradient(ts_values.n_quadrature_points_per_cell());
//       ts_values.get_function_gradients(current_solution, old_gradient);


//       // Reset the cell matrix
//       cell_matrix       = 0;
//       cell_rhs          = 0;

//       // Quadrature sum:
//       for (const unsigned int q : ts_values.quadrature_point_indices())
//       {
//         // coefficient a_n depending on old solution gradient for readability
//         const double coeff = 1. / std::pow(1 + old_gradient[q] * old_gradient[q],power);
//         const double coeff_diff = std::pow(1 + old_gradient[q] * old_gradient[q],power+1);
//         //const double coeff = 1. / 1 + old_gradient[q] * old_gradient[q];

//         // Build the cell matrix and rhs
//         for (const unsigned int i : ts_values.dof_indices())
//         {
          
//           for(const unsigned int j : ts_values.dof_indices())
//           {
//             cell_matrix(i,j) +=
//                     (((ts_values.shape_grad(i, q)      // ((\nabla \phi_i
//                        * coeff                         //   * a_n
//                        * ts_values.shape_grad(j, q))   //   * \nabla \phi_j)
//                       -                                //  -
//                       (2 * power
//                        * ts_values.shape_grad(i, q)    //  (\nabla \phi_i
//                        * coeff_diff                    //   * a_n^3
//                        * (ts_values.shape_grad(j, q)   //   * (\nabla \phi_j
//                             * old_gradient[q])         //       * \nabla u_n)
//                          * old_gradient[q]))           //     * \nabla u_n))
//                      * ts_values.JxW(q));              // * dx                         
//           } // for ( j )

//           cell_rhs(i) -= ((ts_values.shape_grad(i, q)             // \nabla \phi_i
//                               * coeff                             //    * a_n
//                               * old_gradient[q])                  //    * \nabla u_n
//                             * ts_values.JxW(q));                  // * dx
//         } // for ( i )
//       } // for ( q )

      
//       // Add the values to the system
//       std::vector< unsigned int > local_dof_indices =
//               tria.get_IEN_array(cell);
      

//       system_matrix.add(local_dof_indices, cell_matrix);
//       system_rhs.add(local_dof_indices, cell_rhs);
//     } // for ( cell )

//     unsigned int n_global_dofs = tria.n_active_splines();
//     const auto& boundary_dofs = tria.get_boundary_dofs();

//     // Impose zero Dirichlet boundary condition to the system for
//     // all boundary dofs
//     for (const auto& [boundary, dofs] : boundary_dofs){
//       for (const auto& dof : dofs){

//         for (unsigned int i = 0; i < n_global_dofs; i++){
//           system_matrix.set(i, dof, 0.);
//           system_matrix.set(dof, i, 0.);
//         }
//         system_matrix.set(dof, dof, 1.);
//         system_rhs(dof) =  0.;
//         //std::cout << "setting boundary at dof: " << dof << std::endl;
//       }
//     }
//   } // assemble_system


// void Minimal_Test::impose_boundary_condition()
//   {
    
//     std::cout << "Imposing boundary condition ..." << std::endl;
//     const auto& boundary_dofs = tria.get_boundary_dofs();
//     unsigned int n_global_dofs = tria.n_active_splines();
//     // For the case the initial mesh has no boundary dofs
//     if (boundary_dofs.size() == 0)
//       return;
//     // print number of dofs at boundary
//     //std::cout << "boundary_dofs: " << std::endl;
//     //for (const auto& [boundary, dofs] : boundary_dofs)
//     //  std::cout << "Boundary " << boundary << ": " << dofs.size() << std::endl;
//     //std::cout << "number of global dofs: " << tria.n_active_splines() << std::endl;

//     std::vector< unsigned int > degrees = tria.get_degree();
//     degrees[0] = degrees[0]  + 1;
//     degrees[1] = degrees[1]  + 1;

//     std::map<
//         types::global_dof_index, 
//         double
//       > boundary_values;


//     // Firstly set the nonzero Dirichlet boundary values
//     if (boundary_dofs.find(Boundary::Dirichlet_1) 
//         == boundary_dofs.end())
//       return;
//    std::map<
//      types::boundary_id,
//      const Function< 2 >*
//    > boundary_fcns =
//        {{Boundary::Dirichlet_1, &DC_fcn}};
//    tria.project_boundary_values(
//          boundary_fcns,
//          degrees,
//          boundary_values);


//     // Secondly set the zero Dirichlet boundary values
//     if (boundary_dofs.find(Boundary::Dirichlet_0) 
//           != boundary_dofs.end())
//       for (const auto& dof : boundary_dofs.at(Boundary::Dirichlet_0))
//         boundary_values[dof] = 0.;


    

//     //TODO: Für circle sind die dofs noch Null obwohl Dirichlet_c gesetzt ist
//     // Die Punkte, die in DC_c inputted werden, sind irgendwie alle Null???
//     // Obwohl .vtg und matlab das richtige zeigt
//     // print boundary values
//     //Boundary bndry = Boundary::Dirichlet_1;
//     //if (boundary_dofs.find(Boundary::Dirichlet_1) 
//     //      != boundary_dofs.end())
//     //{
//     //  std::cout << "Boundary dofs of Dirichlet_1: " << std::endl;
//     //  const auto& splines = tria.get_splines(); 
//     //  for (const auto& dof : boundary_dofs.at(bndry)){
//     //    const auto& ts = splines.at(dof); 
//     //    const auto& anchor = ts -> get_anchor();
//     //    std::cout << dof << ": " << 0.5 * anchor.first + 0.5*anchor.second << ", value = " << boundary_values.at(dof) << std::endl;
//     //  } // for ( dof )
//     //}
    
//     //// set Dirichlet boundary values to current_solution
//     typename std::map<types::global_dof_index, double>::const_iterator
//       dof  = boundary_values.begin(),
//       endd = boundary_values.end();
//     for (; dof != endd; ++dof)
//       {
//         Assert(dof->first < n_global_dofs, ExcInternalError());
//         current_solution(dof->first) = dof->second;
//       }


//     //MatrixTools::apply_boundary_values(
//     //  boundary_values,
//     //  system_matrix,
//     //  current_solution,
//     //  system_rhs
//     //);
    
//   } // impose_boundary_condition


//   void Minimal_Test::solve_system_constant()
//   {
//     // std::cout << "Solving system ... " << std::endl;

//     //SolverControl            solver_control(1e3,
//     //                             1e-5);
//     //SolverCG<Vector<double>> solver(solver_control);
// //
//     //PreconditionJacobi<SparseMatrix<double>> preconditioner;
//     //preconditioner.initialize(system_matrix);
// //
//     //solver.solve(system_matrix, newton_update, system_rhs, preconditioner);

//     // Solve with direct UMFPACK
//     newton_update = system_rhs;
//     SparseDirectUMFPACK A_direct;
//     A_direct.solve(system_matrix, newton_update);

//     // print newton_update
//     //std::cout << "newton_update: " << std::endl;
//     //for (unsigned int i = 0; i < tria.n_active_splines(); i++)
//     //  std::cout << newton_update[i] << " ";
//     //std::cout << std::endl;

//     // Add the update to the current solution
//     double alpha = determine_step_length_const();
//     current_solution.add(alpha, newton_update);

//     // Add current solution values to the Tsplines data structure
//     const auto& splines = tria.get_splines();
//     for (unsigned int i = 0; i < tria.n_active_splines(); i++)
//         splines[i] -> set_solution(current_solution[i]);
//   } // solve_system


//   double Minimal_Test::compute_residual_for_steplength(double alpha)
//   {
//     // Compute the discrete residual for the current solution
//     double n_global_dofs = tria.n_active_splines();
//     Vector<double> evaluation_points(n_global_dofs);
//     evaluation_points = current_solution;
//     evaluation_points.add(alpha, newton_update);

//     Vector<double> cell_residuals(n_global_dofs);
//     const std::vector<unsigned int>& degrees = tria.get_degree();

//     TSValues<2> ts_values(
//         &tria,
//         degrees,
//         update_values |
//         update_gradients |
//         update_quadrature_points |
//         update_JxW_values);

//     const unsigned int nvc = GeometryInfo<2>::vertices_per_cell;
//     //const unsigned int nvf = GeometryInfo<dim>::vertices_per_face;
    
    
//     // Loop over all cells
//     for (const auto& cell : tria.active_cell_iterators()) 
//     {

//       std::vector< unsigned int > local_dof_indices = tria.get_IEN_array(cell);
//       ts_values.reinit(cell);

//       // get old solution gradient
//       std::vector< Tensor<1, 2> > sol_grad(ts_values.n_quadrature_points_per_cell());
//       ts_values.get_function_gradients(evaluation_points, sol_grad);

//       // Quadrature sum:
//       for (const unsigned int q : ts_values.quadrature_point_indices())
//       {
//         // coefficient a_n depending on old solution gradient for visibilities
//         const double coeff = 1. / std::sqrt(1 + sol_grad[q] * sol_grad[q]);

//         // Fehler aber keine Ahnung warum bei ts_values.dof_indices()
//         for (const unsigned int i : ts_values.dof_indices())
//         {

//           // root and squares cancel out -> res_cell = cell_width^2 * norm^2
//           cell_residuals(i) += (ts_values.shape_grad(i, q)        // \nabla \phi_i
//                                * coeff                        // * a_n
//                                * sol_grad[q])                 // * \nabla u_n
//                                * ts_values.JxW(q);            // * dx
        
//         } // for ( i )

//       } // for ( q )

//     } // for (cell)

//     //zero_constrains.distribute_local_to_global(cell_residuals, 
//     //                                          local_dof_indices,residual_norm);
//     // Set zero Dirichlet boundary condition to the residual according to deal.II step-15
//     const auto& boundary_dofs = tria.get_boundary_dofs();
//     if (boundary_dofs.find(Boundary::Dirichlet_0) 
//           != boundary_dofs.end())
//       for (const auto& dof : boundary_dofs.at(Boundary::Dirichlet_0))
//         cell_residuals[dof] = 0.;


//     // Compute the L2 norm of the residual
//     double residual_norm = 0.;
//     for (unsigned int i = 0; i < n_global_dofs; ++i)
//       residual_norm += cell_residuals(i) * cell_residuals(i);

//     return std::sqrt(residual_norm);
//   } // compute_residual_for_steplength



//   void Minimal_Test::estimate_and_mark()
//   {
//     std::cout << "    Refining grid..." << std::endl;
//     if (refinement_strategy == RefinementStrategy::Uniform) 
//     {
//       tria.coarsen_bezier_elements();
//       tria.refine_global();   
//     } 
//     else // RefinementStrategy::Adaptive
//     {
//       std::vector< TriaIterator<CellAccessor<2, 2>> > mark;
//       const std::vector<unsigned int>& degrees = tria.get_degree();
//       Vector<double>  local_residuals(tria.n_active_cells());

//       tria.minimal_residual_error_estimate(
//                           {degrees[0]*degrees[0] + 1,
//                            degrees[1]*degrees[1] + 1},
//                            current_solution,
//                            local_residuals
//                            );


//       tria.refine_fixed_number(local_residuals, 0.10);
//     }

//     tria.prepare_assembly();
//   } // estimate_and_mark



//   void Minimal_Test::print_numerical_solution(std::string addition)
//   {
//     // Number of points to evaluate per direction
//     const unsigned int N = 25;

//     // Declare a container to store the values 
//     FullMatrix<double> B_num(3,N*N);

//     const auto& splines = tria.get_splines();
//     const auto& mapping = tria.get_IPF();
//     const auto& kv      = data.kv;
//     const unsigned int n0 = kv[0].size() - 1;
//     const unsigned int n1 = kv[1].size() - 1;   

    
//     for (unsigned int i = 0; i < N; i++)
//     {
//       for (unsigned int j = 0; j < N; j++)
//       {
//         const Point<2> P(kv[0][0] + (kv[0][n0] -  kv[0][0]) * i / (N-1.), 
//                           kv[1][0] + (kv[1][n1] -  kv[1][0]) * j / (N-1.));
//         B_num(0, i*N + j) = P[0];
//         B_num(1, i*N + j) = P[1];
//         for (unsigned int k = 0; k < tria.n_active_splines(); k++)
//         {
//           B_num(2, i*N + j) += current_solution[k] * splines[k] -> value(P);
//         }
//       }
//     }

    
    
//     std::string name_num =  problem_out.dat.string() 
//                             + "numsol_" + addition + "l" 
//                             + std::to_string(tria.n_levels() - 1) 
//                             + ".dat";
//     //std::cout << name_num << std::endl;
//     std::filebuf f0;
//     f0.open(name_num.c_str(), std::ios::out);

//     std::ostream out0(&f0);
//     B_num.print_formatted(out0, 16, true, 1, "0");

//     f0.close();


//     GridOutFlags::Svg svg_flags;
//     svg_flags.coloring = GridOutFlags::Svg::Coloring::none;
//     // svg_flags.label_level_number  = true;
//     // svg_flags.label_cell_index    = true;
//     svg_flags.label_boundary_id   = true;
    
//     std::string svg_name = problem_out.svg.string() + "parametric_grid_l"
//                           + std::to_string(tria.n_levels() - 1)
//                           + ".svg";  
//     //std::cout << "... printing to: " << svg_name << std::endl;
//     std::ofstream out(svg_name);
//     GridOut       grid_out;
//     grid_out.set_flags(svg_flags);
//     grid_out.write_svg(tria, out);
//     out.close();
//   } // print_numerical_solution


//   void Minimal_Test::run()
//   {
//     std::cout << "Running benchmark ... " << std::endl;
//     Vector<double> residuals(tria.n_active_splines());
//     setup_system();
//     // Set initial solution u_0 to zero. //TODO: Brauche ich das oder ist schon default?
//     for (unsigned int i = 0; i < current_solution.size(); ++i)
//       current_solution[i] = 0.;

//     // boundary condition on the first solution u_0
//     impose_boundary_condition();
    

//     // Refinement cycle
//     unsigned int cycle = 0;
//     while (tria.n_levels() < this -> ref + 1)
//     {
//       std::cout << "  Refinement cycle " << cycle << ':' << std::endl;

//       if (cycle != 0)
//       {
//         estimate_and_mark();
//         setup_system();
//         residuals.reinit(tria.n_active_splines());
//         tria.transfer_solution(current_solution);
//         impose_boundary_condition();
//       }

//       std::cout << " n_levels: " << tria.n_levels() << std::endl;
//       std::cout << "    Initial residual: " << compute_residual_for_steplength(0) << std::endl;

//       unsigned int newton_iteration = 0;
//       double last_norm = 0., current_norm = 1., norm_diff = 1., last_residual = 0.;
//       do {
//         // solve system with zero boundary condition
//         assemble_system();

//         solve_system_constant();
//         residuals = system_rhs;


//         double current_residual = residuals.l2_norm();
//         std::cout << "Residual of N_It " << newton_iteration+1 
//                   << ":   " << std::fixed << std::setprecision(6) << current_residual;

//         current_norm = newton_update.l2_norm();
//         std::cout << "      L2 ||update_n+1|| " 
//                   << ":   " << std::fixed << std::setprecision(6) << current_norm << std::endl;


//         last_norm = current_norm;
//         last_residual = current_residual;
//         newton_iteration++;

//       } // do ( ... )
//       while (newton_iteration < 1000
//              //&& norm_diff > 1e-4 
//              && current_norm > 5e-8
//             );
//       std::cout << "    Last residual: " << last_residual << std::endl;

//       std::cout << "Outputting system after last newton iteration ..." << std::endl;
//       if (cycle != 0)
//       {
//         print_numerical_solution();
//         problem_out.add_values_to_table(
//           tria.n_levels() - 1,
//           cycle,
//           tria.n_active_cells(),
//           tria.n_active_splines(),
//           newton_iteration,
//           current_norm,
//           last_residual
//         );
//       }
//       cycle++;
//     }
//     std::cout << " ... done!" << std::endl;
//   } // run




//   double Minimal_DC::value(
//     const Point<2>&     p, 
//     const unsigned int /* component */
//   ) const {
//     double out = 0.;
//     double a = 0.5;
//     //if (p[0] >= 0. && p[0] < 0.5)
//     //  out = a*p[0];
//     //else if(p[0] >= 0.5 && p[0] <= 1.)
//     //  out = a*(1 - p[0]);
//     if (p[0] >= 0. && p[0] < 0.25)
//       out = a*p[0];
//     else if (p[0] >= 0.25 && p[0] < 0.5)
//       out = a*(0.5 - p[0]);
//     else if(p[0] >= 0.5 && p[0] <= 0.75)
//       out = a*(p[0] - 0.5);
//     else if (p[0] > 0.75 && p[0] <= 1)
//       out = a*(1 - p[0]);
//     return out;
//   } // value_DC

//   IPF_Data<2, 2> Minimal_Test::get_IPF_data(
//   ) {
//     std::vector< std::vector< double > > kv;
//     std::vector< Point<2 + 1> > cps;
//     std::vector< unsigned int > deg;

//     // define the knot vectors
//     kv = std::vector< std::vector< double > >(2);
//     kv[0] = {0, 0, 1, 1};
//     kv[1] = {0, 0, 1, 1};


//     // define the control points vector
//     cps = std::vector< Point<3> >(4);
//     cps[0] = Point<2 + 1>( 0.,  0.,  1.);
//     cps[1] = Point<2 + 1>( 1.,  0.,  1.);
//     cps[2] = Point<2 + 1>( 0.,  1.,  1.);
//     cps[3] = Point<2 + 1>( 1.,  1.,  1.);

//     // define the degrees
//     deg = {1, 1};

//     return IPF_Data<2, 2>(cps, kv, deg);
//   } // get_IPF_data_square


// } // namespace Minimal_Surface


